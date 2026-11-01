//! Bodenkontakt.
//!
//! Der Boden ist eine Zwangsbedingung, keine Feder: Kontaktpunkte dürfen
//! sich nicht in die Oberfläche hineinbewegen. Pro Schritt wird der kleinste
//! Stoß gesucht, der das erfüllt (Normalstoß ≥ 0, Reibung nach Coulomb),
//! iterativ über alle Kontaktpunkte (Projected Gauss-Seidel). Ein Schiff,
//! das steht, erfährt so genau die Stützkraft, die es trägt; ein Schiff mit
//! genug Schub hebt ab, weil die Stützkraft nicht negativ werden kann.
//!
//! Die Stützkraft geht als äußere Kraft in die Lastberechnung ein. Damit ein
//! Aufprall, der in einem einzigen Schritt abgebremst wird, nicht eine von der
//! Schrittweite abhängige Kraftspitze erzeugt, wird die Kontaktkraft über die
//! Nachgiebigkeit der Struktur [`IMPACT_TIME`] geglättet: Ruhende Last kommt
//! voll an, ein Stoß verteilt sich über diese Zeit.
//!
//! Der Boden ist (noch) die glatte Kugel mit dem mittleren Radius und dreht
//! mit dem Planeten mit.

use crate::loads::ExternalLoad;
use crate::mass::MassProperties;
use crate::vessel::Vessel;
use crate::dynamics::RigidState;
use glam::{DMat3, DVec3};
use sim_core::{CelestialBody, Vec3, VesselFrame};

/// Abstand, ab dem ein Punkt als Kontakt zählt (m).
pub const CONTACT_SLOP: f64 = 0.01;
/// Reibungsbeiwert Teil/Boden.
pub const FRICTION: f64 = 0.8;
/// Zeit, über die sich ein Aufprall in der Struktur verteilt (s).
pub const IMPACT_TIME: f64 = 0.05;
const ITERATIONS: usize = 30;

/// Kontaktpunkte: die Eckpunkte des Polygonmodells aller Teile.
pub fn contact_points(v: &Vessel) -> Vec<(usize, Vec3<VesselFrame>)> {
    v.geometry().contact_points.clone()
}

/// Größter Abstand eines Kontaktpunkts vom Schwerpunkt.
pub fn extent(v: &Vessel, mp: &MassProperties) -> f64 {
    v.geometry().extent_from(mp.com)
}

#[derive(Debug, Clone, Default)]
pub struct ContactResult {
    /// Kontaktkräfte für die Lastberechnung (geglättet), Schiffssystem.
    pub loads: Vec<ExternalLoad>,
    /// Teil, dessen Aufprallgrenze überschritten wurde, mit Aufprallgeschwindigkeit.
    pub crash: Option<(usize, f64)>,
    /// Größte Annäherungsgeschwindigkeit an den Boden in diesem Schritt.
    pub impact_speed: f64,
    /// Mindestens ein Punkt berührt den Boden.
    pub touching: bool,
}

struct Contact {
    part: usize,
    /// Hebelarm vom Schwerpunkt (inertial).
    r: DVec3,
    /// Bodennormale (inertial).
    n: DVec3,
    depth: f64,
    lambda_n: f64,
    lambda_t: DVec3,
}

/// Löst den Bodenkontakt für den Zustand am Ende eines Schritts.
pub fn resolve(v: &Vessel, mp: &MassProperties, s: &mut RigidState, body: &CelestialBody, dt: f64) -> ContactResult {
    let rot: DMat3 = s.attitude.matrix();
    let pos = s.pos.raw();
    let com = mp.com.raw();
    let geo = v.geometry();
    let mut contacts: Vec<Contact> = geo
        .contact_points
        .iter()
        .filter_map(|&(part, p)| {
            let r = rot * (p.raw() - com);
            let world = pos + r;
            let dist = world.length();
            let alt = dist - body.radius;
            (alt < CONTACT_SLOP).then(|| Contact { part, r, n: world / dist, depth: -alt, lambda_n: 0.0, lambda_t: DVec3::ZERO })
        })
        .collect();

    if contacts.is_empty() {
        s.contact_force = 0.0;
        return ContactResult::default();
    }

    let m_inv = 1.0 / mp.mass;
    let i_inv = rot * mp.inertia.0.inverse() * rot.transpose();
    let omega_p = body.angular_velocity().raw();
    let mut vel = s.vel.raw();
    let mut w = rot * s.omega.raw();
    let rel = |vel: DVec3, w: DVec3, c: &Contact| vel + w.cross(c.r) - omega_p.cross(pos + c.r);
    let eff_mass = |c: &Contact, d: DVec3| 1.0 / (m_inv + d.dot((i_inv * c.r.cross(d)).cross(c.r)));

    // Aufprall vor dem Abbremsen
    let mut result = ContactResult { touching: true, ..Default::default() };
    for c in &contacts {
        let approach = -rel(vel, w, c).dot(c.n);
        if approach > result.impact_speed {
            result.impact_speed = approach;
        }
        let tolerance = v.parts[c.part].def.crash_tolerance;
        if approach > tolerance && result.crash.is_none_or(|(_, a)| approach > a) {
            result.crash = Some((c.part, approach));
        }
    }

    // 1) Haften: Kann der Boden das Schiff ganz zur Ruhe (relativ zur
    //    Oberfläche) bringen? Exakt und ohne Iteration; deckt Stehen und
    //    jedes Aufsetzen ohne Rutschen/Kippen ab.
    // 2) Sonst (Rutschen, Kippen, Abheben): iteratives Verfahren.
    let impulses = match stick(&contacts, mp, rot, i_inv, pos, vel, w, omega_p) {
        Some((v_new, w_new, imp)) => {
            vel = v_new;
            w = w_new;
            imp
        }
        None => {
            pgs(&mut contacts, m_inv, i_inv, &mut vel, &mut w, &rel, &eff_mass);
            contacts.iter().map(|c| c.n * c.lambda_n + c.lambda_t).collect()
        }
    };

    // Eindringen beheben (reine Lagekorrektur, ohne Energie zuzuführen)
    if let Some(deepest) = contacts.iter().max_by(|a, b| a.depth.total_cmp(&b.depth))
        && deepest.depth > 0.0
    {
        s.pos = Vec3::from_raw(pos + deepest.n * deepest.depth);
    }
    s.vel = Vec3::from_raw(vel);
    s.omega = Vec3::from_raw(rot.transpose() * w);

    // Kontaktkraft, über die Nachgiebigkeit geglättet
    let total_now: f64 = impulses.iter().map(|j| j.length()).sum::<f64>() / dt;
    let a = (dt / IMPACT_TIME).min(1.0);
    s.contact_force += (total_now - s.contact_force) * a;
    let scale = if total_now > 0.0 { s.contact_force / total_now } else { 0.0 };
    let to_body = rot.transpose();
    result.loads = contacts
        .iter()
        .zip(&impulses)
        .filter(|(_, j)| j.length_squared() > 0.0)
        .map(|(c, j)| ExternalLoad {
            part: c.part,
            force: Vec3::from_raw(to_body * (*j / dt * scale)),
            point: Vec3::from_raw(com + to_body * c.r),
            torque: Vec3::ZERO,
        })
        .collect();
    result
}

/// Haftlösung: Stoß, der das Schiff relativ zum Boden vollständig anhält.
/// Zulässig, wenn der Boden ihn liefern kann: Normalstoß > 0, Reibung im
/// Reibkegel, Druckmittelpunkt innerhalb der Standfläche. Rückgabe: neue
/// Geschwindigkeit, Winkelgeschwindigkeit und Stoß pro Kontaktpunkt.
#[allow(clippy::too_many_arguments)]
fn stick(
    contacts: &[Contact],
    mp: &MassProperties,
    rot: DMat3,
    i_inv: DMat3,
    pos: DVec3,
    vel: DVec3,
    w: DVec3,
    omega_p: DVec3,
) -> Option<(DVec3, DVec3, Vec<DVec3>)> {
    if contacts.len() < 3 {
        return None;
    }
    let inertia = rot * mp.inertia.0 * rot.transpose();
    let v_target = omega_p.cross(pos);
    let j = (v_target - vel) * mp.mass;
    let l = inertia * (omega_p - w);
    let _ = i_inv;

    let n = contacts.iter().map(|c| c.n).sum::<DVec3>().normalize();
    let jn = j.dot(n);
    if jn <= 0.0 {
        return None; // das Schiff drückt nicht auf den Boden
    }
    let jt = j - n * jn;

    // Ebene Koordinaten der Kontaktpunkte um ihren Mittelpunkt
    let c = contacts.iter().map(|c| c.r).sum::<DVec3>() / contacts.len() as f64;
    let e1 = n.any_orthonormal_vector();
    let e2 = n.cross(e1);
    let d: Vec<(f64, f64)> = contacts.iter().map(|k| ((k.r - c).dot(e1), (k.r - c).dot(e2))).collect();

    // Moment um den Mittelpunkt → Druckmittelpunkt und Torsion
    let m = l - c.cross(j);
    let m_n = m.dot(n);
    let cp = n.cross(m - n * m_n) / jn;
    let cp = (cp.dot(e1), cp.dot(e2));
    let rho = d.iter().map(|(x, y)| x.hypot(*y)).sum::<f64>() / d.len() as f64;
    if jt.length() + m_n.abs() / rho.max(1e-6) > FRICTION * jn {
        return None; // rutscht oder dreht sich auf dem Boden
    }
    if !inside_hull(&d, cp, 0.98) {
        return None; // kippt
    }

    // Stoß auf die Punkte verteilen (für die Lastberechnung): linear über die
    // Fläche, so dass der Druckmittelpunkt stimmt.
    let k = d.len() as f64;
    let (mx, my) = d.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0 / k, a.1 + p.1 / k));
    let (mut cxx, mut cxy, mut cyy) = (0.0, 0.0, 0.0);
    for (x, y) in &d {
        cxx += (x - mx) * (x - mx);
        cxy += (x - mx) * (y - my);
        cyy += (y - my) * (y - my);
    }
    let det = cxx * cyy - cxy * cxy;
    let (bx, by) = (cp.0 - mx, cp.1 - my);
    let (ax, ay) = if det.abs() > 1e-12 { ((cyy * bx - cxy * by) / det, (cxx * by - cxy * bx) / det) } else { (0.0, 0.0) };
    let weights: Vec<f64> = d.iter().map(|(x, y)| (1.0 / k + ax * (x - mx) + ay * (y - my)).max(0.0)).collect();
    let sum: f64 = weights.iter().sum();
    let imp = weights.iter().map(|wt| j * (wt / sum)).collect();

    Some((v_target, omega_p, imp))
}

/// Liegt `p` in der konvexen Hülle der Punkte (um `shrink` zum Schwerpunkt
/// der Hülle verkleinert)?
fn inside_hull(points: &[(f64, f64)], p: (f64, f64), shrink: f64) -> bool {
    let mut pts: Vec<(f64, f64)> = points.to_vec();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0);
    let mut hull: Vec<(f64, f64)> = Vec::new();
    for pass in 0..2 {
        let start = hull.len();
        let iter: Box<dyn Iterator<Item = &(f64, f64)>> = if pass == 0 { Box::new(pts.iter()) } else { Box::new(pts.iter().rev()) };
        for &q in iter {
            while hull.len() >= start + 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], q) <= 1e-12 {
                hull.pop();
            }
            hull.push(q);
        }
        hull.pop();
    }
    if hull.len() < 3 {
        return false;
    }
    let (gx, gy) = hull.iter().fold((0.0, 0.0), |a, q| (a.0 + q.0, a.1 + q.1));
    let g = (gx / hull.len() as f64, gy / hull.len() as f64);
    let hull: Vec<(f64, f64)> = hull.iter().map(|q| (g.0 + (q.0 - g.0) * shrink, g.1 + (q.1 - g.1) * shrink)).collect();
    (0..hull.len()).all(|i| cross(hull[i], hull[(i + 1) % hull.len()], p) >= 0.0)
}

/// Projected Gauss-Seidel: Normalstoß ≥ 0, Reibung ≤ μ·Normalstoß.
/// Abwechselnd vorwärts und rückwärts, damit die Reihenfolge der Punkte
/// keine einseitige Restdrehung hinterlässt.
fn pgs(
    contacts: &mut [Contact],
    m_inv: f64,
    i_inv: DMat3,
    vel: &mut DVec3,
    w: &mut DVec3,
    rel: &dyn Fn(DVec3, DVec3, &Contact) -> DVec3,
    eff_mass: &dyn Fn(&Contact, DVec3) -> f64,
) {
    let n = contacts.len();
    for it in 0..ITERATIONS {
        for k in 0..n {
            let c = &mut contacts[if it % 2 == 0 { k } else { n - 1 - k }];
            let vn = rel(*vel, *w, c).dot(c.n);
            let new = (c.lambda_n - vn * eff_mass(c, c.n)).max(0.0);
            let j = c.n * (new - c.lambda_n);
            c.lambda_n = new;
            *vel += j * m_inv;
            *w += i_inv * c.r.cross(j);

            let vr = rel(*vel, *w, c);
            let vt = vr - c.n * vr.dot(c.n);
            let speed = vt.length();
            if speed > 1e-12 {
                let t = vt / speed;
                let mut new_t = c.lambda_t - t * (speed * eff_mass(c, t));
                let limit = FRICTION * c.lambda_n;
                if new_t.length() > limit {
                    new_t = new_t.normalize_or_zero() * limit;
                }
                let jt = new_t - c.lambda_t;
                c.lambda_t = new_t;
                *vel += jt * m_inv;
                *w += i_inv * c.r.cross(jt);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::{Controls, DT_PHYSICS, step};
    use crate::mass::mass_properties;
    use crate::part::PartLibrary;
    use crate::scenario::{com_height, first_rocket, test_planet};
    use sim_core::GeoPosition;

    /// Rakete, Startplatz und Höhe der Unterkante über Grund.
    fn setup(drop_height: f64, sink_speed: f64) -> (Vessel, RigidState, CelestialBody) {
        let body = test_planet();
        let v = first_rocket(&PartLibrary::builtin());
        let mut s = RigidState::at_rest(&body, GeoPosition { latitude: 0.0, longitude: 0.0, altitude: com_height(&v) + drop_height }, 0.0);
        s.vel -= s.pos.normalize_or_zero() * sink_speed;
        (v, s, body)
    }

    fn run(v: &mut Vessel, s: &mut RigidState, body: &CelestialBody, c: &Controls, seconds: f64) -> Vec<crate::StepReport> {
        let n = (seconds / DT_PHYSICS).round() as usize;
        (0..n).map(|_| step(v, s, c, body, DT_PHYSICS)).collect()
    }

    fn bottom_altitude(v: &Vessel, s: &RigidState, body: &CelestialBody) -> f64 {
        s.altitude(body) - com_height(v)
    }

    #[test]
    fn rocket_rests_on_the_ground() {
        let (mut v, mut s, body) = setup(0.0, 0.0);
        let reports = run(&mut v, &mut s, &body, &Controls::default(), 20.0);
        let last = reports.last().unwrap();
        assert!(last.touching_ground && last.crash.is_none());
        assert!(bottom_altitude(&v, &s, &body).abs() < CONTACT_SLOP, "Unterkante {} m", bottom_altitude(&v, &s, &body));
        assert!(s.surface_velocity(&body).length() < 1e-3, "rutscht mit {:?}", s.surface_velocity(&body));

        // Stützkraft = Gewicht minus Zentrifugalanteil am Äquator
        let m = v.total_mass();
        let r = s.pos.length();
        let g_eff = body.mu / (r * r) - body.rotation_rate.powi(2) * r;
        assert!((last.ground_force - m * g_eff).abs() / (m * g_eff) < 0.01, "Stützkraft {} N, erwartet {} N", last.ground_force, m * g_eff);

        // Die Verbindung zum unteren Triebwerk trägt alles darüber
        let engine = v.parts.iter().rposition(|p| p.def.engine().is_some()).unwrap();
        let joint = last.loads.joints.iter().find(|j| j.child == engine).unwrap();
        let above = m - v.parts[engine].mass();
        assert!((-joint.axial - above * g_eff).abs() / (above * g_eff) < 0.02, "Druck {} N, erwartet {} N", -joint.axial, above * g_eff);
    }

    #[test]
    fn rocket_lifts_off_with_enough_thrust() {
        let (mut v, mut s, body) = setup(0.0, 0.0);
        v.activate_next_stage();
        let h0 = bottom_altitude(&v, &s, &body);
        run(&mut v, &mut s, &body, &Controls { throttle: 1.0, ..Default::default() }, 5.0);
        assert!(bottom_altitude(&v, &s, &body) - h0 > 10.0);
    }

    #[test]
    fn rocket_stays_down_with_too_little_thrust() {
        let (mut v, mut s, body) = setup(0.0, 0.0);
        v.activate_next_stage();
        // Schub/Gewicht ≈ 0,5
        run(&mut v, &mut s, &body, &Controls { throttle: 0.3, ..Default::default() }, 3.0);
        assert!(bottom_altitude(&v, &s, &body).abs() < CONTACT_SLOP);
    }

    #[test]
    fn gentle_landing_survives_and_settles() {
        let (mut v, mut s, body) = setup(0.5, 4.0);
        let reports = run(&mut v, &mut s, &body, &Controls::default(), 5.0);
        assert!(reports.iter().all(|r| r.crash.is_none()));
        let impact = reports.iter().map(|r| r.impact_speed).fold(0.0, f64::max);
        assert!(impact > 4.0 && impact < 8.0, "Aufprall {impact} m/s");
        assert!(s.surface_velocity(&body).length() < 1e-2);
        // Kein Abprallen: nach dem Aufsetzen bleibt die Unterkante unten
        let touched = reports.iter().position(|r| r.touching_ground).unwrap();
        assert!(reports[touched..].iter().all(|r| r.touching_ground));
    }

    #[test]
    fn hard_landing_crashes() {
        let (mut v, mut s, body) = setup(0.5, 20.0);
        let reports = run(&mut v, &mut s, &body, &Controls::default(), 1.0);
        let crash = reports.iter().find_map(|r| r.crash).expect("Aufprall mit 20 m/s muss zerstören");
        assert!(v.parts[crash.0].def.engine().is_some(), "das Triebwerk schlägt zuerst auf");
    }

    #[test]
    fn impact_load_does_not_depend_on_step_size() {
        // Gleicher Aufprall, halbe Schrittweite: Spitzenlast fast gleich
        let peak = |dt: f64| {
            let (mut v, mut s, body) = setup(0.2, 6.0);
            let n = (1.0 / dt) as usize;
            (0..n)
                .map(|_| step(&mut v, &mut s, &Controls::default(), &body, dt))
                .map(|r| r.loads.max_utilization().map_or(0.0, |j| j.utilization))
                .fold(0.0, f64::max)
        };
        let (a, b) = (peak(0.02), peak(0.01));
        assert!((a - b).abs() / a < 0.15, "Spitze {a} bei 20 ms, {b} bei 10 ms");
    }

    #[test]
    fn tilted_rocket_tips_over() {
        // 15° geneigt auf der Fußkante abgesetzt: Schwerpunkt außerhalb der
        // Standfläche, sie muss umfallen (Haftlösung darf nicht greifen).
        let (mut v, mut s, body) = setup(0.0, 0.0);
        let up = s.pos.normalize_or_zero();
        s.attitude = s.attitude * sim_core::Rotation::from_axis_angle(Vec3::X, 15f64.to_radians());
        // so anheben, dass die tiefste Fußkante gerade aufliegt
        let mp = mass_properties(&v);
        let low = contact_points(&v)
            .iter()
            .map(|(_, p)| (s.pos + s.attitude.apply(*p - mp.com)).length() - body.radius)
            .fold(f64::INFINITY, f64::min);
        s.pos -= up * low;
        let tilt = |s: &RigidState| s.nose().dot(s.pos.normalize_or_zero()).clamp(-1.0, 1.0).acos().to_degrees();
        let t0 = tilt(&s);
        run(&mut v, &mut s, &body, &Controls::default(), 3.0);
        assert!(tilt(&s) > t0 + 20.0, "Neigung {t0:.1}° → {:.1}°", tilt(&s));
    }

    #[test]
    fn tank_lying_on_its_side_stays_put() {
        let body = test_planet();
        let mut v = crate::vessel::VesselBuilder::new("Tank", PartLibrary::builtin().part("tank")).build();
        // Die Hülle ist ein Vieleck; der Tank liegt auf einer Fläche
        // (zwei Punkte je Stirnseite), nicht balancierend auf einer Kante.
        let eighth = std::f64::consts::PI / crate::geometry::HULL_SIDES as f64;
        let site = GeoPosition { latitude: 0.3, longitude: 1.0, altitude: 0.625 * eighth.cos() };
        let mut s = RigidState::at_rest(&body, site, 0.0);
        // um 90° kippen (Längsachse waagrecht), dann um die Längsachse auf eine Fläche drehen
        s.attitude = s.attitude
            * sim_core::Rotation::from_axis_angle(Vec3::X, std::f64::consts::FRAC_PI_2)
            * sim_core::Rotation::from_axis_angle(Vec3::Z, eighth);
        let start = s.geo_position(&body);
        run(&mut v, &mut s, &body, &Controls::default(), 10.0);
        let end = s.geo_position(&body);
        assert!(body.surface_distance(start, end) < 0.05, "rollt {} m", body.surface_distance(start, end));
        assert!((end.altitude - 0.625 * eighth.cos()).abs() < 0.02, "Höhe {}", end.altitude);
        let _ = mass_properties(&v);
    }
}
