//! Luftkräfte, für alle Teile mit derselben Rechnung.
//!
//! Grundlage ist das Polygonmodell des ganzen Schiffs
//! ([`crate::geometry::HullGeometry`]): welche Flächen frei angeströmt
//! werden und wo das Schiff breiter wird. Teile unterscheiden sich nur durch
//! ihre Hülle (und optionale Abstimmfaktoren), nicht durch ihr Luftmodell.
//!
//! Die Anströmung wird in einen Anteil längs der Schiffsachse und einen
//! Queranteil zerlegt (Modell nach Allen und Perkins, Standard in der
//! Raketen-Aerodynamik). So wird nichts doppelt gezählt:
//!
//! 1. **Längsströmung** (Staudruck ∝ cos²α): Druck auf frei angeströmte
//!    Flächen (modifizierte Newton-Theorie), Sog auf freien Rückflächen
//!    (Heckwiderstand). Daraus folgt der Widerstand im Geradeausflug.
//! 2. **Querströmung** (∝ sin²α): dieselbe Druckrechnung mit dem Queranteil.
//!    Liefert den Querwiderstand bei großem Anstellwinkel.
//! 3. **Auftrieb** (∝ sin 2α) nach der Schlankkörpertheorie: Querkraft
//!    entsteht dort, wo die Spannweite des Schiffs in Strömungsrichtung
//!    zunimmt (Nasenspitze, Vorderkanten der Finnen). Korrigiert nach
//!    Seitenverhältnis (Helmbold). Gilt für Rumpf und Finnen gleichermaßen.
//! 4. **Reibung** auf allen umströmten Flächen.
//!
//! Alle Kräfte sind proportional zur Luftdichte: Vakuum ist Dichte 0.
//! Ein Ruderausschlag dreht die Anströmung, die das Teil sieht; sonst
//! gibt es für Steuerflächen keine Sonderregel.

use crate::geometry::HullGeometry;
use crate::mass::mass_properties;
use crate::vessel::Vessel;
use glam::{DQuat, DVec3};
use sim_core::{AtmosphereSample, Vec3, VesselFrame};

/// Höchster Druckbeiwert unterhalb der Schallgeschwindigkeit (stumpfe Stirnfläche).
pub const CP_MAX_SUBSONIC: f64 = 0.8;
/// Höchster Druckbeiwert bei Hyperschall (modifizierte Newton-Theorie, γ = 1,4).
pub const CP_MAX_HYPERSONIC: f64 = 1.84;
/// Sog auf freien Rückflächen bei Längsströmung (achsensymmetrische
/// Ablösung am Heck, Unterschall).
pub const CP_BASE: f64 = 0.15;
/// Sog auf der Leeseite bei Querströmung (ebene Ablösung wie am quer
/// angeströmten Zylinder, Unterschall). Ergibt Cd ≈ 1,2 für Zylinder.
pub const CP_BASE_CROSS: f64 = 0.85;
/// Reibungsbeiwert (turbulente Grenzschicht).
pub const CF: f64 = 0.004;

fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Staupunktdruck-Beiwert über der Machzahl, mit Anstieg im Transschall.
pub fn cp_max(mach: f64) -> f64 {
    let bump = (-((mach - 1.1) / 0.3).powi(2)).exp();
    CP_MAX_SUBSONIC * (1.0 + 0.6 * bump) + (CP_MAX_HYPERSONIC - CP_MAX_SUBSONIC) * smoothstep((mach - 0.8) / 2.2)
}

/// Sog-Beiwert auf Rückflächen; fällt im Überschall wie 1/M².
pub fn cp_base(base: f64, mach: f64) -> f64 {
    base * (1.0 / (mach * mach)).min(1.0)
}

/// Eine Luftkraft mit Angriffspunkt (Schiffssystem).
#[derive(Debug, Clone, Copy)]
pub struct AeroForce {
    pub force: Vec3<VesselFrame>,
    pub point: Vec3<VesselFrame>,
}

/// Luftkräfte auf Teil `i`.
///
/// `air`: Geschwindigkeit der Luft relativ zum Teil (Schiffssystem), also
/// −(Eigengeschwindigkeit des Teils gegenüber der Luft).
/// `deflection`: Ruderausschlag in rad (Drehung um die Spannweitenachse).
pub fn part_forces(v: &Vessel, geo: &HullGeometry, i: usize, air: Vec3<VesselFrame>, atm: &AtmosphereSample, deflection: f64) -> [AeroForce; 3] {
    let p = &v.parts[i];
    let com = p.com();
    let zero = AeroForce { force: Vec3::ZERO, point: com };
    let speed = air.length();
    if speed == 0.0 || atm.density == 0.0 {
        return [zero; 3];
    }

    // Ausgeschlagenes Ruder: Strömung im gedrehten Teil betrachten, Kraft zurückdrehen
    let hinge = p.orientation.apply(Vec3::X).raw();
    let turn = DQuat::from_axis_angle(hinge, deflection);
    let u = turn.inverse() * air.raw() / speed;

    let q = 0.5 * atm.density * speed * speed;
    let mach = speed / atm.speed_of_sound;
    let ua = u.z;
    let w = DVec3::new(u.x, u.y, 0.0);
    let ws = w.length();

    // Druck: Längsanteil und Queranteil mit derselben Rechnung
    let mut f = DVec3::ZERO;
    let mut pressure = |flow: DVec3, q_part: f64, m_part: f64, base: f64, up: usize, down: usize| {
        for (j, panel) in geo.panels.iter().enumerate().filter(|(_, pn)| pn.part == i) {
            let nf = panel.normal.dot(flow);
            if nf < 0.0 {
                f -= panel.normal * (q_part * cp_max(m_part) * nf * nf * panel.area * geo.visible[up][j]);
            } else if nf > 0.0 {
                f += panel.normal * (q_part * cp_base(base, m_part) * nf * panel.area * geo.visible[down][j]);
            }
        }
    };
    if ua.abs() > 1e-12 {
        let flow = DVec3::Z * ua.signum();
        // Kamera stromaufwärts sieht die angeströmten Flächen
        let (up, down) = if ua < 0.0 { (0, 1) } else { (1, 0) };
        pressure(flow, q * ua * ua, mach * ua.abs(), CP_BASE, up, down);
    }
    if ws > 1e-12 {
        let g = w / ws;
        pressure(g, q * ws * ws, mach * ws, CP_BASE_CROSS, geo.nearest_cross_direction(-g), geo.nearest_cross_direction(g));
    }
    // Reibung entlang der Oberfläche
    for (j, panel) in geo.panels.iter().enumerate().filter(|(_, pn)| pn.part == i) {
        let t = u - panel.normal * panel.normal.dot(u);
        f += t * (q * CF * panel.area * geo.wetted[j]);
    }
    f *= p.def.drag_factor;

    // Auftrieb nach Schlankkörpertheorie: N = q · 2 cosα · sinα · ΔS
    let end = if ua <= 0.0 { 0 } else { 1 };
    let (sx, sy) = (geo.lift[end][0][i], geo.lift[end][1][i]);
    let k = 2.0 * q * ua.abs() * geo.aspect_correction[i] * p.def.lift_factor;
    let lift_x = DVec3::X * (k * u.x * sx.delta_area);
    let lift_y = DVec3::Y * (k * u.y * sy.delta_area);

    let back = |d: DVec3| Vec3::<VesselFrame>::from_raw(turn * d);
    let at = |z: f64| Vec3::new(com.x(), com.y(), z);
    [
        AeroForce { force: back(f), point: com },
        AeroForce { force: back(lift_x), point: at(sx.center_z) },
        AeroForce { force: back(lift_y), point: at(sy.center_z) },
    ]
}

/// Summe der Luftkräfte und Moment um den Schwerpunkt für eine gleichförmige
/// Anströmung (Schiff dreht sich nicht). Für Auswertung und Tests.
pub fn vessel_forces(v: &Vessel, air: Vec3<VesselFrame>, atm: &AtmosphereSample) -> (Vec3<VesselFrame>, Vec3<VesselFrame>) {
    let geo = v.geometry();
    let com = mass_properties(v).com;
    let mut force = Vec3::ZERO;
    let mut moment = Vec3::ZERO;
    for i in 0..v.parts.len() {
        for a in part_forces(v, &geo, i, air, atm, 0.0) {
            force += a.force;
            moment += (a.point - com).cross(a.force);
        }
    }
    (force, moment)
}

/// Druckpunkt auf der Schiffsachse (z im Schiffssystem) bei Anstellwinkel
/// `alpha` (Spitze voran, Anströmung in der Y-Z-Ebene). `None`, wenn keine
/// Querkraft entsteht.
pub fn center_of_pressure(v: &Vessel, atm: &AtmosphereSample, speed: f64, alpha: f64) -> Option<f64> {
    let air = Vec3::new(0.0, -speed * alpha.sin(), -speed * alpha.cos());
    let (f, m) = vessel_forces(v, air, atm);
    if f.y().abs() < 1e-9 {
        return None;
    }
    Some(mass_properties(v).com.z() - m.x() / f.y())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::PartLibrary;
    use crate::scenario::first_rocket;
    use crate::vessel::VesselBuilder;
    use sim_core::Atmosphere;
    use std::f64::consts::PI;

    fn sea_level() -> AtmosphereSample {
        Atmosphere::earthlike().sample(0.0)
    }

    fn rocket(fins: bool) -> Vessel {
        let mut v = first_rocket(&PartLibrary::builtin());
        if !fins {
            v.parts.retain(|p| p.def.control_surface().is_none());
            v.invalidate_geometry();
        }
        v
    }

    /// Normalkraftgradient dCN/dα, bezogen auf q·A (A = π r², r = 0,625 m).
    fn cn_alpha(v: &Vessel, atm: &AtmosphereSample, speed: f64) -> f64 {
        let a = 0.5f64.to_radians();
        let air = Vec3::new(0.0, -speed * a.sin(), -speed * a.cos());
        let (f, _) = vessel_forces(v, air, atm);
        let q = 0.5 * atm.density * speed * speed;
        -f.y() / (q * PI * 0.625 * 0.625 * a)
    }

    #[test]
    fn body_normal_force_matches_munk_and_barrowman() {
        // Kapsel + Tank ohne Finnen: CNα = 2 (Munk/Barrowman, unabhängig von
        // der Nasenform). Kleine Abweichung durch das 12-Eck und die Korrektur.
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Rumpf", lib.part("kapsel"));
        b.stack_below(0, lib.part("tank"));
        let v = b.build();
        let cn = cn_alpha(&v, &sea_level(), 100.0);
        assert!((cn - 2.0).abs() < 0.2, "CNα = {cn}");
    }

    #[test]
    fn nose_center_of_pressure_matches_barrowman() {
        // Rakete ohne Finnen bei sehr kleinem Anstellwinkel: nur die Kapsel
        // trägt. Barrowman: Stirnfläche (r = 0,25 m) bei z = +0,6 mit
        // CNα = 2·(0,25/0,625)², Kegelübergang mit CNα = 2·(1 − 0,16) und
        // Druckpunkt L/3·(1 + (1 − d₁/d₂)/(1 − (d₁/d₂)²)) hinter der Spitze.
        let v = rocket(false);
        let ratio: f64 = 0.25 / 0.625;
        let (cn_face, cn_cone) = (2.0 * ratio * ratio, 2.0 * (1.0 - ratio * ratio));
        let cone_cp = 0.6 - 1.2 / 3.0 * (1.0 + (1.0 - ratio) / (1.0 - ratio * ratio));
        let barrowman = (cn_face * 0.6 + cn_cone * cone_cp) / (cn_face + cn_cone);
        let cp = center_of_pressure(&v, &sea_level(), 100.0, 0.1f64.to_radians()).unwrap();
        assert!((cp - barrowman).abs() < 0.3, "DP {cp:.2} m, Barrowman {barrowman:.2} m");
    }

    #[test]
    fn cylinder_crossflow_drag() {
        // Tank quer angeströmt: Cd ≈ 1,2 (Messwert Zylinder im Unterschall)
        let v = VesselBuilder::new("Tank", PartLibrary::builtin().part("tank")).build();
        let atm = sea_level();
        let speed = 50.0;
        let (f, _) = vessel_forces(&v, Vec3::new(-speed, 0.0, 0.0), &atm);
        let q = 0.5 * atm.density * speed * speed;
        let cd = -f.x() / (q * 1.25 * 1.9);
        assert!(cd > 1.0 && cd < 1.4, "Cd quer = {cd}");
    }

    #[test]
    fn blunt_cylinder_axial_drag() {
        // Stumpfer Zylinder längs angeströmt: Cd ≈ 0,8–1,0 (Messwerte)
        let v = VesselBuilder::new("Tank", PartLibrary::builtin().part("tank")).build();
        let atm = sea_level();
        let speed = 50.0;
        let (f, _) = vessel_forces(&v, Vec3::new(0.0, 0.0, -speed), &atm);
        let q = 0.5 * atm.density * speed * speed;
        let cd = -f.z() / (q * PI * 0.625 * 0.625);
        assert!(cd > 0.8 && cd < 1.1, "Cd längs = {cd}");
    }

    #[test]
    fn stacking_hides_faces_drag_barely_grows() {
        let lib = PartLibrary::builtin();
        let one = VesselBuilder::new("1", lib.part("tank")).build();
        let mut b = VesselBuilder::new("2", lib.part("tank"));
        b.stack_below(0, lib.part("tank"));
        let two = b.build();
        let atm = sea_level();
        let air = Vec3::new(0.0, 0.0, -100.0);
        let d1 = -vessel_forces(&one, air, &atm).0.z();
        let d2 = -vessel_forces(&two, air, &atm).0.z();
        // Nur die Reibung des zweiten Tanks kommt dazu (wenige Prozent)
        assert!(d2 > d1 && d2 < d1 * 1.05, "ein Tank {d1} N, zwei {d2} N");
    }

    #[test]
    fn capsule_nose_first_has_less_drag_than_tail_first() {
        let v = VesselBuilder::new("Kapsel", PartLibrary::builtin().part("kapsel")).build();
        let atm = sea_level();
        let nose = vessel_forces(&v, Vec3::new(0.0, 0.0, -100.0), &atm).0.length();
        let tail = vessel_forces(&v, Vec3::new(0.0, 0.0, 100.0), &atm).0.length();
        assert!(nose < tail, "Spitze voran {nose} N, Boden voran {tail} N");
    }

    #[test]
    fn fins_move_center_of_pressure_aft() {
        let atm = Atmosphere::earthlike().sample(8_000.0);
        let with = rocket(true);
        let without = rocket(false);
        let com_with = mass_properties(&with).com.z();
        let com_without = mass_properties(&without).com.z();
        for deg in [2.0f64, 10.0, 30.0] {
            let cp_with = center_of_pressure(&with, &atm, 300.0, deg.to_radians()).unwrap();
            let cp_without = center_of_pressure(&without, &atm, 300.0, deg.to_radians()).unwrap();
            // Mit Finnen: Druckpunkt deutlich hinter dem Schwerpunkt (stabil)
            assert!(cp_with < com_with - 0.5, "α={deg}°: DP {cp_with:.2} SP {com_with:.2}");
            // Ohne Finnen: Druckpunkt vor dem Schwerpunkt (instabil)
            assert!(cp_without > com_without + 0.5, "α={deg}°: DP {cp_without:.2} SP {com_without:.2}");
        }
    }

    #[test]
    fn forces_vanish_in_vacuum_and_scale_with_density() {
        let v = rocket(true);
        let air = Vec3::new(10.0, 20.0, -300.0);
        let mut atm = Atmosphere::earthlike().sample(10_000.0);
        let f1 = vessel_forces(&v, air, &atm).0;
        atm.density *= 1e-9;
        let f2 = vessel_forces(&v, air, &atm).0;
        assert!((f2 * 1e9 - f1).length() < 1e-9 * f1.length());
        atm.density = 0.0;
        assert_eq!(vessel_forces(&v, air, &atm).0, Vec3::ZERO);
    }

    #[test]
    fn deflected_fin_creates_side_force() {
        let v = rocket(true);
        let geo = v.geometry();
        let fin = v.parts.iter().position(|p| p.def.control_surface().is_some()).unwrap();
        let atm = sea_level();
        let air = Vec3::new(0.0, 0.0, -100.0);
        let sum = |d: f64| part_forces(&v, &geo, fin, air, &atm, d).iter().map(|a| a.force).sum::<Vec3<VesselFrame>>();
        let lateral = |f: Vec3<VesselFrame>| f.x().hypot(f.y());
        assert!(lateral(sum(0.0)) < 1e-6);
        let small = lateral(sum(0.05));
        let large = lateral(sum(0.1));
        assert!(small > 0.0 && (large / small - 2.0).abs() < 0.1, "nicht linear: {small} {large}");
    }

    #[test]
    fn cp_max_is_smooth() {
        let mut last = cp_max(0.0);
        for k in 1..=5_000 {
            let m = k as f64 * 0.001;
            let c = cp_max(m);
            assert!((c - last).abs() < 0.01, "Sprung bei M {m}");
            last = c;
        }
        assert!((cp_max(20.0) - CP_MAX_HYPERSONIC).abs() < 1e-6);
    }
}
