//! Flug im Universum: Wechsel der Einflusssphäre (Patched Conics).
//!
//! Die Physik eines Schritts ist immer dieselbe ([`dynamics::step`]) und
//! rechnet relativ zum aktuellen Körper des Schiffs. Diese Schicht sorgt nur
//! dafür, dass das Schiff zum richtigen Körper gehört:
//!
//! - Verlässt es die Einflusssphäre seines Körpers, wird es dem
//!   Elternkörper zugeordnet.
//! - Betritt es die Einflusssphäre eines Kindkörpers (z. B. Muz), wird es
//!   diesem zugeordnet.
//!
//! Der Wechsel ist eine reine Umrechnung des Bezugspunkts: Position und
//! Geschwindigkeit relativ zum Stern bleiben exakt erhalten, Lage und
//! Drehung ohnehin (alle Körper teilen dieselben Achsen). Die Schrittweite
//! wird so begrenzt, dass kein Schritt über eine Grenze hinwegspringt.

use crate::dynamics::{self, Controls, RigidState, StepReport};
use crate::vessel::Vessel;
use sim_core::kepler::{self, StateVector};
use sim_core::{BodyId, CelestialBody, Universe};

/// Kürzester Schritt beim Annähern an eine Grenze (s). Der Wechsel geschieht
/// dann höchstens so weit hinter der Grenze; die Umrechnung selbst ist exakt.
pub const MIN_BOUNDARY_STEP: f64 = 1.0;

/// Ein Wechsel der Einflusssphäre: (vorher, nachher).
pub type SoiChange = (BodyId, BodyId);

/// Ein Schritt im Universum: Physik relativ zum aktuellen Körper, danach
/// gegebenenfalls Wechsel der Einflusssphäre.
pub fn step(v: &mut Vessel, s: &mut RigidState, c: &Controls, u: &Universe, dt: f64) -> (StepReport, Option<SoiChange>) {
    let report = dynamics::step(v, s, c, u.body(s.body), dt);
    let change = update_soi(s, u);
    (report, change)
}

/// Ordnet das Schiff dem richtigen Körper zu (höchstens ein Wechsel nach
/// außen und dann nach innen pro Aufruf).
pub fn update_soi(s: &mut RigidState, u: &Universe) -> Option<SoiChange> {
    let before = s.body;
    // Nach außen: Einflusssphäre verlassen
    loop {
        let body = u.body(s.body);
        let Some(parent) = body.parent else { break };
        if s.pos.length() <= body.soi_radius {
            break;
        }
        let rel = u.state_rel_parent(body.id, s.time);
        s.pos += rel.r;
        s.vel += rel.v;
        s.body = parent;
        s.contact_force = 0.0;
    }
    // Nach innen: in die Einflusssphäre eines Kindkörpers
    loop {
        let entered = u.children(s.body).find_map(|c| {
            let rel = u.state_rel_parent(c.id, s.time);
            ((s.pos - rel.r).length() < c.soi_radius).then_some((c.id, rel))
        });
        let Some((child, rel)) = entered else { break };
        s.pos -= rel.r;
        s.vel -= rel.v;
        s.body = child;
        s.contact_force = 0.0;
    }
    (s.body != before).then_some((before, s.body))
}

/// Höchstgeschwindigkeit auf einer Keplerbahn (am Pericenter).
fn max_speed(orbit: StateVector<sim_core::Inertial>, mu: f64) -> f64 {
    let el = kepler::elements(orbit, mu);
    let a = el.semi_major_axis;
    let rp = el.periapsis.max(1.0);
    let v2 = if a.is_finite() { mu * (2.0 / rp - 1.0 / a) } else { 2.0 * mu / rp };
    v2.max(0.0).sqrt().max(orbit.v.length())
}

/// Schrittweite wie [`dynamics::suggest_dt`], zusätzlich begrenzt durch die
/// Grenzen der Einflusssphären: Ein Schritt endet spätestens kurz nach dem
/// Erreichen einer Grenze.
pub fn suggest_dt(requested: f64, v: &Vessel, s: &RigidState, last: Option<&StepReport>, u: &Universe) -> f64 {
    let body = u.body(s.body);
    let mut dt = dynamics::suggest_dt(requested, v, s, last, body);
    let orbit = s.orbit();
    let el = kepler::elements(orbit, body.mu);
    let r = s.pos.length();
    let v_max = max_speed(orbit, body.mu);

    // Austritt: nur möglich, wenn die Bahn über die Grenze hinausreicht
    if body.soi_radius.is_finite() && el.apoapsis.is_none_or(|ap| ap >= body.soi_radius) {
        dt = dt.min(((body.soi_radius - r) / v_max).max(MIN_BOUNDARY_STEP));
    }
    // Eintritt in Kindkörper: nur, wenn sich die Radienbereiche überschneiden
    for child in u.children(s.body) {
        if !may_meet(&el, r, child, body.mu) {
            continue;
        }
        let rel = u.state_rel_parent(child.id, s.time);
        let gap = (s.pos - rel.r).length() - child.soi_radius;
        let closing = v_max + max_speed(rel, body.mu);
        dt = dt.min((gap / closing).max(MIN_BOUNDARY_STEP));
    }
    dt
}

/// Können sich Schiffsbahn und Einflusssphäre des Kindkörpers treffen?
fn may_meet(el: &kepler::Elements, r: f64, child: &CelestialBody, _mu: f64) -> bool {
    let Some(o) = child.orbit else { return false };
    let (child_lo, child_hi) = (o.semi_major_axis * (1.0 - o.eccentricity), o.semi_major_axis * (1.0 + o.eccentricity));
    let ship_lo = el.periapsis.min(r);
    let ship_hi = el.apoapsis.unwrap_or(f64::INFINITY).max(r);
    ship_hi >= child_lo - child.soi_radius && ship_lo <= child_hi + child.soi_radius
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::PartLibrary;
    use crate::scenario::{first_rocket, universe};
    use sim_core::{Inertial, Rotation, Vec3};
    use std::f64::consts::PI;

    fn in_orbit(u: &Universe, body: BodyId, r: f64, speed_factor: f64, angle: f64) -> RigidState {
        let mu = u.body(body).mu;
        let dir = Vec3::<Inertial>::new(angle.cos(), angle.sin(), 0.0);
        let tangent = Vec3::<Inertial>::new(-angle.sin(), angle.cos(), 0.0);
        RigidState {
            body,
            time: 0.0,
            pos: dir * r,
            vel: tangent * ((mu / r).sqrt() * speed_factor),
            attitude: Rotation::IDENTITY,
            omega: Vec3::ZERO,
            contact_force: 0.0,
        }
    }

    /// Antriebslos fliegen, bis `stop` erfüllt ist oder `t_max` erreicht.
    fn coast(u: &Universe, s: &mut RigidState, t_max: f64, mut stop: impl FnMut(&RigidState, Option<SoiChange>) -> bool) -> Vec<SoiChange> {
        let mut v = first_rocket(&PartLibrary::builtin());
        let mut last = None;
        let mut changes = Vec::new();
        while s.time < t_max {
            let dt = suggest_dt(t_max - s.time, &v, s, last.as_ref(), u);
            let (r, ch) = step(&mut v, s, &Controls::default(), u, dt);
            last = Some(r);
            if let Some(c) = ch {
                changes.push(c);
            }
            if stop(s, ch) {
                break;
            }
        }
        changes
    }

    #[test]
    fn escaping_zorb_switches_to_zun_without_jump() {
        let u = universe();
        let zorb = u.find("Zorb").unwrap();
        let zun = u.find("Zun").unwrap();
        // Fluchtbahn von Muz weg (sonst fängt Muz das Schiff unterwegs ein,
        // was ebenfalls richtig wäre, aber hier nicht geprüft wird)
        let muz = u.state_rel_parent(u.find("Muz").unwrap(), 0.0).r;
        let away = muz.y().atan2(muz.x()) + PI;
        let mut s = in_orbit(&u, zorb, 700_000.0, 2f64.sqrt() * 1.3, away);
        // Zustand relativ zum Stern vor und nach dem Wechsel vergleichen
        let mut before = None;
        let changes = coast(&u, &mut s, 1.0e7, |st, ch| {
            if ch.is_none() {
                let p = u.state_in_system(st.body, st.time);
                before = Some((st.time, st.pos + p.r, st.vel + p.v));
            }
            ch.is_some()
        });
        assert_eq!(changes, vec![(zorb, zun)]);
        // Unmittelbar nach dem Wechsel: dieselbe Lage relativ zu Zun wie vorher
        // plus ein Schritt Flug (keine Sprünge in Position oder Geschwindigkeit)
        let (t0, r0, v0) = before.unwrap();
        let dt = s.time - t0;
        assert!((s.pos - r0).length() < v0.length() * dt * 1.01 + 1.0);
        assert!((s.vel - v0).length() < 50.0, "Geschwindigkeitssprung {:?}", s.vel - v0);
    }

    #[test]
    fn hohmann_transfer_reaches_muz() {
        let u = universe();
        let (zorb, muz) = (u.find("Zorb").unwrap(), u.find("Muz").unwrap());
        let mu = u.body(zorb).mu;
        let r1 = 700_000.0;
        let r2 = u.body(muz).orbit.unwrap().semi_major_axis;
        let a_t = (r1 + r2) / 2.0;
        let t_transfer = PI * (a_t.powi(3) / mu).sqrt();
        let n_muz = (mu / r2.powi(3)).sqrt();
        // Muz muss beim Start um (π − n·T) vorauslaufen
        let muz_angle = {
            let p = u.state_rel_parent(muz, 0.0).r;
            p.y().atan2(p.x())
        };
        let ship_angle = muz_angle - (PI - n_muz * t_transfer);
        let v_transfer = (mu * (2.0 / r1 - 1.0 / a_t)).sqrt();
        let mut s = in_orbit(&u, zorb, r1, v_transfer / (mu / r1).sqrt(), ship_angle);

        let changes = coast(&u, &mut s, t_transfer * 1.5, |st, _| st.body == muz);
        assert_eq!(changes.first(), Some(&(zorb, muz)), "Muz nicht erreicht");
        // Eintritt nahe am Rand der Einflusssphäre, nicht weit dahinter
        let r = s.pos.length();
        let soi = u.body(muz).soi_radius;
        assert!(r < soi && r > soi * 0.97, "Eintritt bei {r} m (SOI {soi} m)");
        // Vor der Ankunft am Zielpunkt: die Grenze von Muz liegt 2 430 km davor,
        // und nahe dem Apozentrum ist das Schiff langsam.
        assert!(s.time < t_transfer && s.time > t_transfer * 0.6, "{} s statt < {t_transfer} s", s.time);
    }

    #[test]
    fn low_orbit_is_not_slowed_by_soi_checks() {
        // Stabile Bahn weit innerhalb aller Grenzen: Schritte bleiben lang
        let u = universe();
        let zorb = u.find("Zorb").unwrap();
        let s = in_orbit(&u, zorb, 700_000.0, 1.0, 0.0);
        let v = first_rocket(&PartLibrary::builtin());
        let dt = suggest_dt(1.0e6, &v, &s, None, &u);
        assert!(dt >= 1.0e6 * 0.999, "Schritt {dt} s");
    }
}
