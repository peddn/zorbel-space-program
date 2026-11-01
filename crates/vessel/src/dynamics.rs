//! Bewegung eines Schiffs als starrer Körper.
//!
//! Es gibt genau **eine** Bewegungsgleichung und **ein** Schrittverfahren
//! ([`step`]). In jedem Schritt wirken immer alle Kräfte: Gravitation,
//! Schub, Luft, Reaktionsräder. Ob das Schiff am Boden, in der Atmosphäre
//! oder im Vakuum ist, ändert nur die Größe der Kräfte, nie den Codepfad.
//!
//! Verfahren: Kick-Drift-Kick. Die Gravitation des Zentralkörpers wird
//! exakt über die Keplerbahn gelöst ("Drift"), alle anderen Kräfte wirken
//! als halbe Stöße davor und danach. Sind die anderen Kräfte null, ist der
//! Schritt exakt die Keplerbahn; deshalb darf er dann beliebig lang sein.
//! Die Schrittweite ([`suggest_dt`]) ist eine rein numerische Wahl und
//! ändert die Physik nicht.
//!
//! Der Boden ist eine Zwangsbedingung am Ende jedes Schritts
//! ([`crate::ground`]); seine Stützkraft geht wie jede andere äußere Kraft
//! in die Lastberechnung ein.

use crate::aero::{AeroForce, part_forces};
use crate::geometry::HullGeometry;
use crate::ground;
use crate::loads::{ExternalLoad, LoadAnalysis, Motion, joint_loads, motion_from_loads};
use crate::mass::{MassProperties, mass_properties};
use crate::vessel::Vessel;
use sim_core::kepler::{self, StateVector};
use sim_core::{AtmosphereSample, BodyId, CelestialBody, GeoPosition, Inertial, Rotation, Vec3, VesselFrame};

/// Kleinster Schritt bei voller Dynamik (Schub, dichte Luft).
pub const DT_PHYSICS: f64 = 0.02;

#[derive(Debug, Clone, Copy)]
pub struct RigidState {
    /// Körper, auf den sich Position und Geschwindigkeit beziehen (aktuelle
    /// Einflusssphäre). Wechselt in [`crate::flight`].
    pub body: BodyId,
    /// Simulationszeit in s. Bestimmt u. a. die Drehlage des Planeten.
    pub time: f64,
    /// Position des Schwerpunkts relativ zum Zentralkörper.
    pub pos: Vec3<Inertial>,
    pub vel: Vec3<Inertial>,
    pub attitude: Rotation<VesselFrame, Inertial>,
    /// Winkelgeschwindigkeit im Schiffssystem.
    pub omega: Vec3<VesselFrame>,
    /// Geglättete Bodenkontaktkraft in N (siehe [`crate::ground`]).
    pub contact_force: f64,
}

impl RigidState {
    pub fn orbit(&self) -> StateVector<Inertial> {
        StateVector::new(self.pos, self.vel)
    }
    /// Richtung der Spitze im Inertialsystem.
    pub fn nose(&self) -> Vec3<Inertial> {
        self.attitude.apply(Vec3::Z)
    }
    pub fn altitude(&self, body: &CelestialBody) -> f64 {
        self.pos.length() - body.radius
    }
    /// Geschwindigkeit relativ zur mitdrehenden Oberfläche (und zur Luft).
    pub fn surface_velocity(&self, body: &CelestialBody) -> Vec3<Inertial> {
        self.vel - body.corotation_velocity(self.pos)
    }
    pub fn geo_position(&self, body: &CelestialBody) -> GeoPosition {
        body.geo_position(self.pos, self.time)
    }

    /// Ruhend über dem Ort `site` (Höhe = Höhe des Schwerpunkts), Spitze
    /// senkrecht nach oben. "Ruhend" heißt: dreht mit dem Planeten mit,
    /// in Position, Geschwindigkeit und Winkelgeschwindigkeit.
    pub fn at_rest(body: &CelestialBody, site: GeoPosition, time: f64) -> Self {
        let pos = body.orientation(time).apply(body.geo_to_body_fixed(site));
        let up = pos.normalize_or_zero();
        let attitude: Rotation<VesselFrame, Inertial> =
            Rotation::from_quat(glam::DQuat::from_rotation_arc(glam::DVec3::Z, up.raw()));
        RigidState {
            body: body.id,
            time,
            pos,
            vel: body.corotation_velocity(pos),
            attitude,
            omega: attitude.inverse().apply(body.angular_velocity()),
            contact_force: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Controls {
    /// 0..1
    pub throttle: f64,
    /// Drehkommando je Achse, -1..1 (Schiffssystem). Wirkt auf
    /// Reaktionsräder, Schubvektor und Steuerflächen gleichzeitig.
    pub rotation: Vec3<VesselFrame>,
}

#[derive(Debug, Clone)]
pub struct StepReport {
    pub loads: LoadAnalysis,
    /// Gesamtschub in N.
    pub thrust: f64,
    /// Gesamte Luftkraft in N.
    pub aero_force: f64,
    /// Betrag der Eigenbeschleunigung in m/s².
    pub proper_accel: f64,
    pub air: AtmosphereSample,
    pub mach: f64,
    /// Staudruck ½ρv² in Pa.
    pub dynamic_pressure: f64,
    /// Aktive Triebwerke ohne Treibstoff.
    pub flameout: bool,
    /// Bodenkontakt in diesem Schritt.
    pub touching_ground: bool,
    /// Annäherungsgeschwindigkeit an den Boden beim Aufsetzen (m/s).
    pub impact_speed: f64,
    /// Teil, das beim Aufprall zerstört wurde, mit Aufprallgeschwindigkeit.
    pub crash: Option<(usize, f64)>,
    /// Geglättete Bodenkontaktkraft in N.
    pub ground_force: f64,
}

/// Alle Kräfte für einen Zustand.
struct Evaluation {
    ext: Vec<ExternalLoad>,
    burn: Vec<(usize, f64)>,
    flameout: bool,
    motion: Motion,
    thrust: f64,
    aero_force: f64,
    air: AtmosphereSample,
    airspeed: f64,
}

/// Probeausschlag, um die Wirkrichtung einer Steuerfläche zu bestimmen.
const DEFLECTION_PROBE: f64 = 0.05;

fn evaluate(v: &Vessel, geo: &HullGeometry, mp: &MassProperties, s: &RigidState, c: &Controls, body: &CelestialBody) -> Evaluation {
    let air = body.atmosphere.sample(s.altitude(body));
    let cmd = Vec3::<VesselFrame>::new(
        c.rotation.x().clamp(-1.0, 1.0),
        c.rotation.y().clamp(-1.0, 1.0),
        c.rotation.z().clamp(-1.0, 1.0),
    );
    // Die Luft dreht mit dem Planeten; maßgeblich ist die Anströmung.
    let v_body = s.attitude.inverse().apply(s.surface_velocity(body));
    // Drehung relativ zur Luft: die Luft am Ort jedes Teils dreht ebenfalls mit.
    let omega_rel = s.omega - s.attitude.inverse().apply(body.angular_velocity());

    let mut ext = Vec::new();
    let mut burn = Vec::new();
    let mut flameout = false;
    let mut thrust_total = 0.0;
    let mut aero_total = Vec3::<VesselFrame>::ZERO;

    for (i, p) in v.parts.iter().enumerate() {
        // Reaktionsrad
        let wheel = p.def.command_torque();
        if wheel > 0.0 && cmd.length_squared() > 0.0 {
            ext.push(ExternalLoad { part: i, force: Vec3::ZERO, point: p.com(), torque: cmd * wheel });
        }

        // Triebwerk
        if let Some(eng) = p.def.engine()
            && p.engine_active
            && c.throttle > 0.0
        {
            let fuel: f64 = v.fuel_sources(i).iter().map(|&t| v.parts[t].fuel).sum();
            if fuel <= 0.0 {
                flameout = true;
            } else {
                let throttle = c.throttle.clamp(0.0, 1.0);
                let thrust = eng.thrust(throttle, air.pressure);
                // Schubvektor so auslenken, dass das Moment dem Kommando folgt.
                let lever = p.origin - mp.com;
                let side = if lever.z() <= 0.0 { 1.0 } else { -1.0 };
                let g = eng.gimbal_rad.tan() * side;
                let axis = p.orientation.apply(Vec3::Z);
                let dir = (axis + Vec3::new(-cmd.y() * g, cmd.x() * g, 0.0)).normalize_or_zero();
                ext.push(ExternalLoad { part: i, force: dir * thrust, point: p.origin, torque: Vec3::ZERO });
                burn.push((i, eng.mass_flow(throttle)));
                thrust_total += thrust;
            }
        }

        // Luftkräfte des Teils, mit der Anströmung an seinem eigenen Ort
        let lever = p.com() - mp.com;
        let air_rel = -(v_body + omega_rel.cross(lever));
        let forces_at = |defl: f64| part_forces(v, geo, i, air_rel, &air, defl);
        let torque_of = |fs: &[AeroForce; 3]| fs.iter().map(|a| (a.point - mp.com).cross(a.force)).sum::<Vec3<VesselFrame>>();
        let deflection = match p.def.control_surface() {
            Some(cs) if cmd.length_squared() > 0.0 => {
                // Welche Drehung erzeugt ein positiver Ausschlag gerade?
                let d_torque = torque_of(&forces_at(DEFLECTION_PROBE)) - torque_of(&forces_at(0.0));
                let len = d_torque.length();
                if len > 0.0 { cs.max_deflection_rad * cmd.dot(d_torque / len).clamp(-1.0, 1.0) } else { 0.0 }
            }
            _ => 0.0,
        };
        for a in forces_at(deflection) {
            if a.force != Vec3::ZERO {
                aero_total += a.force;
                ext.push(ExternalLoad { part: i, force: a.force, point: a.point, torque: Vec3::ZERO });
            }
        }
    }

    let motion = motion_from_loads(mp, &ext, s.omega);
    Evaluation {
        ext,
        burn,
        flameout,
        motion,
        thrust: thrust_total,
        aero_force: aero_total.length(),
        air,
        airspeed: v_body.length(),
    }
}

/// Ein Zeitschritt. Gilt überall: am Boden, in der Luft, im Vakuum.
pub fn step(v: &mut Vessel, s: &mut RigidState, c: &Controls, body: &CelestialBody, dt: f64) -> StepReport {
    debug_assert_eq!(s.body, body.id, "Zustand bezieht sich auf einen anderen Körper");
    let geo = v.geometry();
    let mp0 = mass_properties(v);
    let ev = evaluate(v, &geo, &mp0, s, c, body);

    // Kick
    s.vel += s.attitude.apply(ev.motion.accel) * (dt / 2.0);
    s.omega += ev.motion.alpha * (dt / 2.0);

    // Drift: exakte Keplerbahn, freie Drehung
    s.time += dt;
    let o = kepler::propagate(s.orbit(), body.mu, dt);
    s.pos = o.r;
    s.vel = o.v;
    if s.omega != Vec3::ZERO {
        s.attitude = (s.attitude * Rotation::from_scaled_axis(s.omega * dt)).renormalize();
    }

    // Treibstoff, anteilig aus allen erreichbaren Tanks
    for &(engine, mdot) in &ev.burn {
        let tanks = v.fuel_sources(engine);
        let total: f64 = tanks.iter().map(|&t| v.parts[t].fuel).sum();
        let draw = (mdot * dt).min(total);
        if total > 0.0 {
            for t in tanks {
                let share = v.parts[t].fuel / total;
                v.parts[t].fuel = (v.parts[t].fuel - draw * share).max(0.0);
            }
        }
    }
    let mp1 = if ev.burn.is_empty() { mp0 } else { mass_properties(v) };
    if !ev.burn.is_empty() {
        // Schwerpunkt wandert innerhalb des Schiffs
        s.pos += s.attitude.apply(mp1.com - mp0.com);
    }

    // Kick
    let ev2 = evaluate(v, &geo, &mp1, s, c, body);
    s.vel += s.attitude.apply(ev2.motion.accel) * (dt / 2.0);
    s.omega += ev2.motion.alpha * (dt / 2.0);

    // Boden: Zwangsbedingung am Schrittende
    let contact = ground::resolve(v, &mp1, s, body, dt);

    // Lasten aus allen äußeren Kräften am Schrittende, inklusive Boden
    let mut ext = ev2.ext;
    ext.extend(contact.loads);
    let motion = motion_from_loads(&mp1, &ext, s.omega);
    let loads = joint_loads(v, &mp1, &motion, &ext);

    let q = 0.5 * ev.air.density * ev.airspeed * ev.airspeed;
    StepReport {
        loads,
        thrust: ev.thrust,
        aero_force: ev.aero_force,
        proper_accel: motion.accel.length().max(ev.motion.accel.length()),
        mach: ev.airspeed / ev.air.speed_of_sound,
        dynamic_pressure: q,
        air: ev.air,
        flameout: ev.flameout,
        touching_ground: contact.touching,
        impact_speed: contact.impact_speed,
        crash: contact.crash,
        ground_force: s.contact_force,
    }
}

/// Schrittweite für den nächsten Schritt, höchstens `requested`.
///
/// - Starke Störkräfte (Schub, dichte Luft) → kurze Schritte, die Länge
///   wächst stetig, je schwächer sie werden.
/// - Schnelle Drehung → kurze Schritte (höchstens ~1° pro Schritt).
/// - Oberhalb der Atmosphäre: Schritt endet spätestens beim Eintritt, damit
///   kein langer Kepler-Schritt durch die Lufthülle "hindurchspringt".
/// - Ebenso endet er, bevor irgendein Teil den Boden erreichen kann; in
///   Bodennähe nur kurze Schritte.
/// - Ohne Vorwissen (erster Schritt) in Luft oder Bodennähe: kurz.
pub fn suggest_dt(requested: f64, v: &Vessel, s: &RigidState, last: Option<&StepReport>, body: &CelestialBody) -> f64 {
    let mut dt = requested;
    let r = s.pos.length();
    let ground = body.radius + ground::extent(v, &mass_properties(v)) + 1.0;
    let quiet = body.quiet_radius().max(ground);
    if r <= quiet && last.is_none() {
        dt = dt.min(DT_PHYSICS);
    }
    if r <= ground {
        dt = dt.min(DT_PHYSICS);
    } else if let Some(t) = kepler::time_until_below(s.orbit(), body.mu, ground, dt) {
        dt = dt.min(t.max(DT_PHYSICS));
    }
    if let Some(r) = last
        && r.proper_accel > 0.0
    {
        dt = dt.min(DT_PHYSICS * (1.0 / r.proper_accel).clamp(1.0, 500.0));
    }
    let w = s.omega.length();
    if w > 0.0 {
        dt = dt.min((0.02 / w).max(DT_PHYSICS));
    }
    if r > quiet
        && let Some(t) = kepler::time_until_below(s.orbit(), body.mu, quiet, dt)
    {
        dt = dt.min(t.max(DT_PHYSICS));
    }
    dt
}

#[derive(Debug, Clone)]
pub struct Separated {
    pub vessel: Vessel,
    pub state: RigidState,
}

/// Trennt den Teilbaum unter `child` ab (Entkoppler oder Bruch).
/// `impulse` (N·s) drückt beide Teile entlang der Verbindungsachse auseinander.
/// Rückgabe: (Rest, abgetrenntes Teil).
pub fn separate(v: &Vessel, s: &RigidState, child: usize, impulse: f64) -> (Separated, Separated) {
    let axis = v.parts[child].joint.expect("Wurzel kann nicht getrennt werden").axis;
    let mp = mass_properties(v);
    let (rest, detached) = v.split_off(child);

    let piece = |vessel: Vessel, push: f64| {
        let mpp = mass_properties(&vessel);
        let d = mpp.com - mp.com;
        let state = RigidState {
            body: s.body,
            time: s.time,
            pos: s.pos + s.attitude.apply(d),
            vel: s.vel + s.attitude.apply(s.omega.cross(d) + axis * (push / mpp.mass)),
            attitude: s.attitude,
            omega: s.omega,
            contact_force: 0.0,
        };
        Separated { vessel, state }
    };
    (piece(rest, -impulse), piece(detached, impulse))
}

/// Löst die nächste Stufe aus: Triebwerke zünden, Entkoppler trennen.
/// `ship`/`s` behalten den oberen Teil; zurück kommen die abgeworfenen Teile.
pub fn activate_stage(ship: &mut Vessel, s: &mut RigidState) -> Vec<Separated> {
    let mut dropped = Vec::new();
    // Von hinten, damit die Indizes der übrigen Entkoppler gültig bleiben.
    for d in ship.activate_next_stage().into_iter().rev() {
        let impulse = ship.parts[d].def.decoupler_impulse().unwrap_or(0.0);
        let (rest, lower) = separate(ship, s, d, impulse);
        *ship = rest.vessel;
        *s = rest.state;
        dropped.push(lower);
    }
    dropped
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::PartLibrary;
    use crate::vessel::VesselBuilder;
    use sim_core::Atmosphere;

    fn two_stage() -> Vessel {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Zweistufig", lib.part("kapsel"));
        let t2 = b.stack_below(0, lib.part("tank"));
        let e2 = b.stack_below(t2, lib.part("triebwerk"));
        let dec = b.stack_below(e2, lib.part("entkoppler"));
        let t1 = b.stack_below(dec, lib.part("tank"));
        let e1 = b.stack_below(t1, lib.part("triebwerk"));
        b.radial_symmetric(t1, lib.part("steuerflaeche"), 4, -0.5);
        b.set_stage(e1, 1).set_stage(dec, 0).set_stage(e2, 0);
        b.build()
    }

    fn some_state() -> RigidState {
        RigidState {
            body: BodyId(0),
            time: 0.0,
            pos: Vec3::new(700_000.0, 0.0, 0.0),
            vel: Vec3::new(0.0, 2_200.0, 0.0),
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(0.3)),
            omega: Vec3::new(0.05, 0.0, 0.2),
            contact_force: 0.0,
        }
    }

    fn planet(atmosphere: Atmosphere) -> CelestialBody {
        CelestialBody {
            id: BodyId(0),
            parent: None,
            orbit: None,
            appearance: sim_core::Appearance::Terran,
            name: "Test".into(),
            mu: 3.5316e12,
            radius: 600_000.0,
            soi_radius: 8.4e7,
            atmosphere,
            rotation_rate: 0.0,
            rotation_offset: 0.0,
        }
    }

    #[test]
    fn staging_order_and_fuel_isolation() {
        let mut v = two_stage();
        assert!(v.activate_next_stage().is_empty()); // Stufe 1: nur Zündung unten
        assert!(v.parts[5].engine_active && !v.parts[2].engine_active);
        assert_eq!(v.fuel_sources(5), vec![4]);
        assert_eq!(v.fuel_sources(2), vec![1]);
        assert_eq!(v.activate_next_stage(), vec![3]);
        assert!(v.parts[2].engine_active);
    }

    #[test]
    fn separation_conserves_mass_and_momentum() {
        let v = two_stage();
        let s = some_state();
        let mp = mass_properties(&v);
        let (a, b) = separate(&v, &s, 3, 2_500.0);
        assert_eq!(a.vessel.len(), 3);
        assert_eq!(b.vessel.len(), 7);
        let (ma, mb) = (a.vessel.total_mass(), b.vessel.total_mass());
        assert!((ma + mb - mp.mass).abs() < 1e-9);
        let p1 = a.state.vel * ma + b.state.vel * mb;
        assert!((p1 - s.vel * mp.mass).length() < 1e-6);
        let c1 = (a.state.pos * ma + b.state.pos * mb) / mp.mass;
        assert!((c1 - s.pos).length() < 1e-9);
        let rel = s.attitude.inverse().apply(b.state.vel - a.state.vel);
        assert!(rel.z() < 0.0);
    }

    #[test]
    fn force_free_step_is_exactly_kepler() {
        // Über der Atmosphäre, ohne Schub, ohne Drehung: ein Schritt beliebiger
        // Länge ist die exakte Keplerbahn.
        let body = planet(Atmosphere::earthlike());
        let mut v = two_stage();
        let r = 700_000.0;
        let mut s = RigidState {
            body: BodyId(0),
            time: 0.0,
            pos: Vec3::new(r, 0.0, 0.0),
            vel: Vec3::new(0.0, (body.mu / r).sqrt(), 0.0),
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(0.3)),
            omega: Vec3::ZERO,
            contact_force: 0.0,
        };
        let expect = kepler::propagate(s.orbit(), body.mu, 600.0);
        step(&mut v, &mut s, &Controls::default(), &body, 600.0);
        assert_eq!(s.pos, expect.r);
        assert_eq!(s.vel, expect.v);
    }

    #[test]
    fn vacuum_is_the_limit_of_thin_air() {
        // Gleicher Flug durch immer dünnere Luft konvergiert stetig gegen
        // den Flug im Vakuum, ohne Sprung.
        let fly = |scale: f64| {
            let mut atm = Atmosphere::earthlike();
            atm.surface_pressure *= scale;
            let body = planet(atm);
            let mut v = two_stage();
            v.activate_next_stage();
            let mut s = RigidState {
                body: BodyId(0),
                time: 0.0,
                pos: Vec3::new(600_000.0 + 2_000.0, 0.0, 0.0),
                vel: Vec3::new(300.0, 100.0, 0.0),
                attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
                omega: Vec3::ZERO,
                contact_force: 0.0,
            };
            let c = Controls { throttle: 1.0, rotation: Vec3::new(0.2, 0.0, 0.0) };
            for _ in 0..500 {
                step(&mut v, &mut s, &c, &body, DT_PHYSICS);
            }
            s.pos
        };
        // Bei vollem Druck fliegt die Rakete deutlich anders (Luftkräfte
        // drehen sie). Je dünner die Luft, desto näher am Vakuumflug, und
        // zwar linear im Druck: kein Sprung, kein Schwellwert.
        let vac = fly(0.0);
        let d0 = (fly(1.0) - vac).length();
        let d1 = (fly(1e-4) - vac).length();
        let d2 = (fly(1e-6) - vac).length();
        let d3 = (fly(1e-8) - vac).length();
        assert!(d0 > 100.0, "{d0}");
        assert!(d1 > 0.0 && d2 < d1 && d3 < d2, "{d1} {d2} {d3}");
        assert!((d2 / d1 - 1e-2).abs() < 1e-3, "{}", d2 / d1);
        assert!((d3 / d2 - 1e-2).abs() < 1e-3, "{}", d3 / d2);
    }

    #[test]
    fn falling_tank_reaches_terminal_velocity() {
        // Atmosphäre mit (fast) konstanter Dichte: sehr große Skalenhöhe
        let mut atm = Atmosphere::earthlike();
        atm.scale_height = 1.0e9;
        atm.top = 1.0e9;
        atm.lapse_rate = 0.0;
        let body = planet(atm);
        let lib = PartLibrary::builtin();
        let mut v = VesselBuilder::new("Tank", lib.part("tank")).build();
        let mut s = RigidState {
            body: BodyId(0),
            time: 0.0,
            pos: Vec3::new(600_000.0 + 50_000.0, 0.0, 0.0),
            vel: Vec3::ZERO,
            // Boden voran: Spitze (+Z) zeigt nach oben
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
            omega: Vec3::ZERO,
            contact_force: 0.0,
        };
        // Lage festhalten: ein kurzer Zylinder würde sonst zu taumeln
        // beginnen; geprüft wird hier nur das Gleichgewicht längs.
        let attitude = s.attitude;
        for _ in 0..8_000 {
            step(&mut v, &mut s, &Controls::default(), &body, DT_PHYSICS);
            s.attitude = attitude;
            s.omega = Vec3::ZERO;
        }
        // Endgeschwindigkeit: dort, wo der Widerstand des Modells das Gewicht trägt
        let air = atm.sample(s.altitude(&body));
        let g = body.mu / s.pos.length_squared();
        let m = v.total_mass();
        let drag = |speed: f64| crate::aero::vessel_forces(&v, Vec3::new(0.0, 0.0, speed), &air).0.length();
        let (mut lo, mut hi) = (1.0, 2_000.0);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if drag(mid) < m * g { lo = mid } else { hi = mid }
        }
        let vt = lo;
        let speed = s.vel.length();
        assert!((speed - vt).abs() / vt < 0.005, "v {speed:.2} m/s, erwartet {vt:.2} m/s");
        // und im Bereich der Messwerte für stumpfe Zylinder (Cd 0,8–1,1)
        let a = std::f64::consts::PI * 0.625 * 0.625;
        let cd = 2.0 * m * g / (air.density * vt * vt * a);
        assert!(cd > 0.8 && cd < 1.1, "Cd = {cd}");
    }

    #[test]
    fn control_surfaces_steer_in_commanded_direction() {
        let body = planet(Atmosphere::earthlike());
        let mut v = two_stage();
        let base = RigidState {
            body: BodyId(0),
            time: 0.0,
            pos: Vec3::new(600_000.0 + 3_000.0, 0.0, 0.0),
            vel: Vec3::new(300.0, 0.0, 0.0),
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
            omega: Vec3::ZERO,
            contact_force: 0.0,
        };
        // Kapsel ohne Reaktionsrad betrachten: nur Steuerflächen wirken
        let mut no_wheel = v.clone();
        let mut def = (*no_wheel.parts[0].def).clone();
        def.modules.clear();
        no_wheel.parts[0].def = std::sync::Arc::new(def);
        for cmd in [Vec3::X, -Vec3::X, Vec3::Y, Vec3::Z] {
            let mut s = base;
            step(&mut no_wheel, &mut s, &Controls { throttle: 0.0, rotation: cmd }, &body, DT_PHYSICS);
            assert!(s.omega.dot(cmd) > 0.0, "Kommando {cmd:?} → ω {:?}", s.omega);
        }
        let _ = &mut v;
    }

    #[test]
    fn corotating_air_exerts_no_force_on_resting_vessel() {
        let mut body = planet(Atmosphere::earthlike());
        body.rotation_rate = CelestialBody::rate_from_period(21_549.425);
        let mut v = two_stage();
        let site = GeoPosition { latitude: 0.5, longitude: 1.0, altitude: 3_000.0 };
        let mut s = RigidState::at_rest(&body, site, 1_234.0);
        assert!((s.vel.length() - body.rotation_rate * (body.radius + 3_000.0) * 0.5_f64.cos()).abs() < 1e-9);
        let r = step(&mut v, &mut s, &Controls::default(), &body, DT_PHYSICS);
        assert!(r.aero_force < 1e-9, "Luftkraft {} N", r.aero_force);
        assert_eq!(r.dynamic_pressure, 0.0);
        // Nach einem Schritt freien Falls: fast dieselbe Lage über Grund
        let g = s.geo_position(&body);
        assert!((g.latitude - 0.5).abs() < 1e-9 && (g.longitude - 1.0).abs() < 1e-9);
        assert!((s.time - 1_234.02).abs() < 1e-12);
    }

    #[test]
    fn rotation_gives_eastward_launch_bonus() {
        let mut body = planet(Atmosphere::VACUUM);
        body.rotation_rate = CelestialBody::rate_from_period(21_549.425);
        let s = RigidState::at_rest(&body, GeoPosition { latitude: 0.0, longitude: 0.0, altitude: 10.0 }, 0.0);
        // Am Äquator: ~175 m/s Bahngeschwindigkeit schon auf der Rampe
        assert!((s.vel.length() - 174.9).abs() < 0.2);
        assert!(s.surface_velocity(&body).length() < 1e-9);
    }

    /// Antriebslos mit 10° Anstellwinkel in 5 km Höhe, ohne Steuerung.
    /// Rückgabe: Anstellwinkel nach `seconds` Sekunden (Grad).
    fn coast_with_aoa(fins: bool, seconds: f64) -> f64 {
        let body = planet(Atmosphere::earthlike());
        let mut v = crate::scenario::first_rocket(&PartLibrary::builtin());
        if !fins {
            v.parts.retain(|p| p.def.control_surface().is_none());
            v.invalidate_geometry();
        }
        let up = Vec3::<Inertial>::new(1.0, 0.0, 0.0);
        let east = Vec3::<Inertial>::new(0.0, 1.0, 0.0);
        let aoa = 10f64.to_radians();
        // Spitze 45° über dem Horizont, Flugrichtung 10° tiefer
        let nose = (up + east).normalize_or_zero();
        let vel_dir = (up * (45f64 - 10.0).to_radians().sin() + east * (45f64 - 10.0).to_radians().cos()).normalize_or_zero();
        let mut s = RigidState {
            body: BodyId(0),
            time: 0.0,
            pos: up * (600_000.0 + 5_000.0),
            vel: vel_dir * 250.0,
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_arc(glam::DVec3::Z, nose.raw())),
            omega: Vec3::ZERO,
            contact_force: 0.0,
        };
        assert!((s.nose().dot(vel_dir).acos() - aoa).abs() < 1e-9);
        let n = (seconds / DT_PHYSICS) as usize;
        for _ in 0..n {
            step(&mut v, &mut s, &Controls::default(), &body, DT_PHYSICS);
        }
        let surf = s.surface_velocity(&body).normalize_or_zero();
        s.nose().dot(surf).clamp(-1.0, 1.0).acos().to_degrees()
    }

    #[test]
    fn rocket_with_fins_weathervanes_without_fins_flips() {
        let with = coast_with_aoa(true, 4.0);
        let without = coast_with_aoa(false, 4.0);
        assert!(with < 5.0, "mit Finnen: Anstellwinkel {with:.1}°");
        assert!(without > 90.0, "ohne Finnen: Anstellwinkel {without:.1}° (muss umschlagen)");
    }

    #[test]
    fn step_size_stops_at_atmosphere_entry() {
        let body = planet(Atmosphere::earthlike());
        // fällt von 200 km senkrecht
        let s = RigidState {
            body: BodyId(0),
            time: 0.0,
            pos: Vec3::new(800_000.0, 0.0, 0.0),
            vel: Vec3::ZERO,
            attitude: Rotation::IDENTITY,
            omega: Vec3::ZERO,
            contact_force: 0.0,
        };
        let dt = suggest_dt(1.0e6, &two_stage(), &s, None, &body);
        let at = kepler::propagate(s.orbit(), body.mu, dt).r.length();
        assert!((at - body.quiet_radius()).abs() < 1e-3, "endet bei {at}");
    }
}
