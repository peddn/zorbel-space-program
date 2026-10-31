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

use crate::aero::{PartAero, aero_geometry, part_force};
use crate::loads::{ExternalLoad, LoadAnalysis, Motion, joint_loads, motion_from_loads};
use crate::mass::{MassProperties, mass_properties};
use crate::vessel::Vessel;
use sim_core::kepler::{self, StateVector};
use sim_core::{AtmosphereSample, CelestialBody, Inertial, Rotation, Vec3, VesselFrame};

/// Kleinster Schritt bei voller Dynamik (Schub, dichte Luft).
pub const DT_PHYSICS: f64 = 0.02;

#[derive(Debug, Clone, Copy)]
pub struct RigidState {
    /// Position des Schwerpunkts relativ zum Zentralkörper.
    pub pos: Vec3<Inertial>,
    pub vel: Vec3<Inertial>,
    pub attitude: Rotation<VesselFrame, Inertial>,
    /// Winkelgeschwindigkeit im Schiffssystem.
    pub omega: Vec3<VesselFrame>,
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

fn evaluate(v: &Vessel, geo: &[PartAero], mp: &MassProperties, s: &RigidState, c: &Controls, body: &CelestialBody) -> Evaluation {
    let air = body.atmosphere.sample(s.altitude(body));
    let cmd = Vec3::<VesselFrame>::new(
        c.rotation.x().clamp(-1.0, 1.0),
        c.rotation.y().clamp(-1.0, 1.0),
        c.rotation.z().clamp(-1.0, 1.0),
    );
    // Luft ruht im Inertialsystem (Planet rotiert noch nicht).
    let v_body = s.attitude.inverse().apply(s.vel);

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

        // Luftkraft am eigenen Ort des Teils
        let point = p.com();
        let lever = point - mp.com;
        let v_part = v_body + s.omega.cross(lever);
        let v_local = p.orientation.inverse().apply(v_part);
        let force_at = |defl: f64| p.orientation.apply(part_force(p, &geo[i], v_local, &air, defl));
        let deflection = match p.def.control_surface() {
            Some(cs) if cmd.length_squared() > 0.0 => {
                // Welche Drehung erzeugt ein positiver Ausschlag gerade?
                let f0 = force_at(0.0);
                let d_torque = lever.cross(force_at(DEFLECTION_PROBE) - f0);
                let len = d_torque.length();
                if len > 0.0 { cs.max_deflection_rad * cmd.dot(d_torque / len).clamp(-1.0, 1.0) } else { 0.0 }
            }
            _ => 0.0,
        };
        let f = force_at(deflection);
        if f != Vec3::ZERO {
            aero_total += f;
            ext.push(ExternalLoad { part: i, force: f, point, torque: Vec3::ZERO });
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
    let geo = aero_geometry(v);
    let mp0 = mass_properties(v);
    let ev = evaluate(v, &geo, &mp0, s, c, body);
    let loads = joint_loads(v, &mp0, &ev.motion, &ev.ext);

    // Kick
    s.vel += s.attitude.apply(ev.motion.accel) * (dt / 2.0);
    s.omega += ev.motion.alpha * (dt / 2.0);

    // Drift: exakte Keplerbahn, freie Drehung
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

    let q = 0.5 * ev.air.density * ev.airspeed * ev.airspeed;
    StepReport {
        loads,
        thrust: ev.thrust,
        aero_force: ev.aero_force,
        proper_accel: ev.motion.accel.length(),
        mach: ev.airspeed / ev.air.speed_of_sound,
        dynamic_pressure: q,
        air: ev.air,
        flameout: ev.flameout,
    }
}

/// Schrittweite für den nächsten Schritt, höchstens `requested`.
///
/// - Starke Störkräfte (Schub, dichte Luft) → kurze Schritte, die Länge
///   wächst stetig, je schwächer sie werden.
/// - Schnelle Drehung → kurze Schritte (höchstens ~1° pro Schritt).
/// - Oberhalb der Atmosphäre: Schritt endet spätestens beim Eintritt, damit
///   kein langer Kepler-Schritt durch die Lufthülle "hindurchspringt".
pub fn suggest_dt(requested: f64, s: &RigidState, last: Option<&StepReport>, body: &CelestialBody) -> f64 {
    let mut dt = requested;
    if let Some(r) = last
        && r.proper_accel > 0.0
    {
        dt = dt.min(DT_PHYSICS * (1.0 / r.proper_accel).clamp(1.0, 500.0));
    }
    let w = s.omega.length();
    if w > 0.0 {
        dt = dt.min((0.02 / w).max(DT_PHYSICS));
    }
    let quiet = body.quiet_radius();
    if s.pos.length() > quiet
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
            pos: s.pos + s.attitude.apply(d),
            vel: s.vel + s.attitude.apply(s.omega.cross(d) + axis * (push / mpp.mass)),
            attitude: s.attitude,
            omega: s.omega,
        };
        Separated { vessel, state }
    };
    (piece(rest, -impulse), piece(detached, impulse))
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
            pos: Vec3::new(700_000.0, 0.0, 0.0),
            vel: Vec3::new(0.0, 2_200.0, 0.0),
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(0.3)),
            omega: Vec3::new(0.05, 0.0, 0.2),
        }
    }

    fn planet(atmosphere: Atmosphere) -> CelestialBody {
        CelestialBody { name: "Test".into(), mu: 3.5316e12, radius: 600_000.0, soi_radius: 8.4e7, atmosphere }
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
            pos: Vec3::new(r, 0.0, 0.0),
            vel: Vec3::new(0.0, (body.mu / r).sqrt(), 0.0),
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(0.3)),
            omega: Vec3::ZERO,
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
                pos: Vec3::new(600_000.0 + 2_000.0, 0.0, 0.0),
                vel: Vec3::new(300.0, 100.0, 0.0),
                attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
                omega: Vec3::ZERO,
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
            pos: Vec3::new(600_000.0 + 50_000.0, 0.0, 0.0),
            vel: Vec3::ZERO,
            // Boden voran: Spitze (+Z) zeigt nach oben
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
            omega: Vec3::ZERO,
        };
        for _ in 0..8_000 {
            step(&mut v, &mut s, &Controls::default(), &body, DT_PHYSICS);
        }
        let air = atm.sample(0.0);
        let g = body.mu / s.pos.length_squared();
        let cda = (CD_FLAT_T + CD_BASE_T) * std::f64::consts::PI * 0.625 * 0.625;
        let m = v.total_mass();
        // v_t = √(2mg / (ρ·CdA·f(M))), Machfaktor selbstkonsistent
        let mut vt = 100.0;
        for _ in 0..50 {
            let f = crate::aero::mach_factor(vt / air.speed_of_sound);
            vt = (2.0 * m * g / (air.density * cda * f)).sqrt();
        }
        let speed = s.vel.length();
        assert!((speed - vt).abs() / vt < 0.005, "v {speed:.2} m/s, erwartet {vt:.2} m/s");
    }
    const CD_FLAT_T: f64 = crate::aero::CD_FLAT;
    const CD_BASE_T: f64 = crate::aero::CD_BASE;

    #[test]
    fn control_surfaces_steer_in_commanded_direction() {
        let body = planet(Atmosphere::earthlike());
        let mut v = two_stage();
        let base = RigidState {
            pos: Vec3::new(600_000.0 + 3_000.0, 0.0, 0.0),
            vel: Vec3::new(300.0, 0.0, 0.0),
            attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
            omega: Vec3::ZERO,
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
    fn step_size_stops_at_atmosphere_entry() {
        let body = planet(Atmosphere::earthlike());
        // fällt von 200 km senkrecht
        let s = RigidState {
            pos: Vec3::new(800_000.0, 0.0, 0.0),
            vel: Vec3::ZERO,
            attitude: Rotation::IDENTITY,
            omega: Vec3::ZERO,
        };
        let dt = suggest_dt(1.0e6, &s, None, &body);
        let at = kepler::propagate(s.orbit(), body.mu, dt).r.length();
        assert!((at - body.quiet_radius()).abs() < 1e-3, "endet bei {at}");
    }
}
