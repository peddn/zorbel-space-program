//! Headless-Testflug durch die Atmosphäre in die Umlaufbahn.
//!
//! Testplanet mit KSP-artigen Kenndaten (600 km Radius,
//! μ = 3,5316·10¹² m³/s², Atmosphäre bis 70 km), nicht rotierend.
//!
//! Der ganze Flug, auch der Zeitraffer und der Absturz der Unterstufe,
//! läuft über dieselbe Funktion `step`. Es gibt keinen Wechsel des
//! Physikmodells beim Verlassen oder Betreten der Atmosphäre; nur die
//! Schrittweite passt sich an die Stärke der Kräfte an.

use sim_core::kepler::{self, Elements};
use sim_core::{Atmosphere, CelestialBody, Inertial, Rotation, Vec3, VesselFrame};
use std::collections::BTreeMap;
use vessel::dynamics::{DT_PHYSICS, separate, step, suggest_dt};
use vessel::mass::mass_properties;
use vessel::{Controls, LoadKind, PartLibrary, RigidState, StepReport, Vessel, VesselBuilder};

const TARGET_ALT: f64 = 100_000.0;
/// Größter erlaubter Anstellwinkel beim Aufstieg (begrenzt Querlasten).
const MAX_AOA_DEG: f64 = 4.0;

fn build_rocket(lib: &PartLibrary) -> Vessel {
    let mut b = VesselBuilder::new("Erstflug", lib.part("kapsel"));
    let t2 = b.stack_below(0, lib.part("tank"));
    let e2 = b.stack_below(t2, lib.part("triebwerk"));
    let dec = b.stack_below(e2, lib.part("entkoppler"));
    let t1a = b.stack_below(dec, lib.part("tank"));
    let t1b = b.stack_below(t1a, lib.part("tank"));
    let t1c = b.stack_below(t1b, lib.part("tank"));
    let e1 = b.stack_below(t1c, lib.part("triebwerk"));
    b.radial_symmetric(t1c, lib.part("steuerflaeche"), 4, -0.5);
    // Höchste Stufe zuerst: 1 = Zündung unten, 0 = Trennen + Zündung oben
    b.set_stage(e1, 1).set_stage(dec, 0).set_stage(e2, 0);
    b.build()
}

/// Lageregler: dreht die Spitze in Richtung `target`, dämpft Rollen.
/// Die Verstärkung wird aus Trägheit und verfügbarem Stellmoment
/// (Reaktionsrad + Schubvektor) berechnet, damit der Regler bei jeder
/// Masse gleich schnell und stabil reagiert.
fn attitude_control(v: &Vessel, s: &RigidState, target: Vec3<Inertial>, thrust: f64) -> Vec3<VesselFrame> {
    const RESPONSE: f64 = 0.4; // s, Zeitkonstante der Drehratenregelung
    let nose = s.nose();
    let axis = nose.cross(target);
    let angle = axis.length().atan2(nose.dot(target));
    let axis_body = s.attitude.inverse().apply(axis.normalize_or_zero());
    let want = axis_body * (angle * 0.8).min(0.2);
    let err = Vec3::<VesselFrame>::new(want.x(), want.y(), 0.0) - s.omega;

    let mp = mass_properties(v);
    let wheel: f64 = v.parts.iter().map(|p| p.def.command_torque()).sum();
    let gimbal: f64 = v
        .parts
        .iter()
        .filter(|p| p.engine_active)
        .filter_map(|p| p.def.engine().map(|e| e.gimbal_rad.sin() * (p.origin.z() - mp.com.z()).abs()))
        .fold(0.0, f64::max)
        * thrust;
    let i = mp.inertia.0;
    let (pitch_yaw, roll) = ((wheel + gimbal).max(1.0), wheel.max(1.0));
    Vec3::new(
        err.x() * i.x_axis.x / (pitch_yaw * RESPONSE),
        err.y() * i.y_axis.y / (pitch_yaw * RESPONSE),
        err.z() * i.z_axis.z / (roll * RESPONSE),
    )
}

/// Dreht `from` höchstens um `max_angle` in Richtung `to`.
fn limit_towards(from: Vec3<Inertial>, to: Vec3<Inertial>, max_angle: f64) -> Vec3<Inertial> {
    let angle = from.dot(to).clamp(-1.0, 1.0).acos();
    if angle <= max_angle {
        return to;
    }
    let perp = to.reject_from_unit(from).normalize_or_zero();
    from * max_angle.cos() + perp * max_angle.sin()
}

fn angle_deg(a: Vec3<Inertial>, b: Vec3<Inertial>) -> f64 {
    a.normalize_or_zero().dot(b.normalize_or_zero()).clamp(-1.0, 1.0).acos().to_degrees()
}

fn fmt_alt(r: f64, body: &CelestialBody) -> String {
    format!("{:.1} km", (r - body.radius) / 1000.0)
}

fn orbit_line(el: &Elements, body: &CelestialBody) -> String {
    let ap = el.apoapsis.map(|a| fmt_alt(a, body)).unwrap_or_else(|| "offen".into());
    format!("Ap {ap}, Pe {}", fmt_alt(el.periapsis, body))
}

#[derive(Default)]
struct Log {
    /// Verbindung → (max. Auslastung, Lastart, Zeitpunkt)
    joints: BTreeMap<String, (f64, LoadKind, f64)>,
    max_g: f64,
    max_q: (f64, f64, f64), // (q, t, Höhe)
    steps: usize,
    steps_in_air: usize,
}

impl Log {
    fn record(&mut self, v: &Vessel, r: &StepReport, t: f64, alt: f64, g0: f64) {
        self.steps += 1;
        if r.air.density > 0.0 {
            self.steps_in_air += 1;
        }
        self.max_g = self.max_g.max(r.proper_accel / g0);
        if r.dynamic_pressure > self.max_q.0 {
            self.max_q = (r.dynamic_pressure, t, alt);
        }
        for j in &r.loads.joints {
            let e = self.joints.entry(v.joint_label(j.child)).or_insert((0.0, j.critical, t));
            if j.utilization > e.0 {
                *e = (j.utilization, j.critical, t);
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Phase {
    Ascent,
    Coast,
    Circularize,
}

fn fail(msg: String) -> ! {
    println!("\n!! {msg}");
    std::process::exit(1);
}

fn main() {
    let body = CelestialBody {
        name: "Testplanet".into(),
        mu: 3.5316e12,
        radius: 600_000.0,
        soi_radius: 8.4e7,
        atmosphere: Atmosphere::earthlike(),
    };
    let top = body.atmosphere.top();
    let lib = PartLibrary::builtin();
    let mut ship = build_rocket(&lib);
    let g0 = body.surface_gravity();

    // Auf dem Startplatz: Spitze zeigt radial nach außen (+X).
    let mp = mass_properties(&ship);
    let bottom = ship.parts.iter().map(|p| p.origin.z() - p.def.shape.stack_length() / 2.0).fold(f64::INFINITY, f64::min);
    let mut s = RigidState {
        pos: Vec3::new(body.radius + (mp.com.z() - bottom) + 0.1, 0.0, 0.0),
        vel: Vec3::ZERO,
        attitude: Rotation::from_quat(glam::DQuat::from_rotation_y(std::f64::consts::FRAC_PI_2)),
        omega: Vec3::ZERO,
    };
    let engine = lib.part("triebwerk").engine().unwrap();
    let sl_thrust = engine.thrust(1.0, body.atmosphere.sample(0.0).pressure);

    println!("== {} auf {} (Atmosphäre bis {:.0} km) ==", ship.name, body.name, top / 1000.0);
    println!("Teile: {}, Startmasse {:.0} kg, Schub am Boden {:.0} kN, Schub/Gewicht {:.2}", ship.len(), mp.mass, sl_thrust / 1000.0, sl_thrust / (mp.mass * g0));
    println!();

    let mut log = Log::default();
    let mut t = 0.0;
    let mut debris: Vec<(f64, RigidState, Vessel)> = Vec::new();
    let mut last: Option<StepReport> = None;
    let mut phase = Phase::Ascent;
    let mut in_air = true;
    let mut next_telemetry = 0.0;

    ship.activate_next_stage();
    println!("T+{t:>6.1}s  Zündung Stufe 1");
    println!("           {:>8} {:>8} {:>6} {:>8} {:>6} {:>6}", "Höhe", "v", "Mach", "q", "Neig.", "AoA");

    loop {
        let alt = s.altitude(&body);
        let up = s.pos.normalize_or_zero();
        let east = Vec3::<Inertial>::Z.cross(up).normalize_or_zero();
        let speed = s.vel.length();
        let prograde = if speed > 1.0 { s.vel / speed } else { up };
        let el = kepler::elements(s.orbit(), body.mu);
        let apo_alt = el.apoapsis.map_or(f64::INFINITY, |a| a - body.radius);

        // --- Führung ---
        let (requested, throttle, target) = match phase {
            Phase::Ascent => {
                if apo_alt >= TARGET_ALT {
                    println!("T+{t:>6.1}s  Brennschluss bei {:.1} km, {:.0} m/s, {}", alt / 1000.0, speed, orbit_line(&el, &body));
                    phase = Phase::Coast;
                    continue;
                }
                let pitch: f64 = if speed < 60.0 { 90.0 } else { 90.0 * (1.0 - (alt / 48_000.0).min(1.0).sqrt()) };
                let program = up * pitch.to_radians().sin() + east * pitch.to_radians().cos();
                let target = if speed < 60.0 { up } else { limit_towards(prograde, program, MAX_AOA_DEG.to_radians()) };
                (DT_PHYSICS, 1.0, target)
            }
            Phase::Coast => {
                // Luftwiderstand senkt die Apoapsis noch; nachschieben.
                let t_apo = el.time_to_apoapsis.unwrap_or(0.0);
                let r_apo = el.apoapsis.unwrap_or(s.pos.length());
                let v_apo = (body.mu * (2.0 / r_apo - 1.0 / el.semi_major_axis)).sqrt();
                let dv = (body.mu / r_apo).sqrt() - v_apo;
                let burn = dv * ship.total_mass() / engine.thrust_vac;
                let burn_start_lead = burn / 2.0;
                if t_apo <= burn_start_lead {
                    println!("T+{t:>6.1}s  Zündung Zirkularisierung bei {:.1} km, geplant Δv {dv:.0} m/s ({burn:.1} s)", alt / 1000.0);
                    phase = Phase::Circularize;
                    continue;
                }
                let topup = alt < top && apo_alt < TARGET_ALT - 200.0;
                let horizontal = s.vel.reject_from_unit(up).normalize_or_zero();
                let requested = if t_apo - burn_start_lead > 60.0 { t_apo - burn_start_lead - 50.0 } else { DT_PHYSICS };
                (if topup { DT_PHYSICS } else { requested }, if topup { 0.3 } else { 0.0 }, if topup { prograde } else { horizontal })
            }
            Phase::Circularize => {
                let horizontal = s.vel.reject_from_unit(up);
                let dv_left = (body.mu / s.pos.length()).sqrt() - horizontal.length();
                let empty = ship.parts.iter().all(|p| p.def.fuel_capacity() == 0.0 || p.fuel <= 0.0);
                if dv_left <= 0.05 || empty {
                    println!("T+{t:>6.1}s  Brennschluss Zirkularisierung, {}", orbit_line(&el, &body));
                    break;
                }
                let throttle = (dv_left * ship.total_mass() / engine.thrust_vac / 2.0).clamp(0.02, 1.0);
                (DT_PHYSICS, throttle, horizontal.normalize_or_zero())
            }
        };

        // --- Ein Schritt, immer dieselbe Physik ---
        let dt = suggest_dt(requested, &s, last.as_ref(), &body);
        // Im Zeitraffer regelt der Autopilot nicht (Lage wird gehalten).
        let thrust_now = last.as_ref().map_or(0.0, |r| r.thrust);
        let rotation = if dt <= 0.1 { attitude_control(&ship, &s, target, thrust_now.max(throttle * engine.thrust_vac * 0.5)) } else { Vec3::ZERO };
        if dt > 1.0 {
            println!("T+{t:>6.1}s  Zeitraffer: ein Schritt über {dt:.1} s");
        }
        let report = step(&mut ship, &mut s, &Controls { throttle, rotation }, &body, dt);
        t += dt;
        let alt = s.altitude(&body);
        log.record(&ship, &report, t, alt, g0);

        if phase == Phase::Ascent && t >= next_telemetry {
            println!(
                "T+{t:>6.1}s  {:>6.1} km {:>6.0} m/s {:>6.2} {:>5.1} kPa {:>5.1}° {:>5.1}°",
                alt / 1000.0,
                s.vel.length(),
                report.mach,
                report.dynamic_pressure / 1000.0,
                90.0 - angle_deg(s.nose(), up),
                angle_deg(s.nose(), s.vel)
            );
            next_telemetry += 10.0;
        }
        if in_air && alt > top {
            println!("T+{t:>6.1}s  Atmosphäre verlassen ({:.0} km) – Physik läuft unverändert weiter", top / 1000.0);
            in_air = false;
        }
        if let Some(j) = report.loads.broken().next() {
            fail(format!("Strukturversagen T+{t:.1}s: {} ({})", ship.joint_label(j.child), j.critical.label()));
        }
        if alt < 0.0 {
            fail(format!("Aufschlag T+{t:.1}s"));
        }

        if report.flameout && ship.next_stage().is_some() {
            println!("T+{t:>6.1}s  Brennschluss Stufe 1 bei {:.1} km, {:.0} m/s", alt / 1000.0, s.vel.length());
            for d in ship.activate_next_stage().into_iter().rev() {
                let impulse = ship.parts[d].def.decoupler_impulse().unwrap_or(0.0);
                let (rest, lower) = separate(&ship, &s, d, impulse);
                println!("T+{t:>6.1}s  Stufentrennung: {} Teile bleiben, {} abgeworfen ({:.0} kg)", rest.vessel.len(), lower.vessel.len(), lower.vessel.total_mass());
                debris.push((t, lower.state, lower.vessel));
                ship = rest.vessel;
                s = rest.state;
            }
            println!("T+{t:>6.1}s  Zündung Stufe 2");
        }
        last = Some(report);
    }

    // --- Ergebnis ---
    let el = kepler::elements(s.orbit(), body.mu);
    let fuel_left: f64 = ship.parts.iter().map(|p| p.fuel).sum();
    println!();
    println!("== Endbahn ==");
    println!("{}, Exzentrizität {:.5}, Umlaufzeit {:.1} min, Resttreibstoff {fuel_left:.0} kg", orbit_line(&el, &body), el.eccentricity, el.period.unwrap_or(f64::NAN) / 60.0);

    // Zehn Umläufe mit derselben step-Funktion (ohne Schub, über der
    // Atmosphäre ergibt das große, exakte Schritte).
    let mut check = s;
    let mut v_check = ship.clone();
    let (mut tt, mut n) = (0.0, 0);
    let horizon = 10.0 * el.period.unwrap();
    let mut prev = None;
    while tt < horizon {
        let dt = suggest_dt(horizon - tt, &check, prev.as_ref(), &body);
        prev = Some(step(&mut v_check, &mut check, &Controls::default(), &body, dt));
        tt += dt;
        n += 1;
    }
    println!("Nach 10 Umläufen ({n} Schritte): {}", orbit_line(&kepler::elements(check.orbit(), body.mu), &body));

    // --- Unterstufe: ungesteuert zurück durch die Atmosphäre ---
    println!();
    println!("== Abgeworfene Teile ==");
    let mut queue = debris;
    while let Some((t0, mut ds, mut dv)) = queue.pop() {
        let mut t = t0;
        let mut prev: Option<StepReport> = None;
        let (mut apogee, mut max_g, mut max_q, mut max_w, mut steps) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0);
        let mut max_load: (f64, String) = (0.0, String::new());
        loop {
            let dt = suggest_dt(600.0, &ds, prev.as_ref(), &body);
            let r = step(&mut dv, &mut ds, &Controls::default(), &body, dt);
            t += dt;
            steps += 1;
            let alt = ds.altitude(&body);
            apogee = apogee.max(alt);
            max_g = max_g.max(r.proper_accel / g0);
            max_q = max_q.max(r.dynamic_pressure);
            max_w = max_w.max(ds.omega.length());
            if let Some(j) = r.loads.max_utilization()
                && j.utilization > max_load.0
            {
                max_load = (j.utilization, format!("{} ({})", dv.joint_label(j.child), j.critical.label()));
            }
            if let Some(j) = r.loads.broken().next() {
                println!("{}: bricht T+{t:.1}s in {:.1} km Höhe an {}", dv.name, alt / 1000.0, dv.joint_label(j.child));
                let (a, b) = separate(&dv, &ds, j.child, 0.0);
                queue.push((t, b.state, b.vessel));
                dv = a.vessel;
                ds = a.state;
            }
            if alt <= 0.0 || steps > 5_000_000 {
                let end = if alt <= 0.0 { format!("Aufschlag T+{t:.1}s mit {:.0} m/s", ds.vel.length()) } else { "Abbruch: Schrittlimit".into() };
                println!(
                    "{} ({} Teile): Gipfel {:.1} km, max. {:.0} kPa / {:.1} g / {:.1} rad/s (taumelt), {end}, {steps} Schritte",
                    dv.name,
                    dv.len(),
                    apogee / 1000.0,
                    max_q / 1000.0,
                    max_g,
                    max_w
                );
                println!("    höchste Auslastung: {:.0} % {}", max_load.0 * 100.0, max_load.1);
                break;
            }
            prev = Some(r);
        }
    }

    println!();
    println!("== Aufstieg ==");
    println!("Max Q: {:.1} kPa bei T+{:.1}s in {:.1} km Höhe", log.max_q.0 / 1000.0, log.max_q.1, log.max_q.2 / 1000.0);
    println!("max. Eigenbeschleunigung: {:.2} g", log.max_g);
    println!("Schritte: {} (davon {} mit Luftkontakt)", log.steps, log.steps_in_air);
    println!();
    println!("== Höchste Auslastung je Verbindung ==");
    for (label, (u, kind, at)) in &log.joints {
        println!("{label:<32} {:>6.1} %  {:<9} T+{at:>6.1}s", u * 100.0, kind.label());
    }
}
