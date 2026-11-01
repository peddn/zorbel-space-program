//! Headless-Testflug durch die Atmosphäre in die Umlaufbahn.
//!
//! Testplanet mit KSP-artigen Kenndaten (600 km Radius,
//! μ = 3,5316·10¹² m³/s², Atmosphäre bis 70 km), rotierend.
//!
//! Der ganze Flug, auch der Zeitraffer und der Absturz der Unterstufe,
//! läuft über dieselbe Funktion `step`. Es gibt keinen Wechsel des
//! Physikmodells beim Verlassen oder Betreten der Atmosphäre; nur die
//! Schrittweite passt sich an die Stärke der Kräfte an. Geflogen wird vom
//! Aufstiegs-Autopiloten aus `vessel::autopilot`.

use sim_core::kepler::{self, Elements};
use sim_core::{CelestialBody, Inertial, Vec3};
use std::collections::BTreeMap;
use vessel::autopilot::{AscentAutopilot, Event, Phase};
use vessel::dynamics::{activate_stage, separate, step, suggest_dt};
use vessel::mass::mass_properties;
use vessel::scenario::{LAUNCH_SITE, first_rocket, on_launch_site, test_planet};
use vessel::{Controls, LoadKind, PartLibrary, RigidState, StepReport, Vessel};

const TARGET_ALT: f64 = 100_000.0;

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

fn fail(msg: String) -> ! {
    println!("\n!! {msg}");
    std::process::exit(1);
}

fn main() {
    let body = test_planet();
    let top = body.atmosphere.top();
    let lib = PartLibrary::builtin();
    let mut ship = first_rocket(&lib);
    let g0 = body.surface_gravity();
    let pad = LAUNCH_SITE;

    // Startplatz am Äquator; die Rakete steht und dreht mit dem Planeten mit.
    let mp = mass_properties(&ship);
    let mut s = on_launch_site(&body, &ship, 0.0);
    let engine = lib.part("triebwerk").engine().unwrap();
    let sl_thrust = engine.thrust(1.0, body.atmosphere.sample(0.0).pressure);

    println!("== {} auf {} (Atmosphäre bis {:.0} km) ==", ship.name, body.name, top / 1000.0);
    println!("Teile: {}, Startmasse {:.0} kg, Schub am Boden {:.0} kN, Schub/Gewicht {:.2}", ship.len(), mp.mass, sl_thrust / 1000.0, sl_thrust / (mp.mass * g0));
    println!("Startplatz am Äquator: {:.1} m/s Bahngeschwindigkeit durch die Planetenrotation", s.vel.length());
    println!();

    let mut log = Log::default();
    let mut debris: Vec<(f64, RigidState, Vessel)> = Vec::new();
    let mut last: Option<StepReport> = None;
    let mut pilot = AscentAutopilot::new(TARGET_ALT);
    let mut in_air = true;
    let mut next_telemetry = 0.0;

    activate_stage(&mut ship, &mut s);
    println!("T+{:>6.1}s  Zündung Stufe 1", s.time);
    println!("           {:>9} {:>10} {:>10} {:>6} {:>9} {:>6} {:>6}", "Höhe", "v Oberfl.", "v Bahn", "Mach", "q", "Neig.", "AoA");

    while pilot.phase() != Phase::Done {
        let (plan, event) = pilot.update(&ship, &s, &body, last.as_ref());
        let t = s.time;
        let alt = s.altitude(&body);
        match event {
            Some(Event::AscentCutoff { elements }) => {
                println!("T+{t:>6.1}s  Brennschluss bei {:.1} km, {:.0} m/s, {}", alt / 1000.0, s.vel.length(), orbit_line(&elements, &body))
            }
            Some(Event::CircularizationStart { delta_v, burn_time }) => println!(
                "T+{t:>6.1}s  Zündung Zirkularisierung bei {:.1} km, geplant Δv {delta_v:.0} m/s ({burn_time:.1} s)",
                alt / 1000.0
            ),
            Some(Event::Done { elements }) => println!("T+{t:>6.1}s  Brennschluss Zirkularisierung, {}", orbit_line(&elements, &body)),
            None => {}
        }
        if plan.stage {
            println!("T+{t:>6.1}s  Brennschluss Stufe 1 bei {:.1} km, {:.0} m/s", alt / 1000.0, s.vel.length());
            for lower in activate_stage(&mut ship, &mut s) {
                println!("T+{t:>6.1}s  Stufentrennung: {} Teile bleiben, {} abgeworfen ({:.0} kg)", ship.len(), lower.vessel.len(), lower.vessel.total_mass());
                debris.push((t, lower.state, lower.vessel));
            }
            println!("T+{t:>6.1}s  Zündung Stufe 2");
        }
        if event.is_some() || plan.stage {
            continue;
        }

        // --- Ein Schritt, immer dieselbe Physik ---
        let dt = suggest_dt(plan.requested_dt, &ship, &s, last.as_ref(), &body);
        let controls: Controls = plan.controls(&ship, &s, dt, last.as_ref().map_or(0.0, |r| r.thrust));
        if dt > 1.0 {
            println!("T+{t:>6.1}s  Zeitraffer: ein Schritt über {dt:.1} s");
        }
        let report = step(&mut ship, &mut s, &controls, &body, dt);
        let t = s.time;
        let alt = s.altitude(&body);
        log.record(&ship, &report, t, alt, g0);

        if pilot.phase() == Phase::Ascent && t >= next_telemetry {
            let surf = s.surface_velocity(&body);
            let up = s.pos.normalize_or_zero();
            println!(
                "T+{t:>6.1}s  {:>6.1} km {:>6.0} m/s {:>6.0} m/s {:>6.2} {:>5.1} kPa {:>5.1}° {:>5.1}°",
                alt / 1000.0,
                surf.length(),
                s.vel.length(),
                report.mach,
                report.dynamic_pressure / 1000.0,
                90.0 - angle_deg(s.nose(), up),
                angle_deg(s.nose(), surf)
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
        if let Some((part, speed)) = report.crash {
            fail(format!("Aufschlag T+{t:.1}s mit {speed:.1} m/s ({})", ship.parts[part].def.name));
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
        let dt = suggest_dt(horizon - tt, &v_check, &check, prev.as_ref(), &body);
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
            let dt = suggest_dt(600.0, &dv, &ds, prev.as_ref(), &body);
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
            let crashed = r.crash;
            let resting = r.touching_ground && ds.surface_velocity(&body).length() < 0.05;
            if crashed.is_none() && !resting && let Some(j) = r.loads.broken().next() {
                println!("{}: bricht T+{t:.1}s in {:.1} km Höhe an {}", dv.name, alt / 1000.0, dv.joint_label(j.child));
                let (a, b) = separate(&dv, &ds, j.child, 0.0);
                queue.push((t, b.state, b.vessel));
                dv = a.vessel;
                ds = a.state;
            }
            if crashed.is_some() || resting || steps > 5_000_000 {
                // Ort im mitdrehenden System, relativ zum Startplatz
                let g = ds.geo_position(&body);
                let place = format!(
                    "{:.1} km {} vom Startplatz",
                    body.surface_distance(pad, g) / 1000.0,
                    if g.longitude >= 0.0 { "östlich" } else { "westlich" }
                );
                let end = match crashed {
                    Some((part, speed)) => format!("zerschellt T+{t:.1}s mit {speed:.0} m/s ({}), {place}", dv.parts[part].def.name),
                    None if resting => format!("liegt T+{t:.1}s, {place}"),
                    None => "Abbruch: Schrittlimit".into(),
                };
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
