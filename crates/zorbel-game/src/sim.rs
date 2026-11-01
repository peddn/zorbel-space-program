//! Spielzustand: alle Objekte im Flug, Steuerung, Zeitraffer.
//!
//! Reine Logik ohne Darstellung. Jedes Objekt (Schiff oder Trümmer) wird
//! mit `vessel::dynamics::step` fortgeschrieben, genau wie im Testflug.

use bevy::prelude::Resource;
use sim_core::kepler::{self, Elements};
use sim_core::{CelestialBody, Inertial, Vec3, VesselFrame};
use vessel::autopilot::{AscentAutopilot, Event, Phase, attitude_control, rate_control};
use vessel::dynamics::{activate_stage, separate, step, suggest_dt};
use vessel::scenario::{first_rocket_with, on_launch_site, test_planet};
use vessel::{Controls, PartLibrary, RigidState, StepReport, Vessel};

static NEXT_GENERATION: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Zeitraffer-Stufen. Ab 10× regelt niemand mehr die Lage (Schritte zu lang).
pub const WARP_LEVELS: [f64; 8] = [1.0, 2.0, 4.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0];
/// Höchstens so viele Physikschritte pro Objekt und Bild.
const MAX_STEPS_PER_FRAME: usize = 4_000;
/// Größte Drehrate bei voller Handsteuerung (rad/s): Nicken/Gieren, Rollen.
const MANUAL_RATE: f64 = 0.26; // ≈ 15°/s
const MANUAL_ROLL_RATE: f64 = 0.52; // ≈ 30°/s

pub struct SimObject {
    /// Eindeutig über die ganze Sitzung; die Darstellung hängt daran.
    pub id: u64,
    pub vessel: Vessel,
    pub state: RigidState,
    pub last: Option<StepReport>,
    /// Zählt Änderungen an der Teileliste (Neuaufbau der Meshes).
    pub version: u64,
    /// Beim Aufprall zerstört: Teil und Geschwindigkeit.
    pub crashed: Option<(String, f64)>,
}

impl SimObject {
    pub fn elements(&self, body: &CelestialBody) -> Elements {
        kepler::elements(self.state.orbit(), body.mu)
    }
}

/// Eingaben des Spielers für das aktive Schiff.
#[derive(Default, Clone, Copy)]
pub struct PilotInput {
    pub throttle: f64,
    /// Drehkommando −1..1 je Achse (Nicken, Gieren, Rollen).
    pub rotation: Vec3<VesselFrame>,
    /// Lage halten, wenn nicht gesteuert wird.
    pub sas: bool,
}

#[derive(Resource)]
pub struct Sim {
    pub body: CelestialBody,
    /// Index 0 ist immer das aktive Schiff (Kamerafokus).
    pub objects: Vec<SimObject>,
    pub input: PilotInput,
    pub autopilot: Option<AscentAutopilot>,
    /// Steht noch auf der Rampe, gehalten von Startklammern.
    pub clamped: bool,
    pub warp_index: usize,
    /// Zeitraffer konnte im letzten Bild nicht voll erreicht werden.
    pub warp_limited: bool,
    /// Gehaltene Richtung für den Lagehalter.
    sas_target: Option<Vec3<Inertial>>,
    next_id: u64,
    pub messages: Vec<(f64, String)>,
    pub ship_lost: bool,
    /// Rakete mit Finnen (sonst die Variante ohne).
    pub fins: bool,
    /// Eindeutig pro Sitzung (Neustart = neue Generation). Die Darstellung
    /// erkennt daran, dass ihre Ansichten zu einer alten Sitzung gehören.
    pub generation: u64,
    /// Zuletzt tatsächlich gegebenes Gas (Spieler oder Autopilot).
    pub applied_throttle: f64,
}

impl Sim {
    pub fn new() -> Self {
        Self::with_fins(true)
    }

    pub fn with_fins(fins: bool) -> Self {
        let body = test_planet();
        let lib = PartLibrary::builtin();
        let vessel = first_rocket_with(&lib, fins);
        let state = on_launch_site(&body, &vessel, 0.0);
        Self {
            body,
            objects: vec![SimObject { id: 0, vessel, state, last: None, version: 0, crashed: None }],
            input: PilotInput { sas: true, ..Default::default() },
            autopilot: None,
            clamped: true,
            warp_index: 0,
            warp_limited: false,
            sas_target: None,
            next_id: 1,
            messages: vec![(0.0, "Bereit. Leertaste: Start".into())],
            ship_lost: false,
            fins,
            generation: NEXT_GENERATION.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            applied_throttle: 0.0,
        }
    }

    pub fn ship(&self) -> &SimObject {
        &self.objects[0]
    }

    pub fn time(&self) -> f64 {
        self.objects[0].state.time
    }

    pub fn warp(&self) -> f64 {
        WARP_LEVELS[self.warp_index]
    }

    fn log(&mut self, text: impl Into<String>) {
        let t = self.time();
        self.messages.push((t, text.into()));
        if self.messages.len() > 6 {
            self.messages.remove(0);
        }
    }

    /// Nächste Stufe; die erste löst die Startklammern.
    pub fn stage(&mut self) {
        if self.ship_lost {
            return;
        }
        if self.clamped {
            self.clamped = false;
            self.log("Startklammern gelöst");
        }
        let ship = &mut self.objects[0];
        if ship.vessel.next_stage().is_none() {
            return;
        }
        let dropped = activate_stage(&mut ship.vessel, &mut ship.state);
        ship.version += 1;
        // Alter Bericht (z. B. "Brennschluss") gilt nicht für die neue Stufe.
        ship.last = None;
        let n = dropped.len();
        for d in dropped {
            let id = self.next_id;
            self.next_id += 1;
            self.objects.push(SimObject { id, vessel: d.vessel, state: d.state, last: None, version: 0, crashed: None });
        }
        self.log(if n > 0 { "Stufentrennung" } else { "Zündung" });
    }

    /// Spieler greift ein (Lage oder Gas): Autopilot aus, aktuelles Gas bleibt.
    pub fn take_manual_control(&mut self) {
        if self.autopilot.take().is_some() {
            self.input.throttle = self.applied_throttle;
            self.log("Autopilot aus (Handsteuerung)");
        }
    }

    pub fn toggle_autopilot(&mut self) {
        if self.autopilot.take().is_some() {
            self.log("Autopilot aus");
        } else {
            self.autopilot = Some(AscentAutopilot::new(100_000.0));
            self.log("Autopilot: Umlaufbahn 100 km");
            if self.clamped {
                self.stage();
            }
        }
    }

    /// Alle Objekte bis `real_dt × Zeitraffer` weiter.
    pub fn advance(&mut self, real_dt: f64) {
        if self.ship_lost {
            return;
        }
        let target = self.time() + real_dt * self.warp();
        self.advance_to(target);
    }

    /// Alle Objekte bis zur Simulationszeit `target`.
    pub fn advance_to(&mut self, target: f64) {
        self.warp_limited = false;
        if self.clamped {
            // Auf der Rampe: Startklammern halten fest, Schiff dreht mit dem Planeten.
            let ship = &mut self.objects[0];
            ship.state = on_launch_site(&self.body, &ship.vessel, target);
            self.applied_throttle = self.input.throttle;
            return;
        }
        let mut i = 0;
        while i < self.objects.len() {
            self.advance_object(i, target);
            i += 1;
        }
        // Beim Aufprall zerstörte Objekte
        let mut msgs = Vec::new();
        if let Some((part, speed)) = self.objects[0].crashed.clone() {
            self.ship_lost = true;
            msgs.push(format!("Aufschlag mit {speed:.1} m/s ({part}) – Schiff zerstört (R: Neustart)"));
        }
        self.objects.retain(|o| {
            if o.id != 0 && let Some((part, speed)) = &o.crashed {
                msgs.push(format!("Trümmer zerschellt mit {speed:.0} m/s ({part})"));
                return false;
            }
            true
        });
        for m in msgs {
            self.log(m);
        }
    }

    fn advance_object(&mut self, i: usize, target: f64) {
        let is_ship = i == 0;
        let mut steps = 0;
        while self.objects[i].state.time < target - 1e-9 && steps < MAX_STEPS_PER_FRAME {
            steps += 1;
            let body = &self.body;
            let o = &self.objects[i];
            let remaining = target - o.state.time;

            // Wünsche von Autopilot oder Spieler
            let mut requested = remaining;
            let mut plan_target = None;
            let mut throttle = 0.0;
            let mut stage = false;
            if is_ship {
                if let Some(ap) = self.autopilot.as_mut() {
                    let (plan, event) = ap.update(&o.vessel, &o.state, body, o.last.as_ref());
                    requested = requested.min(plan.requested_dt);
                    plan_target = Some(plan.target);
                    throttle = plan.throttle;
                    stage = plan.stage;
                    if let Some(e) = event {
                        let msg = match e {
                            Event::AscentCutoff { .. } => "Autopilot: Brennschluss, Freiflug zur Apoapsis".to_string(),
                            Event::CircularizationStart { delta_v, .. } => format!("Autopilot: Zirkularisierung, Δv {delta_v:.0} m/s"),
                            Event::Done { .. } => "Autopilot: Umlaufbahn erreicht".to_string(),
                        };
                        let done = ap.phase() == Phase::Done;
                        self.log(msg);
                        if done {
                            self.autopilot = None;
                            self.input.throttle = 0.0;
                        }
                        continue;
                    }
                } else {
                    throttle = self.input.throttle;
                }
            }
            if stage {
                self.stage();
                continue;
            }

            let o = &self.objects[i];
            let dt = suggest_dt(requested, &o.vessel, &o.state, o.last.as_ref(), body);
            let controls = if is_ship { self.ship_controls(dt, throttle, plan_target) } else { Controls::default() };
            if is_ship {
                self.applied_throttle = controls.throttle;
            }

            let o = &mut self.objects[i];
            let report = step(&mut o.vessel, &mut o.state, &controls, &self.body, dt);
            if let Some((part, speed)) = report.crash {
                o.crashed = Some((o.vessel.parts[part].def.name.clone(), speed));
                o.last = Some(report);
                break;
            }
            let broken = report.loads.broken().next().map(|j| (j.child, o.vessel.joint_label(j.child)));
            o.last = Some(report);

            if let Some((child, label)) = broken {
                let (rest, piece) = separate(&o.vessel, &o.state, child, 0.0);
                o.vessel = rest.vessel;
                o.state = rest.state;
                o.version += 1;
                let id = self.next_id;
                self.next_id += 1;
                self.objects.push(SimObject { id, vessel: piece.vessel, state: piece.state, last: None, version: 0, crashed: None });
                self.log(format!("Bruch: {label}"));
            }
        }
        if is_ship && self.objects[0].state.time < target - 1e-6 {
            self.warp_limited = true;
        }
    }

    fn ship_controls(&mut self, dt: f64, throttle: f64, autopilot_target: Option<Vec3<Inertial>>) -> Controls {
        let o = &self.objects[0];
        let thrust = o.last.as_ref().map_or(0.0, |r| r.thrust);
        // Bei langen Schritten (Zeitraffer) wäre jede Regelung instabil.
        if dt > 0.1 {
            return Controls { throttle, rotation: Vec3::ZERO };
        }
        if let Some(target) = autopilot_target {
            return Controls { throttle, rotation: attitude_control(&o.vessel, &o.state, target, thrust.max(1.0)) };
        }
        // Handsteuerung gibt eine Drehrate vor, der Regler setzt sie mit
        // allen Stellmitteln um. So bleibt das Schiff bei jeder Masse gleich
        // gut steuerbar und dreht nicht unkontrolliert durch.
        let manual = self.input.rotation;
        if manual != Vec3::ZERO {
            self.sas_target = None;
            let want = Vec3::new(manual.x() * MANUAL_RATE, manual.y() * MANUAL_RATE, manual.z() * MANUAL_ROLL_RATE);
            return Controls { throttle, rotation: rate_control(&o.vessel, &o.state, want, thrust.max(1.0)) };
        }
        if !self.input.sas {
            self.sas_target = None;
            return Controls { throttle, rotation: Vec3::ZERO };
        }
        // Lagehalter: erst die Drehung abbremsen, dann die Lage halten, in
        // der das Schiff zum Stehen kommt (kein Zurückschwingen).
        if self.sas_target.is_none() {
            let w = o.state.omega;
            if (w.x().powi(2) + w.y().powi(2)).sqrt() > 0.005 {
                return Controls { throttle, rotation: rate_control(&o.vessel, &o.state, Vec3::ZERO, thrust.max(1.0)) };
            }
            self.sas_target = Some(o.state.nose());
        }
        let hold = self.sas_target.unwrap_or(o.state.nose());
        Controls { throttle, rotation: attitude_control(&o.vessel, &o.state, hold, thrust.max(1.0)) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(s: &mut Sim, seconds: f64) {
        let t0 = s.time();
        let n = (seconds * 10.0).round() as usize;
        for k in 1..=n {
            s.advance_to(t0 + k as f64 * 0.1);
        }
    }

    fn nose_rate_deg(s: &Sim) -> f64 {
        let w = s.ship().state.omega;
        (w.x().powi(2) + w.y().powi(2)).sqrt().to_degrees()
    }

    #[test]
    fn throttle_shows_on_the_pad() {
        let mut s = Sim::new();
        s.input.throttle = 0.5;
        run(&mut s, 0.5);
        assert_eq!(s.applied_throttle, 0.5);
    }

    #[test]
    fn staging_without_throttle_leaves_rocket_standing() {
        let mut s = Sim::new();
        s.stage();
        run(&mut s, 5.0);
        assert!(!s.ship_lost && !s.clamped);
        let r = s.ship().last.as_ref().unwrap();
        assert!(r.touching_ground, "steht nicht auf dem Boden");
        assert!(s.ship().state.surface_velocity(&s.body).length() < 0.01);
    }

    /// Start von Hand, 10 s Steigflug, dann 2 s Gieren und loslassen.
    /// Rückgabe: Höhe vor dem Gieren, Drehrate beim Gieren, Restrate danach.
    fn launch_and_yaw(vacuum: bool) -> (f64, f64, f64) {
        let mut s = Sim::new();
        if vacuum {
            s.body.atmosphere = sim_core::Atmosphere::VACUUM;
        }
        s.input.throttle = 1.0;
        s.stage();
        run(&mut s, 10.0);
        let h = s.ship().state.altitude(&s.body);
        s.input.rotation = Vec3::new(0.0, 1.0, 0.0);
        run(&mut s, 2.0);
        let rate = nose_rate_deg(&s);
        s.input.rotation = Vec3::ZERO;
        run(&mut s, 3.0);
        (h, rate, nose_rate_deg(&s))
    }

    #[test]
    fn manual_launch_and_steering() {
        // Ohne Luft: die gewünschte Drehrate (15°/s) wird erreicht
        let (h, rate, after) = launch_and_yaw(true);
        assert!(h > 200.0, "Höhe nach 10 s: {h} m");
        assert!((rate - 15.0).abs() < 1.5, "Drehrate {rate:.1}°/s");
        assert!(after < 0.5, "Lagehalter stoppt nicht: {after:.1}°/s");

        // In dichter Luft stabilisieren die Finnen: die Rakete dreht langsamer
        let (_, rate_air, after_air) = launch_and_yaw(false);
        assert!(rate_air > 5.0 && rate_air < rate, "Drehrate in Luft {rate_air:.1}°/s");
        assert!(after_air < 1.0, "Lagehalter in Luft: {after_air:.1}°/s");
    }

    /// Vollgas senkrecht, Lagehalter aus: größter Anstellwinkel ab 50 m/s
    /// in den ersten 30 s (Grad).
    fn hands_off_launch(fins: bool) -> f64 {
        let mut s = Sim::with_fins(fins);
        s.input.throttle = 1.0;
        s.input.sas = false;
        s.stage();
        let mut worst: f64 = 0.0;
        for _ in 0..30 {
            run(&mut s, 1.0);
            let st = &s.ship().state;
            let surf = st.surface_velocity(&s.body);
            if surf.length() > 50.0 {
                worst = worst.max(st.nose().dot(surf.normalize_or_zero()).clamp(-1.0, 1.0).acos().to_degrees());
            }
        }
        worst
    }

    #[test]
    fn without_fins_and_sas_the_rocket_flips() {
        let with = hands_off_launch(true);
        let without = hands_off_launch(false);
        assert!(with < 10.0, "mit Finnen: {with:.1}°");
        assert!(without > 90.0, "ohne Finnen: {without:.1}°");
    }

    #[test]
    fn switching_design_restarts() {
        let a = Sim::with_fins(true);
        let b = Sim::with_fins(false);
        assert_ne!(a.generation, b.generation);
        assert!(a.ship().vessel.len() > b.ship().vessel.len());
    }

    #[test]
    fn manual_input_takes_over_from_autopilot() {
        let mut s = Sim::new();
        s.toggle_autopilot();
        run(&mut s, 3.0);
        assert!(s.autopilot.is_some());
        s.take_manual_control();
        assert!(s.autopilot.is_none());
        assert_eq!(s.input.throttle, 1.0, "Gas bleibt, wo der Autopilot es hatte");
    }

    #[test]
    fn falling_back_crashes_and_freezes() {
        let mut s = Sim::new();
        s.input.throttle = 1.0;
        s.stage();
        run(&mut s, 6.0);
        s.input.throttle = 0.0;
        run(&mut s, 60.0);
        assert!(s.ship_lost, "Absturz aus ~ 1 km muss zerstören");
        assert!(s.ship().crashed.is_some());
    }
}
