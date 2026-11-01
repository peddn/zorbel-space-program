//! Spielzustand: alle Objekte im Flug, Steuerung, Zeitraffer.
//!
//! Reine Logik ohne Darstellung. Jedes Objekt (Schiff oder Trümmer) wird
//! mit `vessel::dynamics::step` fortgeschrieben, genau wie im Testflug.

use bevy::prelude::Resource;
use sim_core::kepler::{self, Elements};
use sim_core::{CelestialBody, Inertial, Vec3, VesselFrame};
use vessel::autopilot::{AscentAutopilot, Event, Phase, attitude_control};
use vessel::dynamics::{activate_stage, separate, step, suggest_dt};
use vessel::scenario::{first_rocket, on_launch_site, test_planet};
use vessel::{Controls, PartLibrary, RigidState, StepReport, Vessel};

/// Zeitraffer-Stufen. Ab 10× regelt niemand mehr die Lage (Schritte zu lang).
pub const WARP_LEVELS: [f64; 8] = [1.0, 2.0, 4.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0];
/// Höchstens so viele Physikschritte pro Objekt und Bild.
const MAX_STEPS_PER_FRAME: usize = 4_000;

pub struct SimObject {
    /// Eindeutig über die ganze Sitzung; die Darstellung hängt daran.
    pub id: u64,
    pub vessel: Vessel,
    pub state: RigidState,
    pub last: Option<StepReport>,
    /// Zählt Änderungen an der Teileliste (Neuaufbau der Meshes).
    pub version: u64,
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
    /// Zuletzt tatsächlich gegebenes Gas (Spieler oder Autopilot).
    pub applied_throttle: f64,
}

impl Sim {
    pub fn new() -> Self {
        let body = test_planet();
        let lib = PartLibrary::builtin();
        let vessel = first_rocket(&lib);
        let state = on_launch_site(&body, &vessel, 0.0);
        Self {
            body,
            objects: vec![SimObject { id: 0, vessel, state, last: None, version: 0 }],
            input: PilotInput { sas: true, ..Default::default() },
            autopilot: None,
            clamped: true,
            warp_index: 0,
            warp_limited: false,
            sas_target: None,
            next_id: 1,
            messages: vec![(0.0, "Bereit. Leertaste: Start".into())],
            ship_lost: false,
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
            self.objects.push(SimObject { id, vessel: d.vessel, state: d.state, last: None, version: 0 });
        }
        self.log(if n > 0 { "Stufentrennung" } else { "Zündung" });
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
            // Auf der Rampe: steht fest und dreht mit dem Planeten.
            let ship = &mut self.objects[0];
            ship.state = on_launch_site(&self.body, &ship.vessel, target);
            return;
        }
        let mut i = 0;
        while i < self.objects.len() {
            self.advance_object(i, target);
            i += 1;
        }
        // Abgestürzte Trümmer entfernen
        let body = &self.body;
        let before = self.objects.len();
        let ship_crashed = self.objects[0].state.altitude(body) <= 0.0;
        self.objects.retain(|o| o.id == 0 || o.state.altitude(body) > 0.0);
        let lost = before - self.objects.len();
        if lost > 0 {
            self.log(format!("{lost} Trümmerteil(e) aufgeschlagen"));
        }
        if ship_crashed {
            self.ship_lost = true;
            self.log("Aufschlag! (R: Neustart)");
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
            let dt = suggest_dt(requested, &o.state, o.last.as_ref(), body);
            let controls = if is_ship { self.ship_controls(dt, throttle, plan_target) } else { Controls::default() };
            if is_ship {
                self.applied_throttle = controls.throttle;
            }

            let o = &mut self.objects[i];
            let report = step(&mut o.vessel, &mut o.state, &controls, &self.body, dt);
            let broken = report.loads.broken().next().map(|j| (j.child, o.vessel.joint_label(j.child)));
            o.last = Some(report);

            if let Some((child, label)) = broken {
                let (rest, piece) = separate(&o.vessel, &o.state, child, 0.0);
                o.vessel = rest.vessel;
                o.state = rest.state;
                o.version += 1;
                let id = self.next_id;
                self.next_id += 1;
                self.objects.push(SimObject { id, vessel: piece.vessel, state: piece.state, last: None, version: 0 });
                self.log(format!("Bruch: {label}"));
            }
            if self.objects[i].state.altitude(&self.body) <= 0.0 {
                break;
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
        let manual = self.input.rotation;
        if manual != Vec3::ZERO {
            self.sas_target = None;
            return Controls { throttle, rotation: manual };
        }
        if !self.input.sas {
            self.sas_target = None;
            return Controls { throttle, rotation: Vec3::ZERO };
        }
        let hold = *self.sas_target.get_or_insert(o.state.nose());
        Controls { throttle, rotation: attitude_control(&o.vessel, &o.state, hold, thrust.max(1.0)) }
    }
}
