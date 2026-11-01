//! Aufstiegs-Autopilot: vom Startplatz in eine Kreisbahn.
//!
//! Wird vom Headless-Testflug, vom Screenshot-Modus und im Spiel (Taste)
//! benutzt. Der Autopilot entscheidet nur Schub, Zielrichtung, gewünschte
//! Schrittweite und Stufen; geflogen wird mit derselben Physik wie von Hand.

use crate::dynamics::{Controls, DT_PHYSICS, RigidState, StepReport};
use crate::mass::mass_properties;
use crate::vessel::Vessel;
use sim_core::kepler::{self, Elements};
use sim_core::{CelestialBody, Inertial, Vec3, VesselFrame};

/// Lageregler: dreht die Spitze in Richtung `target`, dämpft Rollen.
/// Die Verstärkung folgt aus Trägheit und verfügbarem Stellmoment
/// (Reaktionsrad + Schubvektor), damit er bei jeder Masse gleich schnell
/// und stabil reagiert.
pub fn attitude_control(v: &Vessel, s: &RigidState, target: Vec3<Inertial>, thrust: f64) -> Vec3<VesselFrame> {
    let nose = s.nose();
    let axis = nose.cross(target);
    let angle = axis.length().atan2(nose.dot(target));
    let axis_body = s.attitude.inverse().apply(axis.normalize_or_zero());
    let want = axis_body * (angle * 0.8).min(0.2);
    rate_control(v, s, Vec3::new(want.x(), want.y(), 0.0), thrust)
}

/// Drehratenregler: Steuerkommando (−1..1 je Achse), das die
/// Winkelgeschwindigkeit auf `want` (rad/s, Schiffssystem) bringt. Nutzt
/// Reaktionsrad, Schubvektor und Steuerflächen gemeinsam.
pub fn rate_control(v: &Vessel, s: &RigidState, want: Vec3<VesselFrame>, thrust: f64) -> Vec3<VesselFrame> {
    const RESPONSE: f64 = 0.4; // s, Zeitkonstante der Drehratenregelung
    let err = want - s.omega;

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

/// Dreht `from` höchstens um `max_angle` in Richtung `to` (Einheitsvektoren).
pub fn limit_towards(from: Vec3<Inertial>, to: Vec3<Inertial>, max_angle: f64) -> Vec3<Inertial> {
    let angle = from.dot(to).clamp(-1.0, 1.0).acos();
    if angle <= max_angle {
        return to;
    }
    let perp = to.reject_from_unit(from).normalize_or_zero();
    from * max_angle.cos() + perp * max_angle.sin()
}

/// Summe des Vakuumschubs aller gezündeten Triebwerke.
pub fn active_thrust_vac(v: &Vessel) -> f64 {
    v.parts.iter().filter(|p| p.engine_active).filter_map(|p| p.def.engine()).map(|e| e.thrust_vac).sum()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Ascent,
    Coast,
    Circularize,
    Done,
}

/// Was der Autopilot im aktuellen Augenblick will.
#[derive(Debug, Clone, Copy)]
pub struct Plan {
    pub throttle: f64,
    /// Gewünschte Ausrichtung der Spitze.
    pub target: Vec3<Inertial>,
    /// Gewünschte Schrittweite (wird von `suggest_dt` weiter begrenzt).
    pub requested_dt: f64,
    /// Nächste Stufe auslösen.
    pub stage: bool,
}

impl Plan {
    /// Steuerbefehle für einen Schritt der Länge `dt`. Bei großen
    /// Schritten (Zeitraffer) regelt der Autopilot die Lage nicht, ein
    /// expliziter Regler wäre dort instabil.
    pub fn controls(&self, v: &Vessel, s: &RigidState, dt: f64, current_thrust: f64) -> Controls {
        let rotation = if dt <= 0.1 {
            attitude_control(v, s, self.target, current_thrust.max(self.throttle * active_thrust_vac(v) * 0.5))
        } else {
            Vec3::ZERO
        };
        Controls { throttle: self.throttle, rotation }
    }
}

/// Ereignisse für Anzeige und Protokoll.
#[derive(Debug, Clone, Copy)]
pub enum Event {
    /// Apoapsis erreicht, Triebwerk aus.
    AscentCutoff { elements: Elements },
    CircularizationStart { delta_v: f64, burn_time: f64 },
    Done { elements: Elements },
}

#[derive(Debug, Clone)]
pub struct AscentAutopilot {
    pub target_altitude: f64,
    /// Größter Anstellwinkel beim Aufstieg in rad (begrenzt Querlasten).
    pub max_aoa: f64,
    /// Höhe, bei der die Flugbahn waagrecht sein soll (Neigungsprogramm).
    pub turn_end_altitude: f64,
    phase: Phase,
}

impl AscentAutopilot {
    pub fn new(target_altitude: f64) -> Self {
        Self { target_altitude, max_aoa: 4f64.to_radians(), turn_end_altitude: 48_000.0, phase: Phase::Ascent }
    }

    pub fn phase(&self) -> Phase {
        self.phase
    }

    /// Plan für den nächsten Schritt; schaltet bei Bedarf die Phase weiter.
    pub fn update(
        &mut self,
        v: &Vessel,
        s: &RigidState,
        body: &CelestialBody,
        last: Option<&StepReport>,
    ) -> (Plan, Option<Event>) {
        let up = s.pos.normalize_or_zero();
        let east = Vec3::<Inertial>::Z.cross(up).normalize_or_zero();
        let speed = s.vel.length();
        let prograde = if speed > 1.0 { s.vel / speed } else { up };
        let surf = s.surface_velocity(body);
        let surf_speed = surf.length();
        let surf_prograde = if surf_speed > 1.0 { surf / surf_speed } else { up };
        let el = kepler::elements(s.orbit(), body.mu);
        let apo_alt = el.apoapsis.map_or(f64::INFINITY, |a| a - body.radius);
        let alt = s.altitude(body);
        let horizontal = s.vel.reject_from_unit(up);
        let stage = last.is_some_and(|r| r.flameout) && v.next_stage().is_some();
        let idle = Plan { throttle: 0.0, target: prograde, requested_dt: f64::INFINITY, stage: false };

        match self.phase {
            Phase::Ascent => {
                if apo_alt >= self.target_altitude {
                    self.phase = Phase::Coast;
                    return (Plan { requested_dt: DT_PHYSICS, ..idle }, Some(Event::AscentCutoff { elements: el }));
                }
                let pitch: f64 = if surf_speed < 60.0 {
                    90.0
                } else {
                    90.0 * (1.0 - (alt / self.turn_end_altitude).clamp(0.0, 1.0).sqrt())
                };
                let program = up * pitch.to_radians().sin() + east * pitch.to_radians().cos();
                let target = if surf_speed < 60.0 { up } else { limit_towards(surf_prograde, program, self.max_aoa) };
                (Plan { throttle: 1.0, target, requested_dt: DT_PHYSICS, stage }, None)
            }
            Phase::Coast => {
                let t_apo = el.time_to_apoapsis.unwrap_or(0.0);
                let r_apo = el.apoapsis.unwrap_or(s.pos.length());
                let v_apo = (body.mu * (2.0 / r_apo - 1.0 / el.semi_major_axis)).sqrt();
                let dv = (body.mu / r_apo).sqrt() - v_apo;
                let thrust = active_thrust_vac(v).max(1.0);
                let burn = dv * v.total_mass() / thrust;
                let lead = burn / 2.0;
                if t_apo <= lead {
                    self.phase = Phase::Circularize;
                    return (
                        Plan { requested_dt: DT_PHYSICS, ..idle },
                        Some(Event::CircularizationStart { delta_v: dv, burn_time: burn }),
                    );
                }
                // Luftwiderstand senkt die Apoapsis noch: nachschieben.
                if alt < body.atmosphere.top() && apo_alt < self.target_altitude - 200.0 {
                    return (Plan { throttle: 0.3, target: prograde, requested_dt: DT_PHYSICS, stage }, None);
                }
                // Bis kurz vor der Zündung frei fliegen, dann ausrichten.
                let requested_dt = if t_apo - lead > 60.0 { t_apo - lead - 50.0 } else { DT_PHYSICS };
                (Plan { target: horizontal.normalize_or_zero(), requested_dt, ..idle }, None)
            }
            Phase::Circularize => {
                let dv_left = (body.mu / s.pos.length()).sqrt() - horizontal.length();
                let empty = v.parts.iter().all(|p| p.def.fuel_capacity() == 0.0 || p.fuel <= 0.0);
                if dv_left <= 0.05 || empty {
                    self.phase = Phase::Done;
                    return (idle, Some(Event::Done { elements: el }));
                }
                let thrust = active_thrust_vac(v).max(1.0);
                let throttle = (dv_left * v.total_mass() / thrust / 2.0).clamp(0.02, 1.0);
                (Plan { throttle, target: horizontal.normalize_or_zero(), requested_dt: DT_PHYSICS, stage }, None)
            }
            Phase::Done => (idle, None),
        }
    }
}
