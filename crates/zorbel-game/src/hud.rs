//! Fluganzeige und Tastatursteuerung.

use crate::sim::{Sim, WARP_LEVELS};
use bevy::input::ButtonState;
use bevy::input::keyboard::{Key, KeyboardInput};
use bevy::prelude::*;
use sim_core::{Vec3 as SimVec3, VesselFrame};
use vessel::autopilot::Phase;

#[derive(Component)]
pub struct FlightText;
#[derive(Component)]
pub struct MessageText;

/// DejaVu Sans Mono (freie Lizenz, siehe assets/fonts): enthält Umlaute
/// und Zeichen wie × · Δ, die Bevys eingebaute Schrift nicht hat.
const FONT: &[u8] = include_bytes!("../../../assets/fonts/DejaVuSansMono.ttf");

pub fn setup_hud(
    mut commands: Commands,
    camera: Single<Entity, With<crate::camera::MainCamera>>,
    mut fonts: ResMut<Assets<Font>>,
) {
    let cam = *camera;
    let handle = fonts.add(Font::from_bytes(FONT.to_vec()));
    let font = TextFont { font: handle.into(), font_size: FontSize::Px(14.0), ..default() };
    commands.spawn((
        FlightText,
        Text::new(""),
        font.clone(),
        TextColor(Color::srgb(0.85, 0.95, 1.0)),
        Node { position_type: PositionType::Absolute, top: px(10), left: px(12), ..default() },
        UiTargetCamera(cam),
    ));
    commands.spawn((
        MessageText,
        Text::new(""),
        font,
        TextColor(Color::srgb(1.0, 0.9, 0.6)),
        Node { position_type: PositionType::Absolute, bottom: px(10), left: px(12), ..default() },
        UiTargetCamera(cam),
    ));
}

fn km(m: f64) -> String {
    if m.abs() >= 10_000.0 { format!("{:.1} km", m / 1000.0) } else { format!("{m:.0} m") }
}

pub fn update_hud(sim: Res<Sim>, mut flight: Query<&mut Text, With<FlightText>>, mut msgs: Query<&mut Text, (With<MessageText>, Without<FlightText>)>) {
    let ship = sim.ship();
    let body = &sim.body;
    let s = &ship.state;
    let el = ship.elements(body);
    let r = ship.last.as_ref();
    let fuel: f64 = ship.vessel.parts.iter().map(|p| p.fuel).sum();
    let ap = el.apoapsis.map_or("—".into(), |a| km(a - body.radius));
    let pe = km(el.periapsis - body.radius);
    let autopilot = match sim.autopilot.as_ref().map(|a| a.phase()) {
        Some(Phase::Ascent) => "Aufstieg",
        Some(Phase::Coast) => "Freiflug",
        Some(Phase::Circularize) => "Zirkularisierung",
        Some(Phase::Done) | None => "aus",
    };
    let status = if sim.ship_lost {
        "ZERSTÖRT".to_string()
    } else if sim.clamped {
        "Startklammern halten".to_string()
    } else if r.is_some_and(|r| r.touching_ground) {
        format!("Am Boden, Stützkraft {:.0} kN", r.map_or(0.0, |r| r.ground_force / 1000.0))
    } else {
        "Im Flug".to_string()
    };
    let warp = if sim.warp_limited { format!("{}× (begrenzt)", sim.warp()) } else { format!("{}×", sim.warp()) };
    let text = format!(
        "ZORBEL SPACE PROGRAM\n\
         T+ {:.1} s   Zeitraffer {warp}\n\
         Höhe        {}\n\
         v Oberfl.   {:.0} m/s\n\
         v Bahn      {:.0} m/s\n\
         Apoapsis    {ap}\n\
         Periapsis   {pe}\n\
         Mach {:.2}   q {:.1} kPa\n\
         Schub {:.0} kN   Gas {:.0} %\n\
         Treibstoff  {:.0} kg   Teile {}\n\
         Lagehalter  {}   Autopilot {autopilot}\n\
         {status}",
        s.time,
        km(s.altitude(body)),
        s.surface_velocity(body).length(),
        s.vel.length(),
        r.map_or(0.0, |r| r.mach),
        r.map_or(0.0, |r| r.dynamic_pressure / 1000.0),
        r.map_or(0.0, |r| r.thrust / 1000.0),
        sim.applied_throttle * 100.0,
        fuel,
        ship.vessel.len(),
        if sim.input.sas { "an" } else { "aus" },
    );
    for mut t in &mut flight {
        t.0 = text.clone();
    }
    let mut m: String = sim.messages.iter().map(|(t, s)| format!("T+{t:>7.1}  {s}\n")).collect();
    m.push_str("\nLeer Stufe · Shift/Strg Gas · Z/X voll/aus · WASD QE Lage (Eingriff schaltet Autopilot ab) · T Lagehalter · P Autopilot\n, . Zeitraffer · M Karte · rechte Maus Kamera · Rad Zoom · R Neustart");
    for mut t in &mut msgs {
        t.0 = m.clone();
    }
}

/// Befehl einer Taste, die einmal pro Druck wirkt.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Command {
    Stage,
    Autopilot,
    Sas,
    WarpUp,
    WarpDown,
    FullThrottle,
    CutThrottle,
    Map,
    Restart,
}

/// Befehlstasten nach ihrer Beschriftung (logische Taste), damit z. B. „Z“
/// auf QWERTZ und QWERTY dieselbe Taste ist. Gehaltene Tasten (Lage, Gas)
/// werden dagegen nach ihrer Position abgefragt, wie in Spielen üblich.
pub fn command_for(ev: &KeyboardInput) -> Option<Command> {
    if ev.state != ButtonState::Pressed || ev.repeat {
        return None;
    }
    if ev.key_code == KeyCode::Space {
        return Some(Command::Stage);
    }
    let Key::Character(c) = &ev.logical_key else { return None };
    Some(match c.to_lowercase().as_str() {
        "p" => Command::Autopilot,
        "t" => Command::Sas,
        "." | ":" => Command::WarpUp,
        "," | ";" => Command::WarpDown,
        "z" => Command::FullThrottle,
        "x" => Command::CutThrottle,
        "m" => Command::Map,
        "r" => Command::Restart,
        _ => return None,
    })
}

pub fn keyboard(
    mut events: MessageReader<KeyboardInput>,
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    mut sim: ResMut<Sim>,
    mut cam: ResMut<crate::camera::OrbitCamera>,
) {
    let dt = time.delta_secs_f64();
    let commands: Vec<Command> = events.read().filter_map(command_for).collect();
    for cmd in commands {
        match cmd {
            Command::Restart => {
                *sim = Sim::new();
                return;
            }
            Command::Stage => sim.stage(),
            Command::Autopilot => sim.toggle_autopilot(),
            Command::Sas => sim.input.sas = !sim.input.sas,
            Command::WarpUp => sim.warp_index = (sim.warp_index + 1).min(WARP_LEVELS.len() - 1),
            Command::WarpDown => sim.warp_index = sim.warp_index.saturating_sub(1),
            Command::FullThrottle => {
                sim.take_manual_control();
                sim.input.throttle = 1.0;
            }
            Command::CutThrottle => {
                sim.take_manual_control();
                sim.input.throttle = 0.0;
            }
            Command::Map => cam.distance = if cam.distance < 1.0e5 { 3.5e6 } else { 32.0 },
        }
    }

    let up = keys.pressed(KeyCode::ShiftLeft) || keys.pressed(KeyCode::ShiftRight);
    let down = keys.pressed(KeyCode::ControlLeft) || keys.pressed(KeyCode::ControlRight);
    if up || down {
        sim.take_manual_control();
        let delta = (up as i32 - down as i32) as f64 * dt;
        sim.input.throttle = (sim.input.throttle + delta).clamp(0.0, 1.0);
    }

    let axis = |pos: KeyCode, neg: KeyCode| (keys.pressed(pos) as i32 - keys.pressed(neg) as i32) as f64;
    let rotation =
        SimVec3::<VesselFrame>::new(axis(KeyCode::KeyS, KeyCode::KeyW), axis(KeyCode::KeyD, KeyCode::KeyA), axis(KeyCode::KeyE, KeyCode::KeyQ));
    if rotation != SimVec3::ZERO {
        sim.take_manual_control();
    }
    sim.input.rotation = rotation;
}

pub fn advance_sim(time: Res<Time>, mut sim: ResMut<Sim>) {
    sim.advance(time.delta_secs_f64().min(0.1));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn press(code: KeyCode, text: &str) -> KeyboardInput {
        KeyboardInput {
            key_code: code,
            logical_key: Key::Character(text.into()),
            state: ButtonState::Pressed,
            text: None,
            repeat: false,
            window: Entity::PLACEHOLDER,
        }
    }

    #[test]
    fn commands_follow_the_key_label() {
        // QWERTZ: die Taste mit der Aufschrift Z sitzt auf der Position "Y"
        assert_eq!(command_for(&press(KeyCode::KeyY, "z")), Some(Command::FullThrottle));
        // mit gedrückter Shift-Taste (Gas auf) kommt ein großes Z
        assert_eq!(command_for(&press(KeyCode::KeyZ, "Z")), Some(Command::FullThrottle));
        assert_eq!(command_for(&press(KeyCode::KeyZ, "y")), None);
        assert_eq!(command_for(&press(KeyCode::Period, ":")), Some(Command::WarpUp));
        let mut repeat = press(KeyCode::KeyP, "p");
        repeat.repeat = true;
        assert_eq!(command_for(&repeat), None);
    }
}
