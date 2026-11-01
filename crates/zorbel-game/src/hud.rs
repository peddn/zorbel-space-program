//! Fluganzeige und Tastatursteuerung.

use crate::sim::{Sim, WARP_LEVELS};
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
         Lagehalter  {}   Autopilot {autopilot}",
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
    m.push_str("\nLeer Stufe · Shift/Strg Gas · Z/X voll/aus · WASD QE Lage · T Lagehalter · P Autopilot\n, . Zeitraffer · M Karte · rechte Maus Kamera · Rad Zoom · R Neustart");
    for mut t in &mut msgs {
        t.0 = m.clone();
    }
}

pub fn keyboard(keys: Res<ButtonInput<KeyCode>>, time: Res<Time>, mut sim: ResMut<Sim>) {
    let dt = time.delta_secs_f64();
    if keys.just_pressed(KeyCode::KeyR) {
        *sim = Sim::new();
        return;
    }
    if keys.just_pressed(KeyCode::Space) {
        sim.stage();
    }
    if keys.just_pressed(KeyCode::KeyP) {
        sim.toggle_autopilot();
    }
    if keys.just_pressed(KeyCode::KeyT) {
        sim.input.sas = !sim.input.sas;
    }
    if keys.just_pressed(KeyCode::Period) {
        sim.warp_index = (sim.warp_index + 1).min(WARP_LEVELS.len() - 1);
    }
    if keys.just_pressed(KeyCode::Comma) {
        sim.warp_index = sim.warp_index.saturating_sub(1);
    }
    let mut throttle = sim.input.throttle;
    if keys.pressed(KeyCode::ShiftLeft) {
        throttle += dt;
    }
    if keys.pressed(KeyCode::ControlLeft) {
        throttle -= dt;
    }
    if keys.just_pressed(KeyCode::KeyZ) {
        throttle = 1.0;
    }
    if keys.just_pressed(KeyCode::KeyX) {
        throttle = 0.0;
    }
    sim.input.throttle = throttle.clamp(0.0, 1.0);
    let axis = |pos: KeyCode, neg: KeyCode| (keys.pressed(pos) as i32 - keys.pressed(neg) as i32) as f64;
    sim.input.rotation =
        SimVec3::<VesselFrame>::new(axis(KeyCode::KeyS, KeyCode::KeyW), axis(KeyCode::KeyD, KeyCode::KeyA), axis(KeyCode::KeyE, KeyCode::KeyQ));
}

pub fn advance_sim(time: Res<Time>, mut sim: ResMut<Sim>) {
    sim.advance(time.delta_secs_f64().min(0.1));
}
