//! Verfolgerkamera um das aktive Schiff. "Oben" ist die lokale Senkrechte,
//! damit der Horizont waagrecht bleibt.

use crate::render::render_vec;
use crate::sim::Sim;
use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll, MouseScrollUnit};
use bevy::prelude::*;
use sim_core::{Inertial, Vec3 as SimVec3};

#[derive(Resource, Clone, Copy)]
pub struct OrbitCamera {
    /// Winkel um die Senkrechte, 0 = Blick nach Norden.
    pub yaw: f32,
    /// Höhenwinkel über dem Horizont.
    pub pitch: f32,
    pub distance: f32,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self { yaw: -0.7, pitch: 0.15, distance: 32.0 } // Sonne im Rücken
    }
}

#[derive(Component)]
pub struct MainCamera;

pub fn camera_bundle() -> impl Bundle {
    (
        MainCamera,
        Camera3d::default(),
        Projection::Perspective(PerspectiveProjection { fov: 0.9, near: 0.1, far: 1.0e9, ..default() }),
        Transform::default(),
    )
}

pub fn camera_input(
    mut cam: ResMut<OrbitCamera>,
    buttons: Res<ButtonInput<MouseButton>>,
    motion: Res<AccumulatedMouseMotion>,
    scroll: Res<AccumulatedMouseScroll>,
) {
    if buttons.pressed(MouseButton::Right) {
        cam.yaw -= motion.delta.x * 0.005;
        cam.pitch = (cam.pitch + motion.delta.y * 0.005).clamp(-1.5, 1.5);
    }
    // Mausräder melden Rastungen (Zeilen), Touchpads und manche Systeme
    // Pixel; beides auf "Rastungen" umrechnen.
    let notches = match scroll.unit {
        MouseScrollUnit::Line => scroll.delta.y,
        MouseScrollUnit::Pixel => scroll.delta.y / 50.0,
    };
    if notches != 0.0 {
        cam.distance = (cam.distance * 0.88f32.powf(notches)).clamp(4.0, 5.0e7);
    }
}

pub fn place_camera(sim: Res<Sim>, cam: Res<OrbitCamera>, mut q: Query<&mut Transform, With<MainCamera>>) {
    let pos = sim.ship().state.pos;
    let up = render_vec(pos.normalize_or_zero());
    let east = render_vec(SimVec3::<Inertial>::Z.cross(pos).normalize_or_zero());
    let east = if east.length_squared() < 0.5 { Vec3::X } else { east };
    let north = up.cross(east);
    let horizontal = north * cam.yaw.cos() + east * cam.yaw.sin();
    let offset = (horizontal * cam.pitch.cos() + up * cam.pitch.sin()) * cam.distance;
    for mut tr in &mut q {
        *tr = Transform::from_translation(offset).looking_at(Vec3::ZERO, up);
    }
}
