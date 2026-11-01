//! Zorbel Space Program – Darstellung und Steuerung.
//!
//! `zorbel-game`               Spiel im Fenster
//! `zorbel-game --screenshots` ohne Fenster einen festen Flug abspielen und
//!                             Bilder nach `screenshots/` schreiben

// Bevy-Abfragen mit Filtern sind lang; das ist dort üblich.
#![allow(clippy::type_complexity)]

mod camera;
mod hud;
mod render;
mod screenshot;
mod sim;

use bevy::prelude::*;

/// Alles, was Spiel und Screenshot-Modus gemeinsam haben.
pub struct ZorbelPlugin;

impl Plugin for ZorbelPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(sim::Sim::new())
            .insert_resource(camera::OrbitCamera::default())
            .insert_resource(render::Focus::default())
            .insert_resource(ClearColor(Color::BLACK))
            .insert_resource(render::ShowCenters::default())
            .init_gizmo_group::<render::MarkerGizmos>()
            .add_systems(Startup, (render::setup_bodies, render::setup_part_materials, render::setup_marker_gizmos))
            .add_systems(PostStartup, hud::setup_hud)
            .add_systems(
                Update,
                (
                    render::update_focus,
                    render::sync_bodies,
                    render::sync_sky_stars,
                    render::sync_ground_patch,
                    render::sync_vessels,
                    render::draw_orbits,
                    render::draw_centers,
                    camera::place_camera,
                    hud::update_hud,
                )
                    .chain(),
            );
    }
}

fn main() {
    if std::env::args().any(|a| a == "--screenshots") {
        screenshot::run("screenshots");
        return;
    }
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window { title: "Zorbel Space Program".into(), ..default() }),
            ..default()
        }))
        .add_plugins(ZorbelPlugin)
        .add_systems(Startup, |mut commands: Commands| {
            commands.spawn(camera::camera_bundle());
        })
        .add_systems(
            Update,
            (hud::keyboard, camera::camera_input, hud::advance_sim).chain().before(render::update_focus),
        )
        .run();
}
