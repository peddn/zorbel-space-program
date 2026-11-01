//! Screenshot-Modus: ohne Fenster einen festen Flug abspielen und an
//! markanten Punkten Bilder speichern. Dient zur visuellen Prüfung
//! (auch ohne GPU, z. B. mit Mesas Software-Treiber) und als visueller
//! Regressionstest.

use crate::ZorbelPlugin;
use crate::camera::{OrbitCamera, camera_bundle};
use crate::sim::Sim;
use bevy::app::SubApps;
use bevy::asset::RenderAssetUsages;
use bevy::camera::RenderTarget;
use bevy::prelude::*;
use bevy::render::RenderPlugin;
use bevy::render::render_resource::{PollType, TextureFormat};
use bevy::render::renderer::RenderDevice;
use bevy::render::view::screenshot::{Screenshot, save_to_disk};
use bevy::window::ExitCondition;
use bevy::winit::WinitPlugin;

const WIDTH: u32 = 1280;
const HEIGHT: u32 = 720;

struct Shot {
    name: &'static str,
    /// Vorbereitung vor dem Abwarten (z. B. neues Szenario).
    prepare: fn(&mut Sim),
    /// Simulation läuft weiter, bis das erfüllt ist.
    until: fn(&Sim, f64) -> bool,
    camera: OrbitCamera,
}

fn shots() -> Vec<Shot> {
    let cam = |distance: f32, pitch: f32, yaw: f32| OrbitCamera { yaw, pitch, distance };
    vec![
        Shot { prepare: |_| {}, name: "01-startrampe", until: |_, _| true, camera: cam(32.0, 0.12, -0.7) },
        Shot { prepare: |_| {}, name: "02-abheben", until: |s, _| s.time() >= 8.0, camera: cam(45.0, 0.05, -0.5) },
        Shot { prepare: |_| {}, name: "03-max-q", until: |s, _| s.time() >= 57.0, camera: cam(40.0, 0.2, -1.1) },
        Shot { prepare: |_| {}, name: "04-stufentrennung", until: |s, t0| s.objects.len() > 1 && s.time() >= t0 + 2.0, camera: cam(55.0, 0.25, -0.4) },
        Shot {
            prepare: |_| {},
            name: "05-weltraum",
            until: |s, _| s.ship().state.altitude(&s.body) > 75_000.0,
            camera: cam(30.0, 0.35, -0.9),
        },
        Shot { prepare: |_| {}, name: "06-umlaufbahn", until: |s, _| s.autopilot.is_none() && s.time() > 200.0, camera: cam(22.0, 0.15, -0.6) },
        Shot { prepare: |_| {}, name: "07-karte", until: |_, _| true, camera: cam(3.2e6, 0.7, 0.4) },
        // Neuer Start von Hand: Klammern lösen ohne Gas, die Rakete muss stehen bleiben.
        Shot {
            prepare: |s| {
                *s = Sim::new();
                s.stage();
            },
            name: "08-ohne-gas-am-boden",
            until: |s, _| s.time() >= 3.0,
            camera: cam(32.0, 0.12, -0.7),
        },
    ]
}

pub fn run(dir: &str) {
    std::fs::create_dir_all(dir).expect("Verzeichnis anlegen");
    let mut app = App::new();
    app.add_plugins(
        DefaultPlugins
            .set(WindowPlugin { primary_window: None, exit_condition: ExitCondition::DontExit, ..default() })
            .set(RenderPlugin { synchronous_pipeline_compilation: true, ..default() })
            .disable::<WinitPlugin>(),
    )
    .add_plugins(ZorbelPlugin);

    let mut target = Image::new_target_texture(WIDTH, HEIGHT, TextureFormat::Rgba8UnormSrgb, None);
    target.asset_usage = RenderAssetUsages::RENDER_WORLD;
    let handle = app.world_mut().resource_mut::<Assets<Image>>().add(target);
    app.world_mut().spawn((camera_bundle(), RenderTarget::from(handle.clone())));

    app.finish();
    app.cleanup();
    let mut apps: SubApps = std::mem::take(app.sub_apps_mut());

    let update = |apps: &mut SubApps| {
        apps.update();
        apps.main
            .world()
            .resource::<RenderDevice>()
            .wgpu_device()
            .poll(PollType::Wait { submission_index: None, timeout: None })
            .unwrap();
    };

    // Startup-Systeme laufen beim ersten Update.
    update(&mut apps);
    let mut started = false;
    for shot in shots() {
        (shot.prepare)(&mut apps.main.world_mut().resource_mut::<Sim>());
        let t0 = apps.main.world().resource::<Sim>().time();
        loop {
            let world = apps.main.world_mut();
            let mut sim = world.resource_mut::<Sim>();
            if (shot.until)(&sim, t0) || sim.ship_lost {
                break;
            }
            if !started {
                sim.toggle_autopilot();
                started = true;
            }
            let t = sim.time();
            sim.warp_index = 0;
            // Freiflug darf in großen Schritten laufen, sonst halbe Sekunden.
            let chunk = if sim.autopilot.as_ref().is_some_and(|a| a.phase() == vessel::autopilot::Phase::Coast) { 20.0 } else { 0.5 };
            sim.advance_to(t + chunk);
        }
        *apps.main.world_mut().resource_mut::<OrbitCamera>() = shot.camera;
        // Ein paar Bilder, bis Meshes, Bodenstück und Texturen bereitstehen
        for _ in 0..4 {
            update(&mut apps);
        }
        let path = format!("{dir}/{}.png", shot.name);
        apps.main.world_mut().spawn(Screenshot::image(handle.clone())).observe(save_to_disk(path.clone()));
        for _ in 0..3 {
            update(&mut apps);
        }
        let sim = apps.main.world().resource::<Sim>();
        println!("{path}  T+{:.1}s  Höhe {:.1} km", sim.time(), sim.ship().state.altitude(&sim.body) / 1000.0);
    }
}
