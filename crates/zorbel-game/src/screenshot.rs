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
    /// Schwer- und Druckpunkt einblenden.
    centers: bool,
}

fn shot(name: &'static str, until: fn(&Sim, f64) -> bool, distance: f32, pitch: f32, yaw: f32) -> Shot {
    Shot { name, prepare: |_| {}, until, camera: OrbitCamera { yaw, pitch, distance }, centers: false }
}

fn shots() -> Vec<Shot> {
    vec![
        shot("01-startrampe", |_, _| true, 32.0, 0.12, -0.7),
        shot("02-abheben", |s, _| s.time() >= 8.0, 45.0, 0.05, -0.5),
        Shot { centers: true, ..shot("03-max-q", |s, _| s.time() >= 57.0, 40.0, 0.2, -1.1) },
        shot("04-stufentrennung", |s, t0| s.objects.len() > 1 && s.time() >= t0 + 2.0, 55.0, 0.25, -0.4),
        shot("05-weltraum", |s, _| s.ship().state.altitude(s.body()) > 75_000.0, 30.0, 0.35, -0.9),
        shot("06-umlaufbahn", |s, _| s.autopilot.is_none() && s.time() > 200.0, 22.0, 0.15, -0.6),
        shot("07-karte", |_, _| true, 3.2e6, 0.7, 0.4),
        // Das ganze System Zol: Zun, Bahnen von Zorb und Muz
        shot("07b-system-zol", |_, _| true, 4.0e10, 1.3, 0.4),
        // Niedrige Bahn um Muz, Zorb im Hintergrund
        Shot { prepare: orbit_muz, ..shot("07c-muz", |_, _| true, 9.0e5, 0.25, 1.4) },
        // Neuer Start von Hand: Klammern lösen ohne Gas, die Rakete muss stehen bleiben.
        Shot {
            prepare: |s| {
                *s = Sim::new();
                s.stage();
            },
            ..shot("08-ohne-gas-am-boden", |s, _| s.time() >= 3.0, 32.0, 0.12, -0.7)
        },
        // Rakete ohne Finnen, Vollgas, Lagehalter aus: Druckpunkt vor dem Schwerpunkt.
        Shot {
            prepare: |s| {
                *s = Sim::with_fins(false);
                s.input.throttle = 1.0;
                s.input.sas = false;
                s.stage();
            },
            centers: true,
            ..shot("09-ohne-finnen-steigflug", |s, _| s.time() >= 8.0, 40.0, 0.15, -0.9)
        },
        // Ohne Stoß von außen: schon die Seitenströmung durch die Planetendrehung
        // reicht, die instabile Rakete schlägt um.
        Shot { centers: true, ..shot("10-ohne-finnen-umgeschlagen", |s, _| s.time() >= 22.0, 45.0, 0.15, -0.9) },
    ]
}

/// Schiff in eine Kreisbahn 30 km über Muz setzen (antriebslos).
fn orbit_muz(s: &mut Sim) {
    let muz = s.universe.find("Muz").expect("Muz fehlt");
    let body = s.universe.body(muz);
    let r = body.radius + 30_000.0;
    let v = (body.mu / r).sqrt();
    let t = s.time();
    s.autopilot = None;
    s.clamped = false;
    s.input.throttle = 0.0;
    let st = &mut s.objects[0].state;
    st.body = muz;
    st.time = t;
    st.pos = sim_core::Vec3::new(r, 0.0, 0.0);
    st.vel = sim_core::Vec3::new(0.0, v, 0.0);
    st.omega = sim_core::Vec3::ZERO;
    st.contact_force = 0.0;
    s.objects[0].last = None;
    s.objects.truncate(1);
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
        apps.main.world_mut().resource_mut::<crate::render::ShowCenters>().0 = shot.centers;
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
        println!("{path}  T+{:.1}s  Höhe {:.1} km", sim.time(), sim.ship().state.altitude(sim.body()) / 1000.0);
    }
}
