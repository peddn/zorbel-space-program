//! Darstellung: kamerarelatives Rendering aus der f64-Simulation.
//!
//! Die Simulation rechnet in f64 im Inertialsystem. Für die Grafik wird in
//! jedem Bild alles relativ zum Kamerafokus (dem aktiven Schiff) in f32
//! umgerechnet. Der Fokus liegt immer im Ursprung der Grafikwelt; einen
//! "Floating Origin", der die Welt in Sprüngen verschiebt, gibt es nicht.

use crate::sim::Sim;
use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use sim_core::{BodyFixed, CelestialBody, Inertial, Rotation as SimRotation, Vec3 as SimVec3};
use std::f32::consts::FRAC_PI_2;
use vessel::Vessel;
use vessel::mass::mass_properties;
use vessel::part::Shape;

// ---------- Koordinaten ----------

/// Inertialsystem (Z = Drehachse) → Grafikwelt (Y = oben):
/// Drehung um −90° um X, also (x, y, z) ↦ (x, z, −y).
pub fn render_vec(v: SimVec3<Inertial>) -> Vec3 {
    Vec3::new(v.x() as f32, v.z() as f32, -v.y() as f32)
}

fn quat_f32(q: glam::DQuat) -> Quat {
    Quat::from_xyzw(q.x as f32, q.y as f32, q.z as f32, q.w as f32)
}

/// Lage eines Systems, dessen Achsen per `rot` ins Inertialsystem zeigen.
pub fn render_rot<F: sim_core::Frame>(rot: SimRotation<F, Inertial>) -> Quat {
    Quat::from_rotation_x(-FRAC_PI_2) * quat_f32(rot.quat())
}

/// Inertiale Position des Kamerafokus.
#[derive(Resource, Default, Clone, Copy)]
pub struct Focus(pub SimVec3<Inertial>);

pub fn update_focus(sim: Res<Sim>, mut focus: ResMut<Focus>) {
    focus.0 = sim.ship().state.pos;
}

// ---------- Planet ----------

#[derive(Component)]
pub struct PlanetRoot;
#[derive(Component)]
pub struct AtmosphereShell;
#[derive(Component)]
pub struct LaunchPad;
#[derive(Component)]
pub struct Sun;

/// Wertrauschen auf der Kugel (nahtlos, weil im 3D-Raum ausgewertet).
fn hash(i: i32, j: i32, k: i32) -> f32 {
    let mut h = (i as u32).wrapping_mul(374_761_393) ^ (j as u32).wrapping_mul(668_265_263) ^ (k as u32).wrapping_mul(2_147_483_647);
    h = (h ^ (h >> 13)).wrapping_mul(1_274_126_177);
    ((h ^ (h >> 16)) & 0xffff) as f32 / 65_535.0
}

fn value_noise(p: Vec3) -> f32 {
    let f = p.floor();
    let t = p - f;
    let s = t * t * (Vec3::splat(3.0) - 2.0 * t);
    let (i, j, k) = (f.x as i32, f.y as i32, f.z as i32);
    let mut acc = 0.0;
    for (dx, wx) in [(0, 1.0 - s.x), (1, s.x)] {
        for (dy, wy) in [(0, 1.0 - s.y), (1, s.y)] {
            for (dz, wz) in [(0, 1.0 - s.z), (1, s.z)] {
                acc += wx * wy * wz * hash(i + dx, j + dy, k + dz);
            }
        }
    }
    acc
}

/// Oberflächenfarbe an einem Punkt (Einheitsvektor im körperfesten System).
pub fn surface_color(dir: Vec3) -> [f32; 4] {
    let mut n = 0.0;
    let mut amp = 0.5;
    let mut freq = 2.0;
    for _ in 0..6 {
        n += amp * value_noise(dir * freq + Vec3::splat(11.0)); // Startplatz liegt im Grasland
        amp *= 0.5;
        freq *= 2.1;
    }
    let lat = dir.z.abs();
    let (r, g, b) = if lat > 0.93 {
        (0.92, 0.94, 0.97) // Eiskappen
    } else if n < 0.47 {
        let d = (0.47 - n) * 2.0;
        (0.05, 0.16 + 0.1 * (1.0 - d), 0.38 + 0.15 * (1.0 - d)) // Ozean
    } else if n < 0.6 {
        (0.22, 0.42 + (n - 0.47), 0.16) // Grasland
    } else {
        (0.45, 0.38, 0.26) // Hochland
    };
    [r, g, b, 1.0]
}

fn planet_texture() -> Image {
    let (w, h) = (1024u32, 512u32);
    let mut data = Vec::with_capacity((w * h * 4) as usize);
    for y in 0..h {
        // Bevys UV-Kugel: v = 0 am Nordpol (+Z), u = 0 bei Länge 0 (+X)
        let lat = std::f32::consts::PI * (0.5 - (y as f32 + 0.5) / h as f32);
        for x in 0..w {
            let lon = std::f32::consts::TAU * (x as f32 + 0.5) / w as f32;
            let dir = Vec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
            let c = surface_color(dir);
            data.extend(c.iter().map(|v| (v.clamp(0.0, 1.0) * 255.0) as u8));
        }
    }
    Image::new(
        Extent3d { width: w, height: h, depth_or_array_layers: 1 },
        TextureDimension::D2,
        data,
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::RENDER_WORLD,
    )
}

/// Die große Kugel liegt etwas tiefer als die wahre Oberfläche; in
/// Kameranähe liegt das genaue Bodenstück darüber (siehe `GroundPatch`).
const SPHERE_SINK: f32 = 60.0;

pub fn setup_planet(
    mut commands: Commands,
    sim: Res<Sim>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut images: ResMut<Assets<Image>>,
) {
    let r = sim.body.radius as f32;
    let texture = images.add(planet_texture());
    commands
        .spawn((PlanetRoot, Transform::default(), Visibility::default()))
        .with_children(|p| {
            p.spawn((
                Mesh3d(meshes.add(Sphere::new(r - SPHERE_SINK).mesh().uv(256, 128))),
                MeshMaterial3d(materials.add(StandardMaterial {
                    base_color_texture: Some(texture),
                    perceptual_roughness: 0.9,
                    ..default()
                })),
                // Bevys UV-Kugel hat ihre Pole auf ±Z und u = Länge ab +X:
                // genau das körperfeste System, keine Drehung nötig.
                Transform::default(),
            ));
        });

    let top = sim.body.atmosphere.top() as f32;
    commands.spawn((
        AtmosphereShell,
        Mesh3d(meshes.add(Sphere::new(r + top).mesh().uv(128, 64))),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::srgba(0.45, 0.65, 1.0, 0.05),
            alpha_mode: AlphaMode::Blend,
            unlit: true,
            ..default()
        })),
        Transform::default(),
    ));

    commands.spawn((
        LaunchPad,
        Mesh3d(meshes.add(Cuboid::new(14.0, 1.0, 14.0))),
        MeshMaterial3d(materials.add(Color::srgb(0.55, 0.55, 0.55))),
        Transform::default(),
    ));

    commands.spawn((
        Sun,
        DirectionalLight { illuminance: 12_000.0, shadow_maps_enabled: false, ..default() },
        Transform::default(),
    ));
    commands.insert_resource(GlobalAmbientLight { brightness: 600.0, ..default() });
}

/// Richtung zur Sonne im Inertialsystem (fest; der Startplatz liegt bei
/// t = 0 auf der Tagseite).
fn sun_direction() -> SimVec3<Inertial> {
    SimVec3::new(0.8, -0.45, 0.3).normalize_or_zero()
}

pub fn sync_planet(
    sim: Res<Sim>,
    focus: Res<Focus>,
    mut planet: Query<&mut Transform, (With<PlanetRoot>, Without<AtmosphereShell>, Without<LaunchPad>, Without<Sun>)>,
    mut shell: Query<&mut Transform, (With<AtmosphereShell>, Without<LaunchPad>, Without<Sun>)>,
    mut pad: Query<&mut Transform, (With<LaunchPad>, Without<Sun>)>,
    mut sun: Query<&mut Transform, With<Sun>>,
    mut clear: ResMut<ClearColor>,
) {
    let t = sim.time();
    let body = &sim.body;
    let center = render_vec(-focus.0);
    let orient = body.orientation(t);
    for mut tr in &mut planet {
        tr.translation = center;
        tr.rotation = render_rot(orient);
    }
    for mut tr in &mut shell {
        tr.translation = center;
    }
    for mut tr in &mut pad {
        let site = body.geo_to_body_fixed(vessel::scenario::LAUNCH_SITE);
        let up = site.normalize_or_zero();
        let pos = orient.apply(site - up * 0.5);
        tr.translation = render_vec(pos - focus.0);
        // Würfel-Y → radial nach außen
        let local = glam::DQuat::from_rotation_arc(glam::DVec3::Y, up.raw());
        tr.rotation = render_rot(orient) * quat_f32(local);
    }
    for mut tr in &mut sun {
        *tr = Transform::default().looking_to(-render_vec(sun_direction()), Vec3::Y);
    }
    // Himmel: blau in dichter Luft, schwarz im Vakuum.
    let alt = sim.ship().state.altitude(body);
    let p = body.atmosphere.sample(alt).pressure / body.atmosphere.surface_pressure.max(1.0);
    let k = (p.max(0.0)).powf(0.25) as f32;
    clear.0 = Color::srgb(0.01 + 0.36 * k, 0.01 + 0.56 * k, 0.03 + 0.85 * k);
}

// ---------- Genaues Bodenstück unter der Kamera ----------

#[derive(Component)]
pub struct GroundPatch {
    /// Mittelpunkt (Einheitsvektor, körperfest), für den das Mesh gebaut ist.
    center: SimVec3<BodyFixed>,
}

const PATCH_HALF: f64 = 15_000.0;
const PATCH_N: u32 = 96;
const PATCH_MAX_ALT: f64 = 60_000.0;

/// Gitter auf der exakten Kugel, Koordinaten relativ zum Mittelpunkt
/// (darum in f32 zentimetergenau).
fn patch_mesh(body: &CelestialBody, center: SimVec3<BodyFixed>) -> Mesh {
    let r = body.radius;
    let c = center * r;
    let east = SimVec3::<BodyFixed>::Z.cross(center).normalize_or_zero();
    let east = if east == SimVec3::ZERO { SimVec3::X } else { east };
    let north = center.cross(east);
    let mut pos = Vec::new();
    let mut nrm = Vec::new();
    let mut col = Vec::new();
    for j in 0..=PATCH_N {
        for i in 0..=PATCH_N {
            let u = (i as f64 / PATCH_N as f64 * 2.0 - 1.0) * PATCH_HALF;
            let v = (j as f64 / PATCH_N as f64 * 2.0 - 1.0) * PATCH_HALF;
            let dir = (center * r + east * u + north * v).normalize_or_zero();
            let p = dir * r - c;
            pos.push([p.x() as f32, p.y() as f32, p.z() as f32]);
            nrm.push([dir.x() as f32, dir.y() as f32, dir.z() as f32]);
            // Vertexfarben sind linear, die Planetentextur ist sRGB.
            let c = surface_color(Vec3::new(dir.x() as f32, dir.y() as f32, dir.z() as f32));
            let lin = Color::srgb(c[0], c[1], c[2]).to_linear();
            col.push([lin.red, lin.green, lin.blue, 1.0]);
        }
    }
    let mut idx = Vec::new();
    let w = PATCH_N + 1;
    for j in 0..PATCH_N {
        for i in 0..PATCH_N {
            let a = j * w + i;
            idx.extend([a, a + 1, a + w, a + 1, a + w + 1, a + w]);
        }
    }
    Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD)
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, pos)
        .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, nrm)
        .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, col)
        .with_inserted_indices(Indices::U32(idx))
}

pub fn sync_ground_patch(
    mut commands: Commands,
    sim: Res<Sim>,
    focus: Res<Focus>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut patch: Query<(Entity, &GroundPatch, &mut Transform, &mut Visibility)>,
) {
    let body = &sim.body;
    let t = sim.time();
    let orient = body.orientation(t);
    let below = body.to_body_fixed(focus.0, t).normalize_or_zero();
    let alt = focus.0.length() - body.radius;

    let needs_new = match patch.single() {
        Ok((_, p, _, _)) => (p.center - below).length() * body.radius > PATCH_HALF * 0.4,
        Err(_) => true,
    };
    if needs_new && alt < PATCH_MAX_ALT {
        if let Ok((e, _, _, _)) = patch.single() {
            commands.entity(e).despawn();
        }
        commands.spawn((
            GroundPatch { center: below },
            Mesh3d(meshes.add(patch_mesh(body, below))),
            MeshMaterial3d(materials.add(StandardMaterial { perceptual_roughness: 0.95, ..default() })),
            Transform::default(),
            Visibility::default(),
        ));
        return;
    }
    if let Ok((_, p, mut tr, mut vis)) = patch.single_mut() {
        let c = orient.apply(p.center * body.radius);
        tr.translation = render_vec(c - focus.0);
        tr.rotation = render_rot(orient);
        *vis = if alt < PATCH_MAX_ALT { Visibility::Inherited } else { Visibility::Hidden };
    }
}

// ---------- Schiffe ----------

#[derive(Component)]
pub struct VesselView {
    generation: u64,
    id: u64,
    version: u64,
}

#[derive(Component)]
pub struct PartView(usize);

#[derive(Component)]
pub struct Flame(usize);

/// Materialien für die Teile, einmal angelegt.
#[derive(Resource)]
pub struct PartMaterials {
    capsule: Handle<StandardMaterial>,
    tank: Handle<StandardMaterial>,
    engine: Handle<StandardMaterial>,
    decoupler: Handle<StandardMaterial>,
    fin: Handle<StandardMaterial>,
    flame: Handle<StandardMaterial>,
}

pub fn setup_part_materials(mut commands: Commands, mut materials: ResMut<Assets<StandardMaterial>>) {
    let mut m = |c: Color, rough: f32, metal: f32| {
        materials.add(StandardMaterial { base_color: c, perceptual_roughness: rough, metallic: metal, ..default() })
    };
    let capsule = m(Color::srgb(0.75, 0.76, 0.78), 0.5, 0.6);
    let tank = m(Color::srgb(0.93, 0.92, 0.88), 0.6, 0.1);
    let engine = m(Color::srgb(0.25, 0.25, 0.27), 0.4, 0.8);
    let decoupler = m(Color::srgb(0.2, 0.2, 0.22), 0.7, 0.3);
    let fin = m(Color::srgb(0.85, 0.3, 0.12), 0.6, 0.1);
    let flame = materials.add(StandardMaterial {
        base_color: Color::srgba(1.0, 0.7, 0.3, 0.8),
        emissive: LinearRgba::rgb(8.0, 4.0, 1.2),
        unlit: true,
        alpha_mode: AlphaMode::Add,
        ..default()
    });
    commands.insert_resource(PartMaterials { capsule, tank, engine, decoupler, fin, flame });
}

/// Mesh und Achskorrektur eines Teils. Bevys Rotationskörper liegen auf
/// +Y, unsere Teile auf +Z.
fn part_mesh(shape: Shape) -> (Mesh, Quat) {
    let y_to_z = Quat::from_rotation_x(FRAC_PI_2);
    match shape {
        Shape::Cylinder { radius, length } => (Cylinder::new(radius as f32, length as f32).mesh().resolution(32).into(), y_to_z),
        Shape::Cone { radius_bottom, radius_top, length } => (
            ConicalFrustum { radius_top: radius_top as f32, radius_bottom: radius_bottom as f32, height: length as f32 }
                .mesh()
                .resolution(32)
                .into(),
            y_to_z,
        ),
        Shape::Plate { span, chord, thickness } => (Cuboid::new(span as f32, thickness as f32, chord as f32).into(), Quat::IDENTITY),
    }
}

#[allow(clippy::too_many_arguments)]
fn spawn_vessel(
    commands: &mut Commands,
    meshes: &mut Assets<Mesh>,
    mats: &PartMaterials,
    generation: u64,
    id: u64,
    version: u64,
    v: &Vessel,
) {
    commands
        .spawn((VesselView { generation, id, version }, Transform::default(), Visibility::default()))
        .with_children(|root| {
            for (i, p) in v.parts.iter().enumerate() {
                let (mesh, fix) = part_mesh(p.def.hull);
                let mat = if p.def.engine().is_some() {
                    &mats.engine
                } else if p.def.is_decoupler() {
                    &mats.decoupler
                } else if p.def.command_torque() > 0.0 {
                    &mats.capsule
                } else if p.def.fuel_capacity() > 0.0 {
                    &mats.tank
                } else {
                    &mats.fin
                };
                root.spawn((PartView(i), Mesh3d(meshes.add(mesh)), MeshMaterial3d(mat.clone()), Transform::from_rotation(fix)));
                if let (Some(_), Shape::Cone { radius_bottom, .. }) = (p.def.engine(), p.def.hull) {
                    let flame = ConicalFrustum { radius_top: radius_bottom as f32 * 0.85, radius_bottom: 0.08, height: 1.0 };
                    root.spawn((
                        Flame(i),
                        Mesh3d(meshes.add(flame.mesh().resolution(24))),
                        MeshMaterial3d(mats.flame.clone()),
                        Transform::default(),
                        Visibility::Hidden,
                    ));
                }
            }
        });
}

#[allow(clippy::too_many_arguments)]
pub fn sync_vessels(
    mut commands: Commands,
    sim: Res<Sim>,
    focus: Res<Focus>,
    mats: Res<PartMaterials>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut roots: Query<(Entity, &VesselView, &mut Transform, &Children)>,
    mut parts: Query<(&PartView, &mut Transform), Without<VesselView>>,
    mut flames: Query<(&Flame, &mut Transform, &mut Visibility), (Without<VesselView>, Without<PartView>)>,
) {
    // Veraltete oder verschwundene Ansichten entfernen
    let current = |view: &VesselView| {
        sim.objects.iter().find(|o| view.generation == sim.generation && o.id == view.id && o.version == view.version)
    };
    let mut present = Vec::new();
    for (e, view, _, _) in &roots {
        match current(view) {
            Some(o) => present.push(o.id),
            None => commands.entity(e).despawn(),
        }
    }
    for o in &sim.objects {
        if !present.contains(&o.id) {
            spawn_vessel(&mut commands, &mut meshes, &mats, sim.generation, o.id, o.version, &o.vessel);
        }
    }

    for (_, view, mut tr, children) in &mut roots {
        let Some(o) = current(view) else { continue };
        tr.translation = render_vec(o.state.pos - focus.0);
        tr.rotation = render_rot(o.state.attitude);
        let com = mass_properties(&o.vessel).com;
        let thrust = o.last.as_ref().map_or(0.0, |r| r.thrust);
        let max_thrust: f64 = o.vessel.parts.iter().filter(|p| p.engine_active).filter_map(|p| p.def.engine()).map(|e| e.thrust_vac).sum();
        let frac = if max_thrust > 0.0 { (thrust / max_thrust).clamp(0.0, 1.2) } else { 0.0 } as f32;
        for child in children.iter() {
            if let Ok((pv, mut ptr)) = parts.get_mut(child) {
                let p = &o.vessel.parts[pv.0];
                let d = p.origin - com;
                ptr.translation = Vec3::new(d.x() as f32, d.y() as f32, d.z() as f32);
                let fix = part_mesh_fix(p.def.hull);
                ptr.rotation = quat_f32(p.orientation.quat()) * fix;
            }
            if let Ok((fl, mut ftr, mut vis)) = flames.get_mut(child) {
                let p = &o.vessel.parts[fl.0];
                let len = p.def.hull.stack_length();
                let flame_len = 1.5 + 9.0 * frac;
                let d = p.origin - com;
                ftr.translation = Vec3::new(d.x() as f32, d.y() as f32, (d.z() - len / 2.0) as f32 - flame_len / 2.0);
                // Flamme zeigt nach −Z: breites Ende oben am Triebwerk
                ftr.rotation = Quat::from_rotation_x(FRAC_PI_2);
                ftr.scale = Vec3::new(1.0, flame_len, 1.0);
                *vis = if p.engine_active && frac > 0.01 { Visibility::Inherited } else { Visibility::Hidden };
            }
        }
    }
}

fn part_mesh_fix(shape: Shape) -> Quat {
    match shape {
        Shape::Plate { .. } => Quat::IDENTITY,
        _ => Quat::from_rotation_x(FRAC_PI_2),
    }
}

// ---------- Bahnlinie ----------

pub fn draw_orbits(sim: Res<Sim>, focus: Res<Focus>, cam: Res<crate::camera::OrbitCamera>, mut gizmos: Gizmos) {
    // Wie in KSP: Bahnlinien nur in der Übersicht, nicht in der Nahansicht.
    if cam.distance < 5_000.0 {
        return;
    }
    let body = &sim.body;
    for (k, o) in sim.objects.iter().enumerate() {
        let s = o.state.orbit();
        let el = sim_core::kepler::elements(s, body.mu);
        // Bis zum Aufschlag oder einen Umlauf (offene Bahn: 2 h)
        let span = el.period.unwrap_or(7_200.0);
        let span = sim_core::kepler::time_until_below(s, body.mu, body.radius, span).unwrap_or(span);
        let n = 256;
        let pts = (0..=n).map(|i| {
            let st = sim_core::kepler::propagate(s, body.mu, span * i as f64 / n as f64);
            render_vec(st.r - focus.0)
        });
        let color = if k == 0 { Color::srgb(0.3, 0.9, 1.0) } else { Color::srgba(0.8, 0.6, 0.4, 0.5) };
        gizmos.linestrip(pts, color);
    }
}


// ---------- Schwerpunkt und Druckpunkt ----------

/// Markierungen für Schwer- und Druckpunkt anzeigen (Taste C).
#[derive(Resource, Default)]
pub struct ShowCenters(pub bool);

/// Eigene Gizmo-Gruppe: wird über allem gezeichnet, auch durch die Rakete.
#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct MarkerGizmos;

pub fn setup_marker_gizmos(mut store: ResMut<GizmoConfigStore>) {
    let (config, _) = store.config_mut::<MarkerGizmos>();
    config.depth_bias = -1.0;
    config.line.width = 3.0;
}

pub struct Centers {
    pub com: SimVec3<Inertial>,
    pub cop: SimVec3<Inertial>,
    /// Abstand Druckpunkt hinter Schwerpunkt (m); positiv = stabil.
    pub margin: f64,
    /// Anstellwinkel über 90°: das Schiff fliegt mit dem Heck voran.
    pub tail_first: bool,
}

/// Schwer- und Druckpunkt des aktiven Schiffs. Der Druckpunkt gilt für den
/// aktuellen Anstellwinkel (mindestens 3°); im Vakuum wird mit Luft auf
/// Meereshöhe gerechnet, damit man die Rakete trotzdem beurteilen kann.
pub fn centers(sim: &Sim) -> Option<Centers> {
    let ship = sim.ship();
    let s = &ship.state;
    let body = &sim.body;
    let mut atm = body.atmosphere.sample(s.altitude(body));
    if atm.density <= 0.0 {
        atm = body.atmosphere.sample(0.0);
    }
    let surf = s.surface_velocity(body);
    let speed = surf.length().max(50.0);
    let aoa = if surf.length() > 1.0 { s.nose().dot(surf / surf.length()).clamp(-1.0, 1.0).acos() } else { 0.0 };
    let alpha = aoa.clamp(3f64.to_radians(), 60f64.to_radians());
    let com = mass_properties(&ship.vessel).com;
    let z = vessel::aero::center_of_pressure(&ship.vessel, &atm, speed, alpha)?;
    let cop = s.pos + s.attitude.apply(sim_core::Vec3::new(0.0, 0.0, z - com.z()));
    Some(Centers { com: s.pos, cop, margin: com.z() - z, tail_first: aoa > std::f64::consts::FRAC_PI_2 })
}

pub fn draw_centers(sim: Res<Sim>, focus: Res<Focus>, show: Res<ShowCenters>, mut gizmos: Gizmos<MarkerGizmos>) {
    if !show.0 {
        return;
    }
    let Some(c) = centers(&sim) else { return };
    let com = render_vec(c.com - focus.0);
    let cop = render_vec(c.cop - focus.0);
    let yellow = Color::srgb(1.0, 0.85, 0.1);
    let cyan = Color::srgb(0.2, 0.8, 1.0);
    gizmos.sphere(Isometry3d::from_translation(com), 0.45, yellow);
    gizmos.sphere(Isometry3d::from_translation(cop), 0.45, cyan);
    gizmos.line(com, cop, if c.margin > 0.0 { Color::srgb(0.3, 1.0, 0.3) } else { Color::srgb(1.0, 0.25, 0.2) });
}
