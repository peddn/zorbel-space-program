//! Schiff als Baum von Teilen.
//!
//! Invariante: der Index eines Elternteils ist immer kleiner als der seiner
//! Kinder (Teil 0 ist die Wurzel). Damit lassen sich Teilbäume in einem
//! einzigen Durchlauf von hinten nach vorn aufsummieren.

use crate::part::{AttachKind, JointStrength, PartDef};
use glam::DMat3;
use sim_core::{PartLocal, Rotation, Vec3, VesselFrame};
use crate::geometry::HullGeometry;
use std::sync::{Arc, OnceLock};

/// Verbindung eines Teils zu seinem Elternteil.
#[derive(Debug, Clone, Copy)]
pub struct Joint {
    /// Verbindungspunkt im Schiffssystem.
    pub point: Vec3<VesselFrame>,
    /// Einheitsvektor vom Eltern- zum Kindteil.
    pub axis: Vec3<VesselFrame>,
    pub strength: JointStrength,
}

#[derive(Debug, Clone)]
pub struct PartInstance {
    pub def: Arc<PartDef>,
    /// Geometrische Mitte des Teils im Schiffssystem.
    pub origin: Vec3<VesselFrame>,
    pub orientation: Rotation<PartLocal, VesselFrame>,
    pub parent: Option<usize>,
    pub joint: Option<Joint>,
    /// Aktueller Treibstoff in kg.
    pub fuel: f64,
    /// Stufe, in der das Teil aktiviert wird (höchste zuerst, wie in KSP).
    pub stage: Option<u32>,
    pub engine_active: bool,
}

impl PartInstance {
    pub fn mass(&self) -> f64 {
        self.def.dry_mass + self.fuel
    }
    pub fn com(&self) -> Vec3<VesselFrame> {
        self.origin + self.orientation.apply(self.def.hull.com_local())
    }
    /// Eigenträgheit um den eigenen Schwerpunkt, in Schiffsachsen.
    pub fn inertia(&self) -> DMat3 {
        let r = self.orientation.matrix();
        r * self.def.hull.inertia_local(self.mass()) * r.transpose()
    }
}

#[derive(Debug, Clone)]
pub struct Vessel {
    pub name: String,
    pub parts: Vec<PartInstance>,
    /// Physikgeometrie (Polygone, Verdeckung, Spannweiten), einmal pro
    /// Teileanordnung berechnet. Nach Änderungen an Lage oder Hülle der
    /// Teile [`Vessel::invalidate_geometry`] aufrufen.
    geometry: OnceLock<Arc<HullGeometry>>,
}

impl Vessel {
    pub fn new(name: impl Into<String>, parts: Vec<PartInstance>) -> Self {
        Self { name: name.into(), parts, geometry: OnceLock::new() }
    }

    pub fn geometry(&self) -> Arc<HullGeometry> {
        self.geometry.get_or_init(|| Arc::new(HullGeometry::build(self))).clone()
    }

    pub fn invalidate_geometry(&mut self) {
        self.geometry = OnceLock::new();
    }

    pub fn len(&self) -> usize {
        self.parts.len()
    }
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    pub fn children(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.parts
            .iter()
            .enumerate()
            .skip(i + 1)
            .filter(move |(_, p)| p.parent == Some(i))
            .map(|(j, _)| j)
    }

    /// Alle Teile im Teilbaum unter `i` (inklusive), aufsteigend sortiert.
    pub fn subtree(&self, i: usize) -> Vec<usize> {
        let mut inside = vec![false; self.parts.len()];
        inside[i] = true;
        for j in i + 1..self.parts.len() {
            if let Some(p) = self.parts[j].parent {
                inside[j] = inside[p];
            }
        }
        (0..self.parts.len()).filter(|&j| inside[j]).collect()
    }

    /// Trennt den Teilbaum unter `child` ab. Rückgabe: (Rest, Abgetrenntes).
    /// Beide behalten das alte Schiffssystem; Positionen bleiben gültig.
    pub fn split_off(&self, child: usize) -> (Vessel, Vessel) {
        assert!(child != 0, "die Wurzel kann nicht abgetrennt werden");
        let sub = self.subtree(child);
        let mut in_sub = vec![false; self.parts.len()];
        for &j in &sub {
            in_sub[j] = true;
        }
        let build = |keep: &dyn Fn(usize) -> bool, name: String| {
            let mut remap = vec![usize::MAX; self.parts.len()];
            let mut parts = Vec::new();
            for (j, p) in self.parts.iter().enumerate() {
                if !keep(j) {
                    continue;
                }
                remap[j] = parts.len();
                let mut p = p.clone();
                match p.parent {
                    Some(par) if keep(par) => p.parent = Some(remap[par]),
                    _ => {
                        p.parent = None;
                        p.joint = None;
                    }
                }
                parts.push(p);
            }
            Vessel::new(name, parts)
        };
        let rest = build(&|j| !in_sub[j], self.name.clone());
        let detached = build(&|j| in_sub[j], format!("{} (Trümmer)", self.name));
        (rest, detached)
    }

    /// Nächste noch nicht ausgelöste Stufe.
    pub fn next_stage(&self) -> Option<u32> {
        self.parts.iter().filter_map(|p| p.stage).max()
    }

    /// Löst die nächste Stufe aus: Triebwerke zünden. Gibt die Indizes
    /// der auszulösenden Entkoppler zurück (Trennung macht der Aufrufer,
    /// weil sie zwei neue Körper mit eigenem Bewegungszustand erzeugt).
    pub fn activate_next_stage(&mut self) -> Vec<usize> {
        let Some(stage) = self.next_stage() else {
            return Vec::new();
        };
        let mut decouplers = Vec::new();
        for (i, p) in self.parts.iter_mut().enumerate() {
            if p.stage == Some(stage) {
                p.stage = None;
                if p.def.engine().is_some() {
                    p.engine_active = true;
                }
                if p.def.is_decoupler() {
                    decouplers.push(i);
                }
            }
        }
        decouplers
    }

    /// Tanks, aus denen ein Triebwerk ziehen kann: alles, was über
    /// Verbindungen erreichbar ist, ohne einen Entkoppler zu durchqueren.
    pub fn fuel_sources(&self, engine: usize) -> Vec<usize> {
        let n = self.parts.len();
        let mut seen = vec![false; n];
        let mut stack = vec![engine];
        seen[engine] = true;
        let mut tanks = Vec::new();
        while let Some(i) = stack.pop() {
            if self.parts[i].def.fuel_capacity() > 0.0 {
                tanks.push(i);
            }
            let neighbours = self.parts[i].parent.into_iter().chain(self.children(i));
            for j in neighbours.collect::<Vec<_>>() {
                if !seen[j] && !self.parts[j].def.is_decoupler() {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        tanks.sort_unstable();
        tanks
    }

    pub fn total_mass(&self) -> f64 {
        self.parts.iter().map(PartInstance::mass).sum()
    }

    /// Lesbarer Name einer Verbindung, z. B. "Tank → Triebwerk".
    pub fn joint_label(&self, child: usize) -> String {
        let c = &self.parts[child];
        match c.parent {
            Some(p) => format!("{} #{p} → {} #{child}", self.parts[p].def.name, c.def.name),
            None => format!("{} #{child} (Wurzel)", c.def.name),
        }
    }
}

/// Baut Schiffe durch Anhängen von Teilen.
pub struct VesselBuilder {
    name: String,
    parts: Vec<PartInstance>,
}

impl VesselBuilder {
    /// Wurzelteil mit Mitte im Ursprung des Schiffssystems.
    pub fn new(name: impl Into<String>, root: Arc<PartDef>) -> Self {
        let fuel = root.fuel_capacity();
        Self {
            name: name.into(),
            parts: vec![PartInstance {
                def: root,
                origin: Vec3::ZERO,
                orientation: Rotation::IDENTITY,
                parent: None,
                joint: None,
                fuel,
                stage: None,
                engine_active: false,
            }],
        }
    }

    fn push(&mut self, def: Arc<PartDef>, origin: Vec3<VesselFrame>, orientation: Rotation<PartLocal, VesselFrame>, parent: usize, point: Vec3<VesselFrame>, axis: Vec3<VesselFrame>) -> usize {
        let strength = def.node_strength.weakest(&self.parts[parent].def.node_strength);
        let fuel = def.fuel_capacity();
        self.parts.push(PartInstance {
            def,
            origin,
            orientation,
            parent: Some(parent),
            joint: Some(Joint { point, axis, strength }),
            fuel,
            stage: None,
            engine_active: false,
        });
        self.parts.len() - 1
    }

    /// Hängt ein Stapelteil unter `parent` an.
    pub fn stack_below(&mut self, parent: usize, def: Arc<PartDef>) -> usize {
        let p = &self.parts[parent];
        assert_eq!(p.def.attach, AttachKind::Stack, "{} hat keinen unteren Knoten", p.def.name);
        assert_eq!(def.attach, AttachKind::Stack, "{} ist kein Stapelteil", def.name);
        let down = Vec3::new(0.0, 0.0, -1.0);
        let node = p.origin + down * (p.def.hull.stack_length() / 2.0);
        let origin = node + down * (def.hull.stack_length() / 2.0);
        self.push(def, origin, Rotation::IDENTITY, parent, node, down)
    }

    /// Befestigt ein Radialteil auf der Mantelfläche von `parent`
    /// unter dem Winkel `angle` (rad, um die Längsachse) und mit
    /// Höhenversatz `z_offset` relativ zur Mitte des Elternteils.
    pub fn radial(&mut self, parent: usize, def: Arc<PartDef>, angle: f64, z_offset: f64) -> usize {
        let p = &self.parts[parent];
        assert_eq!(def.attach, AttachKind::Radial, "{} ist kein Radialteil", def.name);
        let outward = Vec3::new(angle.cos(), angle.sin(), 0.0);
        let node = p.origin + outward * p.def.hull.surface_radius() + Vec3::new(0.0, 0.0, z_offset);
        // Radialteile zeigen mit ihrer lokalen +X-Achse nach außen.
        let half_span = match def.hull {
            crate::part::Shape::Plate { span, .. } => span / 2.0,
            s => s.surface_radius(),
        };
        let origin = node + outward * half_span;
        let orientation = Rotation::from_quat(glam::DQuat::from_rotation_z(angle));
        self.push(def, origin, orientation, parent, node, outward)
    }

    /// `count` Kopien gleichmäßig verteilt (Symmetrie wie im KSP-Editor).
    pub fn radial_symmetric(&mut self, parent: usize, def: Arc<PartDef>, count: usize, z_offset: f64) -> Vec<usize> {
        (0..count)
            .map(|k| {
                let angle = std::f64::consts::TAU * k as f64 / count as f64;
                self.radial(parent, def.clone(), angle, z_offset)
            })
            .collect()
    }

    pub fn set_stage(&mut self, part: usize, stage: u32) -> &mut Self {
        self.parts[part].stage = Some(stage);
        self
    }

    pub fn build(self) -> Vessel {
        Vessel::new(self.name, self.parts)
    }
}
