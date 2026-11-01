//! Bauteil-Definitionen. Werden aus RON-Dateien geladen, damit neue Teile
//! ohne Neukompilieren hinzukommen können (Grundlage für Modding).

use glam::DMat3;
use serde::Deserialize;
use sim_core::{PartLocal, Vec3};
use std::collections::BTreeMap;
use std::fmt;
use std::path::Path;
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
pub enum AttachKind {
    /// Oben/unten auf der Längsachse gestapelt.
    Stack,
    /// Seitlich auf der Mantelfläche eines Stapelteils.
    Radial,
}

/// Geometrie, aus der Schwerpunkt, Trägheit und Knotenpositionen folgen.
/// Lokale Achsen: +Z = Längsachse (zur Spitze).
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
pub enum Shape {
    Cylinder { radius: f64, length: f64 },
    /// Kegelstumpf, unten breiter oder schmaler als oben.
    Cone { radius_bottom: f64, radius_top: f64, length: f64 },
    /// Flache Platte: X = Spannweite, Z = Profiltiefe, Y = Dicke.
    Plate { span: f64, chord: f64, thickness: f64 },
}

impl Shape {
    /// Länge entlang der Stapelachse.
    pub fn stack_length(&self) -> f64 {
        match *self {
            Shape::Cylinder { length, .. } | Shape::Cone { length, .. } => length,
            Shape::Plate { chord, .. } => chord,
        }
    }

    /// Schwerpunkt relativ zur geometrischen Mitte (homogene Dichte).
    pub fn com_local(&self) -> Vec3<PartLocal> {
        match *self {
            Shape::Cone { radius_bottom: r1, radius_top: r2, length: h } => {
                // Schwerpunkt eines Kegelstumpfs, gemessen von der Unterseite
                let from_bottom =
                    h * (r1 * r1 + 2.0 * r1 * r2 + 3.0 * r2 * r2) / (4.0 * (r1 * r1 + r1 * r2 + r2 * r2));
                Vec3::new(0.0, 0.0, from_bottom - h / 2.0)
            }
            _ => Vec3::ZERO,
        }
    }

    /// Trägheitstensor um den eigenen Schwerpunkt, in lokalen Achsen.
    pub fn inertia_local(&self, mass: f64) -> DMat3 {
        let (ixx, iyy, izz) = match *self {
            Shape::Cylinder { radius: r, length: l } => {
                let t = mass * (3.0 * r * r + l * l) / 12.0;
                (t, t, 0.5 * mass * r * r)
            }
            Shape::Cone { radius_bottom: r1, radius_top: r2, length: l } => {
                // Näherung: Zylinder mit volumengleichem Radius.
                let r2eq = (r1 * r1 + r1 * r2 + r2 * r2) / 3.0;
                let t = mass * (3.0 * r2eq + l * l) / 12.0;
                (t, t, 0.5 * mass * r2eq)
            }
            Shape::Plate { span: x, chord: z, thickness: y } => (
                mass * (y * y + z * z) / 12.0,
                mass * (x * x + z * z) / 12.0,
                mass * (x * x + y * y) / 12.0,
            ),
        };
        DMat3::from_diagonal(glam::DVec3::new(ixx, iyy, izz))
    }

    /// Radien der Stirnflächen (unten, oben); `None` für Platten.
    pub fn end_radii(&self) -> Option<(f64, f64)> {
        match *self {
            Shape::Cylinder { radius, .. } => Some((radius, radius)),
            Shape::Cone { radius_bottom, radius_top, .. } => Some((radius_bottom, radius_top)),
            Shape::Plate { .. } => None,
        }
    }

    /// Radius der Mantelfläche für radiale Befestigung.
    pub fn surface_radius(&self) -> f64 {
        match *self {
            Shape::Cylinder { radius, .. } => radius,
            Shape::Cone { radius_bottom, radius_top, .. } => 0.5 * (radius_bottom + radius_top),
            Shape::Plate { thickness, .. } => thickness / 2.0,
        }
    }
}

/// Belastbarkeit einer Verbindung (N bzw. N·m).
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
pub struct JointStrength {
    pub tension: f64,
    pub compression: f64,
    pub shear: f64,
    pub bending: f64,
    pub torsion: f64,
}

impl JointStrength {
    /// Eine Verbindung ist so stark wie ihr schwächerer Partner.
    pub fn weakest(&self, o: &Self) -> Self {
        Self {
            tension: self.tension.min(o.tension),
            compression: self.compression.min(o.compression),
            shear: self.shear.min(o.shear),
            bending: self.bending.min(o.bending),
            torsion: self.torsion.min(o.torsion),
        }
    }
}

/// Funktionsmodule eines Teils. Ein Teil kann mehrere haben.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
pub enum Module {
    /// Schub und Isp im Vakuum; Isp bei 1 atm. Der Massenstrom ist fest,
    /// der Schub skaliert mit dem Isp, also mit dem Umgebungsdruck.
    Engine { thrust_vac: f64, isp_vac: f64, isp_sl: f64, gimbal_deg: f64 },
    FuelTank { capacity: f64 },
    Command { torque: f64 },
    /// Bewegliche Fläche: dreht sich um ihre Spannweitenachse (lokal X).
    /// Ihre Luftkräfte folgen wie bei jedem Teil aus der Hülle.
    ControlSurface { max_deflection_deg: f64 },
    Decoupler { eject_impulse: f64 },
}

/// Referenzdruck für `isp_sl` (1 atm).
pub const P_REF: f64 = 101_325.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EngineSpec {
    pub thrust_vac: f64,
    pub isp_vac: f64,
    pub isp_sl: f64,
    pub gimbal_rad: f64,
}

impl EngineSpec {
    /// Massenstrom in kg/s.
    pub fn mass_flow(&self, throttle: f64) -> f64 {
        self.thrust_vac * throttle / (self.isp_vac * sim_core::G0)
    }
    /// Spezifischer Impuls beim Umgebungsdruck `pressure` (Pa), linear
    /// zwischen Vakuum und 1 atm; Vakuum ist einfach `pressure = 0`.
    pub fn isp(&self, pressure: f64) -> f64 {
        (self.isp_vac + (self.isp_sl - self.isp_vac) * pressure / P_REF).max(0.0)
    }
    pub fn thrust(&self, throttle: f64, pressure: f64) -> f64 {
        self.mass_flow(throttle) * self.isp(pressure) * sim_core::G0
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ControlSurfaceSpec {
    pub max_deflection_rad: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct PartDef {
    pub id: String,
    pub name: String,
    /// Leermasse in kg.
    pub dry_mass: f64,
    pub attach: AttachKind,
    /// Einfaches Polygonmodell für die Physik (Luft, Boden, Verdeckung).
    /// Die Darstellung benutzt eigene Modelle.
    pub hull: Shape,
    pub node_strength: JointStrength,
    /// Größte Aufprallgeschwindigkeit, die das Teil am Boden übersteht (m/s).
    #[serde(default = "default_crash_tolerance")]
    pub crash_tolerance: f64,
    /// Abstimmfaktoren für die Spielbalance. Die Luftkräfte folgen für alle
    /// Teile aus derselben Rechnung; diese Faktoren skalieren sie nur.
    #[serde(default = "one")]
    pub drag_factor: f64,
    #[serde(default = "one")]
    pub lift_factor: f64,
    #[serde(default)]
    pub modules: Vec<Module>,
}

fn default_crash_tolerance() -> f64 {
    8.0
}

fn one() -> f64 {
    1.0
}

impl PartDef {
    pub fn engine(&self) -> Option<EngineSpec> {
        self.modules.iter().find_map(|m| match *m {
            Module::Engine { thrust_vac, isp_vac, isp_sl, gimbal_deg } => Some(EngineSpec {
                thrust_vac,
                isp_vac,
                isp_sl,
                gimbal_rad: gimbal_deg.to_radians(),
            }),
            _ => None,
        })
    }
    pub fn control_surface(&self) -> Option<ControlSurfaceSpec> {
        self.modules.iter().find_map(|m| match *m {
            Module::ControlSurface { max_deflection_deg } => {
                Some(ControlSurfaceSpec { max_deflection_rad: max_deflection_deg.to_radians() })
            }
            _ => None,
        })
    }
    pub fn fuel_capacity(&self) -> f64 {
        self.modules
            .iter()
            .map(|m| match *m {
                Module::FuelTank { capacity } => capacity,
                _ => 0.0,
            })
            .sum()
    }
    pub fn command_torque(&self) -> f64 {
        self.modules
            .iter()
            .map(|m| match *m {
                Module::Command { torque } => torque,
                _ => 0.0,
            })
            .sum()
    }
    pub fn decoupler_impulse(&self) -> Option<f64> {
        self.modules.iter().find_map(|m| match *m {
            Module::Decoupler { eject_impulse } => Some(eject_impulse),
            _ => None,
        })
    }
    pub fn is_decoupler(&self) -> bool {
        self.decoupler_impulse().is_some()
    }
}

#[derive(Debug)]
pub enum LoadError {
    Io(std::io::Error),
    Parse { source: String, error: Box<ron::error::SpannedError> },
    DuplicateId(String),
}

impl fmt::Display for LoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LoadError::Io(e) => write!(f, "E/A-Fehler: {e}"),
            LoadError::Parse { source, error } => write!(f, "{source}: {error}"),
            LoadError::DuplicateId(id) => write!(f, "Teile-ID doppelt vergeben: {id}"),
        }
    }
}
impl std::error::Error for LoadError {}

/// Alle bekannten Teile. `BTreeMap` statt `HashMap`: feste Reihenfolge,
/// damit die Simulation deterministisch bleibt.
#[derive(Debug, Clone, Default)]
pub struct PartLibrary {
    parts: BTreeMap<String, Arc<PartDef>>,
}

impl PartLibrary {
    /// Die mitgelieferten Teile, zur Compile-Zeit eingebettet.
    pub fn builtin() -> Self {
        const SOURCES: [(&str, &str); 5] = [
            ("kapsel.ron", include_str!("../../../assets/parts/kapsel.ron")),
            ("tank.ron", include_str!("../../../assets/parts/tank.ron")),
            ("triebwerk.ron", include_str!("../../../assets/parts/triebwerk.ron")),
            ("steuerflaeche.ron", include_str!("../../../assets/parts/steuerflaeche.ron")),
            ("entkoppler.ron", include_str!("../../../assets/parts/entkoppler.ron")),
        ];
        let mut lib = Self::default();
        for (name, src) in SOURCES {
            lib.add_ron(name, src).expect("eingebettete Teiledefinition fehlerhaft");
        }
        lib
    }

    /// Lädt zusätzlich alle `*.ron` aus einem Verzeichnis (z. B. Mods).
    pub fn load_dir(&mut self, dir: &Path) -> Result<usize, LoadError> {
        let mut paths: Vec<_> = std::fs::read_dir(dir)
            .map_err(LoadError::Io)?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "ron"))
            .collect();
        paths.sort();
        for p in &paths {
            let src = std::fs::read_to_string(p).map_err(LoadError::Io)?;
            self.add_ron(&p.display().to_string(), &src)?;
        }
        Ok(paths.len())
    }

    pub fn add_ron(&mut self, source: &str, ron_text: &str) -> Result<(), LoadError> {
        let def: PartDef = ron::from_str(ron_text)
            .map_err(|error| LoadError::Parse { source: source.to_string(), error: Box::new(error) })?;
        if self.parts.contains_key(&def.id) {
            return Err(LoadError::DuplicateId(def.id));
        }
        self.parts.insert(def.id.clone(), Arc::new(def));
        Ok(())
    }

    pub fn get(&self, id: &str) -> Option<Arc<PartDef>> {
        self.parts.get(id).cloned()
    }

    /// Wie `get`, bricht aber mit klarer Meldung ab.
    pub fn part(&self, id: &str) -> Arc<PartDef> {
        self.get(id).unwrap_or_else(|| panic!("unbekanntes Teil: {id}"))
    }

    pub fn ids(&self) -> impl Iterator<Item = &str> {
        self.parts.keys().map(String::as_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_library_has_five_parts() {
        let lib = PartLibrary::builtin();
        let ids: Vec<_> = lib.ids().collect();
        assert_eq!(ids, ["entkoppler", "kapsel", "steuerflaeche", "tank", "triebwerk"]);
        assert!(lib.part("triebwerk").engine().is_some());
        assert_eq!(lib.part("tank").fuel_capacity(), 2000.0);
        assert!(lib.part("entkoppler").is_decoupler());
        assert!(lib.part("steuerflaeche").control_surface().is_some());
        assert_eq!(lib.part("steuerflaeche").attach, AttachKind::Radial);
    }

    #[test]
    fn thrust_follows_ambient_pressure() {
        let e = PartLibrary::builtin().part("triebwerk").engine().unwrap();
        assert!((e.thrust(1.0, 0.0) - 215_000.0).abs() < 1e-6);
        let sl = e.thrust(1.0, P_REF);
        assert!((sl - 215_000.0 * 280.0 / 320.0).abs() < 1e-6);
        // stetig dazwischen
        let mid = e.thrust(1.0, P_REF / 2.0);
        assert!(mid > sl && mid < 215_000.0);
    }

    #[test]
    fn cone_com_is_closer_to_wide_end() {
        let s = Shape::Cone { radius_bottom: 1.0, radius_top: 0.0, length: 4.0 };
        // Vollkegel: Schwerpunkt bei h/4 über der Basis → -1.0 relativ zur Mitte
        assert!((s.com_local().z() - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn parse_error_names_the_source() {
        let mut lib = PartLibrary::default();
        let err = lib.add_ron("kaputt.ron", "(id: \"x\")").unwrap_err();
        assert!(err.to_string().starts_with("kaputt.ron"));
    }
}
