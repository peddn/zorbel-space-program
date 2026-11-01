//! Das Universum: Sternsysteme und ihre Himmelskörper.
//!
//! Aufbau: Galaxie → Sternsysteme → Stern → Planeten → Monde. Jeder Körper
//! außer den Sternen fliegt auf einer Keplerbahn um seinen Elternkörper.
//! Alle Systeme teilen dieselbe Achsenausrichtung (galaktische Achsen).
//!
//! **Genauigkeit:** Innerhalb eines Systems wird in `f64` relativ zum Stern
//! bzw. zum aktuellen Körper gerechnet; das reicht bis weit über die äußeren
//! Planeten auf Millimeter. Sterne haben ganzzahlige Galaxie-Positionen in
//! Mikrometern ([`GalacticPos`]); Abstände zwischen Systemen sind damit
//! exakt. Zwischen Systemen wird nur gesprungen, es gibt keinen Flug durch
//! den interstellaren Raum.
//!
//! Die Daten kommen aus einer RON-Datei (`assets/universe.ron`).

use crate::atmosphere::Atmosphere;
use crate::body::{Appearance, CelestialBody};
use crate::frame::{Inertial, Vec3};
use crate::kepler::StateVector;
use glam::{DQuat, DVec3};
use serde::{Deserialize, Serialize};
use std::f64::consts::TAU;
use std::fmt;

/// Lichtjahr in Metern.
pub const LIGHT_YEAR: f64 = 9.460_730_472_580_8e15;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct BodyId(pub usize);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct SystemId(pub usize);

/// Klassische Bahnelemente relativ zum Elternkörper (Winkel in rad).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Orbit {
    pub semi_major_axis: f64,
    #[serde(default)]
    pub eccentricity: f64,
    #[serde(default)]
    pub inclination: f64,
    /// Länge des aufsteigenden Knotens.
    #[serde(default)]
    pub ascending_node: f64,
    #[serde(default)]
    pub argument_of_periapsis: f64,
    #[serde(default)]
    pub mean_anomaly_at_epoch: f64,
    #[serde(default)]
    pub epoch: f64,
}

impl Orbit {
    pub fn period(&self, mu_parent: f64) -> f64 {
        TAU * (self.semi_major_axis.powi(3) / mu_parent).sqrt()
    }

    /// Zustand zur Zeit `t` relativ zum Elternkörper (nur geschlossene Bahnen).
    pub fn state_at(&self, mu_parent: f64, t: f64) -> StateVector<Inertial> {
        let (a, e) = (self.semi_major_axis, self.eccentricity);
        let n = (mu_parent / a.powi(3)).sqrt();
        let m = (self.mean_anomaly_at_epoch + n * (t - self.epoch)).rem_euclid(TAU);
        // Kepler-Gleichung E − e·sin E = M (Newton)
        let mut ecc = if e < 0.8 { m } else { std::f64::consts::PI };
        for _ in 0..50 {
            let step = (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
            ecc -= step;
            if step.abs() < 1e-14 {
                break;
            }
        }
        let (se, ce) = ecc.sin_cos();
        let b = a * (1.0 - e * e).sqrt();
        let r = a * (1.0 - e * ce);
        let pos = DVec3::new(a * (ce - e), b * se, 0.0);
        let vel = DVec3::new(-a * se, b * ce, 0.0) * ((mu_parent * a).sqrt() / r / a);
        let rot = DQuat::from_rotation_z(self.ascending_node)
            * DQuat::from_rotation_x(self.inclination)
            * DQuat::from_rotation_z(self.argument_of_periapsis);
        StateVector::new(Vec3::from_raw(rot * pos), Vec3::from_raw(rot * vel))
    }
}

/// Ganzzahlige Position in der Galaxie (Mikrometer). Reicht für das ganze
/// beobachtbare Universum, Addieren und Subtrahieren sind exakt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub struct GalacticPos {
    pub x: i128,
    pub y: i128,
    pub z: i128,
}

impl GalacticPos {
    pub const MICROMETER: f64 = 1e-6;

    pub fn from_light_years(ly: (f64, f64, f64)) -> Self {
        let um = |v: f64| (v * LIGHT_YEAR / Self::MICROMETER).round() as i128;
        Self { x: um(ly.0), y: um(ly.1), z: um(ly.2) }
    }

    /// Abstandsvektor `other − self` in Metern.
    pub fn offset_to(&self, other: &GalacticPos) -> DVec3 {
        let d = |a: i128, b: i128| (b - a) as f64 * Self::MICROMETER;
        DVec3::new(d(self.x, other.x), d(self.y, other.y), d(self.z, other.z))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StarSystem {
    pub id: SystemId,
    pub name: String,
    pub star: BodyId,
    pub position: GalacticPos,
    /// Eigenbewegung des Sterns in der Galaxie (m/s).
    pub velocity: DVec3,
    /// Ring der Sprungorbits um den Stern: (innerer, äußerer) Radius in m.
    pub jump_zone: (f64, f64),
}

#[derive(Debug, Clone)]
pub struct Universe {
    pub systems: Vec<StarSystem>,
    pub bodies: Vec<CelestialBody>,
    system_of: Vec<SystemId>,
}

#[derive(Debug)]
pub enum UniverseError {
    Parse(Box<ron::error::SpannedError>),
    Invalid(String),
}

impl fmt::Display for UniverseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UniverseError::Parse(e) => write!(f, "Universumsdatei: {e}"),
            UniverseError::Invalid(m) => write!(f, "Universumsdatei: {m}"),
        }
    }
}
impl std::error::Error for UniverseError {}

// ---------- Dateiformat ----------

fn vacuum() -> Atmosphere {
    Atmosphere::VACUUM
}

#[derive(Deserialize)]
struct BodyDef {
    name: String,
    mu: f64,
    radius: f64,
    /// Siderische Drehdauer in s (0 = dreht sich nicht).
    #[serde(default)]
    rotation_period: f64,
    #[serde(default)]
    rotation_offset: f64,
    #[serde(default = "vacuum")]
    atmosphere: Atmosphere,
    #[serde(default)]
    orbit: Option<Orbit>,
    #[serde(default)]
    appearance: Appearance,
    #[serde(default)]
    children: Vec<BodyDef>,
}

#[derive(Deserialize)]
struct SystemDef {
    name: String,
    /// Lage in der Galaxie in Lichtjahren.
    position_ly: (f64, f64, f64),
    /// Eigenbewegung in m/s.
    #[serde(default)]
    velocity: (f64, f64, f64),
    jump_zone: (f64, f64),
    star: BodyDef,
}

#[derive(Deserialize)]
struct UniverseDef {
    systems: Vec<SystemDef>,
}

impl Universe {
    /// Das mitgelieferte Universum (Zol und Nachbarsystem).
    pub fn builtin() -> Self {
        Self::from_ron(include_str!("../../../assets/universe.ron")).expect("eingebettetes Universum fehlerhaft")
    }

    pub fn from_ron(text: &str) -> Result<Self, UniverseError> {
        let def: UniverseDef = ron::from_str(text).map_err(|e| UniverseError::Parse(Box::new(e)))?;
        let mut u = Universe { systems: Vec::new(), bodies: Vec::new(), system_of: Vec::new() };
        for (si, sys) in def.systems.into_iter().enumerate() {
            let sid = SystemId(si);
            let star = u.add_body(sid, None, sys.star)?;
            u.systems.push(StarSystem {
                id: sid,
                name: sys.name,
                star,
                position: GalacticPos::from_light_years(sys.position_ly),
                velocity: DVec3::new(sys.velocity.0, sys.velocity.1, sys.velocity.2),
                jump_zone: sys.jump_zone,
            });
        }
        Ok(u)
    }

    fn add_body(&mut self, system: SystemId, parent: Option<BodyId>, def: BodyDef) -> Result<BodyId, UniverseError> {
        let id = BodyId(self.bodies.len());
        if parent.is_some() != def.orbit.is_some() {
            return Err(UniverseError::Invalid(format!("{}: Bahn genau dann, wenn es einen Elternkörper gibt", def.name)));
        }
        // Einflusssphäre nach Laplace: a · (m / M)^(2/5)
        let soi_radius = match (parent, def.orbit) {
            (Some(p), Some(o)) => o.semi_major_axis * (def.mu / self.bodies[p.0].mu).powf(0.4),
            _ => f64::INFINITY,
        };
        self.bodies.push(CelestialBody {
            id,
            name: def.name,
            parent,
            orbit: def.orbit,
            mu: def.mu,
            radius: def.radius,
            soi_radius,
            atmosphere: def.atmosphere,
            rotation_rate: if def.rotation_period > 0.0 { CelestialBody::rate_from_period(def.rotation_period) } else { 0.0 },
            rotation_offset: def.rotation_offset,
            appearance: def.appearance,
        });
        self.system_of.push(system);
        for child in def.children {
            self.add_body(system, Some(id), child)?;
        }
        Ok(id)
    }

    pub fn body(&self, id: BodyId) -> &CelestialBody {
        &self.bodies[id.0]
    }

    pub fn system(&self, id: SystemId) -> &StarSystem {
        &self.systems[id.0]
    }

    pub fn system_of(&self, id: BodyId) -> SystemId {
        self.system_of[id.0]
    }

    pub fn find(&self, name: &str) -> Option<BodyId> {
        self.bodies.iter().find(|b| b.name == name).map(|b| b.id)
    }

    pub fn find_system(&self, name: &str) -> Option<SystemId> {
        self.systems.iter().find(|s| s.name == name).map(|s| s.id)
    }

    pub fn children(&self, id: BodyId) -> impl Iterator<Item = &CelestialBody> + '_ {
        self.bodies.iter().filter(move |b| b.parent == Some(id))
    }

    /// Alle Körper eines Systems.
    pub fn bodies_in(&self, system: SystemId) -> impl Iterator<Item = &CelestialBody> + '_ {
        self.bodies.iter().filter(move |b| self.system_of[b.id.0] == system)
    }

    /// Zustand des Körpers relativ zu seinem Elternkörper (Stern: null).
    pub fn state_rel_parent(&self, id: BodyId, t: f64) -> StateVector<Inertial> {
        let b = self.body(id);
        match (b.parent, b.orbit) {
            (Some(p), Some(o)) => o.state_at(self.body(p).mu, t),
            _ => StateVector::new(Vec3::ZERO, Vec3::ZERO),
        }
    }

    /// Zustand des Körpers relativ zum Stern seines Systems.
    pub fn state_in_system(&self, id: BodyId, t: f64) -> StateVector<Inertial> {
        let mut s = StateVector::new(Vec3::ZERO, Vec3::ZERO);
        let mut cur = Some(id);
        while let Some(c) = cur {
            let rel = self.state_rel_parent(c, t);
            s.r += rel.r;
            s.v += rel.v;
            cur = self.body(c).parent;
        }
        s
    }

    /// Zustand von `to` relativ zu `from` (beide im selben System).
    pub fn relative_state(&self, from: BodyId, to: BodyId, t: f64) -> StateVector<Inertial> {
        assert_eq!(self.system_of(from), self.system_of(to), "Körper in verschiedenen Systemen");
        let a = self.state_in_system(from, t);
        let b = self.state_in_system(to, t);
        StateVector::new(b.r - a.r, b.v - a.v)
    }

    /// Abstand zweier Systeme in Metern.
    pub fn system_distance(&self, a: SystemId, b: SystemId) -> f64 {
        self.system(a).position.offset_to(&self.system(b).position).length()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kepler;

    fn zol() -> Universe {
        Universe::builtin()
    }

    #[test]
    fn names_and_structure() {
        let u = zol();
        let zol = u.find_system("Zol").unwrap();
        let zun = u.find("Zun").unwrap();
        let zorb = u.find("Zorb").unwrap();
        let muz = u.find("Muz").unwrap();
        assert_eq!(u.system(zol).star, zun);
        assert_eq!(u.body(zorb).parent, Some(zun));
        assert_eq!(u.body(muz).parent, Some(zorb));
        assert_eq!(u.systems.len(), 2);
        assert!(u.bodies_in(zol).count() >= 3);
    }

    #[test]
    fn spheres_of_influence_match_ksp_scale() {
        let u = zol();
        let zorb = u.body(u.find("Zorb").unwrap());
        let muz = u.body(u.find("Muz").unwrap());
        assert!((zorb.soi_radius - 84_159_286.0).abs() / 84_159_286.0 < 0.005, "{}", zorb.soi_radius);
        assert!((muz.soi_radius - 2_429_559.0).abs() / 2_429_559.0 < 0.005, "{}", muz.soi_radius);
        assert!(u.body(u.find("Zun").unwrap()).soi_radius.is_infinite());
    }

    #[test]
    fn second_system_is_3_7_light_years_away() {
        let u = zol();
        let d = u.system_distance(SystemId(0), SystemId(1)) / LIGHT_YEAR;
        assert!((d - 3.7).abs() < 1e-5, "{d} Lj");
    }

    #[test]
    fn orbit_state_follows_kepler() {
        let o = Orbit {
            semi_major_axis: 1.2e7,
            eccentricity: 0.3,
            inclination: 0.4,
            ascending_node: 1.1,
            argument_of_periapsis: 0.7,
            mean_anomaly_at_epoch: 2.0,
            epoch: 0.0,
        };
        let mu = 3.5316e12;
        let s0 = o.state_at(mu, 0.0);
        for t in [100.0, 5_000.0, 77_777.0] {
            let a = o.state_at(mu, t);
            let b = kepler::propagate(s0, mu, t);
            assert!((a.r - b.r).length() < 1e-3, "t={t}: {:?}", a.r - b.r);
            assert!((a.v - b.v).length() < 1e-6);
        }
        let el = kepler::elements(s0, mu);
        assert!((el.eccentricity - 0.3).abs() < 1e-9 && (el.inclination - 0.4).abs() < 1e-9);
    }

    #[test]
    fn moon_is_tidally_locked() {
        // Muz zeigt Zorb immer dieselbe Seite
        let u = zol();
        let (zorb, muz) = (u.find("Zorb").unwrap(), u.find("Muz").unwrap());
        let face = |t: f64| {
            let to_planet = u.relative_state(muz, zorb, t).r;
            u.body(muz).to_body_fixed(to_planet, t).normalize_or_zero()
        };
        let f0 = face(0.0);
        for t in [10_000.0, 50_000.0, 123_456.0] {
            assert!((face(t) - f0).length() < 1e-6, "t={t}");
        }
    }
}
