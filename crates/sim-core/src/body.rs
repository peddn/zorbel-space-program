//! Himmelskörper.

use crate::atmosphere::Atmosphere;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CelestialBody {
    pub name: String,
    /// Gravitationsparameter μ = G·M in m³/s².
    pub mu: f64,
    /// Mittlerer Radius in m.
    pub radius: f64,
    /// Radius der Einflusssphäre (Patched Conics) in m.
    pub soi_radius: f64,
    /// Immer vorhanden; ohne Lufthülle ist es [`Atmosphere::VACUUM`].
    pub atmosphere: Atmosphere,
}

impl CelestialBody {
    /// Schwerebeschleunigung an der Oberfläche.
    pub fn surface_gravity(&self) -> f64 {
        self.mu / (self.radius * self.radius)
    }

    /// Kreisbahngeschwindigkeit in Höhe `altitude` über der Oberfläche.
    pub fn circular_speed(&self, altitude: f64) -> f64 {
        (self.mu / (self.radius + altitude)).sqrt()
    }

    /// Radius, oberhalb dessen außer der Gravitation keine Umgebungskräfte
    /// wirken: Oberrand der Atmosphäre, bei Vakuum die Oberfläche.
    pub fn quiet_radius(&self) -> f64 {
        self.radius + self.atmosphere.top()
    }
}
