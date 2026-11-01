//! Himmelskörper.

use crate::atmosphere::Atmosphere;
use crate::frame::{BodyFixed, Inertial, Rotation, Vec3};
use crate::universe::{BodyId, Orbit};
use serde::{Deserialize, Serialize};
use std::f64::consts::TAU;

/// Darstellungshinweis für einen Himmelskörper (für die Physik ohne Belang).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Appearance {
    /// Leuchtender Stern mit Farbe (linear RGB).
    Star { color: (f32, f32, f32) },
    /// Planet mit Meeren, Land und Eiskappen.
    #[default]
    Terran,
    /// Atmosphärloser, grauer Felskörper.
    Barren,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CelestialBody {
    /// Kennung im Universum.
    pub id: BodyId,
    pub name: String,
    /// Körper, um den dieser kreist (`None` für Sterne).
    pub parent: Option<BodyId>,
    /// Bahn um den Elternkörper.
    pub orbit: Option<Orbit>,
    /// Gravitationsparameter μ = G·M in m³/s².
    pub mu: f64,
    /// Mittlerer Radius in m.
    pub radius: f64,
    /// Radius der Einflusssphäre (Patched Conics) in m. Wird beim Laden des
    /// Universums aus den Massen berechnet; für Sterne unendlich.
    pub soi_radius: f64,
    /// Immer vorhanden; ohne Lufthülle ist es [`Atmosphere::VACUUM`].
    pub atmosphere: Atmosphere,
    /// Siderische Drehrate um die inertiale +Z-Achse in rad/s
    /// (0 = dreht sich nicht; kein Sonderfall im Code).
    pub rotation_rate: f64,
    /// Drehwinkel zur Zeit 0 in rad.
    pub rotation_offset: f64,
    pub appearance: Appearance,
}

/// Geografische Lage über einem Himmelskörper.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeoPosition {
    /// Breite in rad (+ = Norden).
    pub latitude: f64,
    /// Länge in rad (+ = Osten), −π..π.
    pub longitude: f64,
    /// Höhe über dem mittleren Radius in m.
    pub altitude: f64,
}

impl CelestialBody {
    /// Drehrate aus der siderischen Rotationsdauer in s.
    pub fn rate_from_period(period: f64) -> f64 {
        TAU / period
    }

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

    /// Lage des körperfesten Systems zur Zeit `t`.
    pub fn orientation(&self, t: f64) -> Rotation<BodyFixed, Inertial> {
        Rotation::from_quat(glam::DQuat::from_rotation_z(self.rotation_offset + self.rotation_rate * t))
    }

    /// Winkelgeschwindigkeit des Körpers.
    pub fn angular_velocity(&self) -> Vec3<Inertial> {
        Vec3::new(0.0, 0.0, self.rotation_rate)
    }

    /// Geschwindigkeit eines Punktes, der mit dem Körper mitdreht, also
    /// auch die Geschwindigkeit der Luft dort.
    pub fn corotation_velocity(&self, pos: Vec3<Inertial>) -> Vec3<Inertial> {
        self.angular_velocity().cross(pos)
    }

    pub fn to_body_fixed(&self, pos: Vec3<Inertial>, t: f64) -> Vec3<BodyFixed> {
        self.orientation(t).inverse().apply(pos)
    }

    pub fn geo_position(&self, pos: Vec3<Inertial>, t: f64) -> GeoPosition {
        let p = self.to_body_fixed(pos, t);
        let r = p.length();
        GeoPosition {
            latitude: (p.z() / r).clamp(-1.0, 1.0).asin(),
            longitude: p.y().atan2(p.x()),
            altitude: r - self.radius,
        }
    }

    pub fn geo_to_body_fixed(&self, g: GeoPosition) -> Vec3<BodyFixed> {
        let r = self.radius + g.altitude;
        let (sl, cl) = g.latitude.sin_cos();
        let (so, co) = g.longitude.sin_cos();
        Vec3::new(r * cl * co, r * cl * so, r * sl)
    }

    /// Großkreisentfernung auf der Oberfläche zwischen zwei Orten.
    pub fn surface_distance(&self, a: GeoPosition, b: GeoPosition) -> f64 {
        let ua = self.geo_to_body_fixed(GeoPosition { altitude: 0.0, ..a }).normalize_or_zero();
        let ub = self.geo_to_body_fixed(GeoPosition { altitude: 0.0, ..b }).normalize_or_zero();
        self.radius * ua.cross(ub).length().atan2(ua.dot(ub))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn planet() -> CelestialBody {
        CelestialBody {
            id: BodyId(0),
            name: "Test".into(),
            parent: None,
            orbit: None,
            mu: 3.5316e12,
            radius: 600_000.0,
            soi_radius: 8.4e7,
            atmosphere: Atmosphere::VACUUM,
            rotation_rate: CelestialBody::rate_from_period(21_549.425),
            rotation_offset: 0.3,
            appearance: Appearance::Terran,
        }
    }

    #[test]
    fn equator_moves_east_at_expected_speed() {
        let b = planet();
        let p = Vec3::<Inertial>::new(b.radius, 0.0, 0.0);
        let v = b.corotation_velocity(p);
        assert!((v.length() - 174.94).abs() < 0.01, "{}", v.length());
        assert!(v.y() > 0.0); // ostwärts
    }

    #[test]
    fn one_rotation_returns_to_start() {
        let b = planet();
        let period = TAU / b.rotation_rate;
        let p = Vec3::<BodyFixed>::new(1.0, 2.0, 3.0);
        let a = b.orientation(0.0).apply(p);
        let c = b.orientation(period).apply(p);
        assert!((a - c).length() < 1e-9);
    }

    #[test]
    fn geo_roundtrip_and_fixed_point_stays_put() {
        let b = planet();
        let g = GeoPosition { latitude: 0.4, longitude: -2.0, altitude: 1234.0 };
        // Ein mitdrehender Punkt hat zu jeder Zeit dieselben Koordinaten.
        for t in [0.0, 100.0, 5_000.0] {
            let inertial = b.orientation(t).apply(b.geo_to_body_fixed(g));
            let back = b.geo_position(inertial, t);
            assert!((back.latitude - g.latitude).abs() < 1e-12);
            assert!((back.longitude - g.longitude).abs() < 1e-12);
            assert!((back.altitude - g.altitude).abs() < 1e-6);
        }
    }

    #[test]
    fn quarter_circumference() {
        let b = planet();
        let a = GeoPosition { latitude: 0.0, longitude: 0.0, altitude: 0.0 };
        let c = GeoPosition { latitude: 0.0, longitude: std::f64::consts::FRAC_PI_2, altitude: 500.0 };
        assert!((b.surface_distance(a, c) - b.radius * std::f64::consts::FRAC_PI_2).abs() < 1e-6);
    }
}
