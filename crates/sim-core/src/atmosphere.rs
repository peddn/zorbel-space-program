//! Atmosphärenmodell.
//!
//! Jeder Himmelskörper hat eine Atmosphäre. Ein Körper ohne Lufthülle ist
//! kein Sonderfall, sondern eine Atmosphäre mit Bodendruck 0
//! ([`Atmosphere::VACUUM`]). Alle Größen laufen stetig und glatt gegen
//! null; am Oberrand verschwindet der Druck mit stetiger erster und
//! zweiter Ableitung. Beim Ein- oder Austritt gibt es daher keinen Sprung
//! in irgendeiner Kraft.

use serde::{Deserialize, Serialize};

/// Universelle Gaskonstante in J/(mol·K).
pub const R_GAS: f64 = 8.314_462_618;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Atmosphere {
    /// Bodendruck in Pa. 0 = Vakuum.
    pub surface_pressure: f64,
    /// Skalenhöhe des exponentiellen Druckabfalls in m.
    pub scale_height: f64,
    /// Höhe, in der der Druck genau 0 erreicht (glatt ausgeblendet).
    pub top: f64,
    /// Bodentemperatur in K.
    pub surface_temperature: f64,
    /// Temperaturgradient in K/m.
    pub lapse_rate: f64,
    /// Untergrenze der Temperatur in K.
    pub min_temperature: f64,
    /// Molare Masse des Gases in kg/mol.
    pub molar_mass: f64,
    /// Adiabatenexponent.
    pub gamma: f64,
}

/// Zustand der Luft an einem Ort.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AtmosphereSample {
    pub pressure: f64,
    pub density: f64,
    pub temperature: f64,
    pub speed_of_sound: f64,
}

/// Glatte Stufe 0→1 auf [0,1] mit verschwindender 1. und 2. Ableitung an den Enden.
fn smootherstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

/// Anteil der Höhe, über den der Druck am Oberrand ausgeblendet wird.
const FADE_FRACTION: f64 = 0.25;

impl Atmosphere {
    /// Lufthülle mit Bodendruck 0: dasselbe Modell, nur ohne Luft.
    pub const VACUUM: Self = Self {
        surface_pressure: 0.0,
        scale_height: 1.0,
        top: 0.0,
        surface_temperature: 3.0,
        lapse_rate: 0.0,
        min_temperature: 3.0,
        molar_mass: 0.028_964_4,
        gamma: 1.4,
    };

    /// Luftähnliche Atmosphäre mit KSP-artigen Kenndaten (Oberrand 70 km).
    pub fn earthlike() -> Self {
        Self {
            surface_pressure: 101_325.0,
            scale_height: 5_600.0,
            top: 70_000.0,
            surface_temperature: 288.15,
            lapse_rate: 0.0065,
            min_temperature: 200.0,
            molar_mass: 0.028_964_4,
            gamma: 1.4,
        }
    }

    /// Höhe, ab der keine Luftkräfte mehr auftreten.
    pub fn top(&self) -> f64 {
        if self.surface_pressure > 0.0 { self.top } else { 0.0 }
    }

    pub fn sample(&self, altitude: f64) -> AtmosphereSample {
        let temperature = (self.surface_temperature - self.lapse_rate * altitude.max(0.0)).max(self.min_temperature);
        let speed_of_sound = (self.gamma * R_GAS * temperature / self.molar_mass).sqrt();

        let fade = if self.top > 0.0 {
            let start = self.top * (1.0 - FADE_FRACTION);
            1.0 - smootherstep((altitude - start) / (self.top - start))
        } else {
            0.0
        };
        let pressure = self.surface_pressure * (-altitude / self.scale_height).exp() * fade;
        let density = pressure * self.molar_mass / (R_GAS * temperature);
        AtmosphereSample { pressure, density, temperature, speed_of_sound }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sea_level_density_is_plausible() {
        let s = Atmosphere::earthlike().sample(0.0);
        assert!((s.density - 1.225).abs() < 0.01, "{}", s.density);
        assert!((s.speed_of_sound - 340.3).abs() < 0.5);
    }

    #[test]
    fn vacuum_is_the_zero_pressure_case() {
        let s = Atmosphere::VACUUM.sample(1234.0);
        assert_eq!(s.pressure, 0.0);
        assert_eq!(s.density, 0.0);
        assert!(s.speed_of_sound > 0.0); // Machzahl bleibt definiert
    }

    #[test]
    fn pressure_vanishes_smoothly_at_the_top() {
        let a = Atmosphere::earthlike();
        assert_eq!(a.sample(a.top).pressure, 0.0);
        assert_eq!(a.sample(a.top + 1.0).pressure, 0.0);
        // Stetig: knapp unterhalb fast null, Ableitung ebenfalls
        let h = 1.0;
        let p1 = a.sample(a.top - h).pressure;
        let p2 = a.sample(a.top - 2.0 * h).pressure;
        assert!(p1 < 1e-12 * a.surface_pressure);
        assert!((p2 - p1).abs() / h < 1e-11 * a.surface_pressure);
        // Monoton fallend über die ganze Höhe
        let mut last = f64::INFINITY;
        for k in 0..=700 {
            let p = a.sample(k as f64 * 100.0).pressure;
            assert!(p <= last);
            last = p;
        }
    }
}
