//! Simulationskern.
//!
//! Grundregeln:
//! - Alles in `f64` und SI-Einheiten (m, s, kg, N).
//! - Jeder Vektor trägt sein Bezugssystem im Typ ([`frame::Vec3`]).
//!   Vektoren aus verschiedenen Systemen lassen sich nicht mischen;
//!   Umrechnung nur über eine explizite [`frame::Rotation`].
//! - Keine Abhängigkeit von Rendering oder Engine: alles ist headless testbar.

pub mod atmosphere;
pub mod body;
pub mod frame;
pub mod integrate;
pub mod kepler;

pub use atmosphere::{Atmosphere, AtmosphereSample};
pub use body::{CelestialBody, GeoPosition};
pub use frame::{BodyFixed, Frame, Inertial, PartLocal, Rotation, Vec3, VesselFrame};
pub use kepler::{Elements, StateVector};

/// Normfallbeschleunigung, Basis für den spezifischen Impuls.
pub const G0: f64 = 9.806_65;
