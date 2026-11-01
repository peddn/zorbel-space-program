//! Schiffe aus Bauteilen.
//!
//! Physikalisches Modell (siehe Projektnotizen): Ein Schiff ist **ein**
//! starrer Körper. Die Verbindungen zwischen Teilen sind keine Federn und
//! kommen im Bewegungs-Solver nicht vor. Sie existieren nur als Graph,
//! auf dem Staging, Treibstofffluss und die Lastberechnung arbeiten.
//! Überschreitet die Last einer Verbindung ihre Festigkeit, zerfällt das
//! Schiff in zwei neue starre Körper.
//!
//! Luftkräfte wirken an jedem Teil einzeln und sind proportional zur
//! Luftdichte; Vakuum ist der Fall Dichte = 0, kein eigener Modus.

pub mod aero;
pub mod autopilot;
pub mod dynamics;
pub mod geometry;
pub mod ground;
pub mod loads;
pub mod mass;
pub mod part;
pub mod scenario;
pub mod vessel;

pub use dynamics::{Controls, RigidState, Separated, StepReport};
pub use loads::{ExternalLoad, JointLoad, LoadAnalysis, LoadKind, Motion};
pub use mass::{InertiaTensor, MassProperties};
pub use part::{AttachKind, JointStrength, Module, PartDef, PartLibrary, Shape};
pub use vessel::{Joint, PartInstance, Vessel, VesselBuilder};
