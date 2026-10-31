//! Aggregation von Masse, Schwerpunkt und Trägheit über alle Teile.

use crate::vessel::Vessel;
use glam::{DMat3, DVec3};
use sim_core::{Vec3, VesselFrame};

/// Trägheitstensor, Achsen des Schiffssystems.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InertiaTensor(pub DMat3);

impl InertiaTensor {
    pub fn apply(&self, v: Vec3<VesselFrame>) -> Vec3<VesselFrame> {
        Vec3::from_raw(self.0 * v.raw())
    }
    pub fn solve(&self, torque: Vec3<VesselFrame>) -> Vec3<VesselFrame> {
        Vec3::from_raw(self.0.inverse() * torque.raw())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct MassProperties {
    pub mass: f64,
    pub com: Vec3<VesselFrame>,
    /// Um den Gesamtschwerpunkt.
    pub inertia: InertiaTensor,
}

/// Steiner-Anteil m·(|d|²·E − d·dᵀ).
fn parallel_axis(m: f64, d: DVec3) -> DMat3 {
    let outer = DMat3::from_cols(d * d.x, d * d.y, d * d.z);
    (DMat3::IDENTITY * d.length_squared() - outer) * m
}

pub fn mass_properties(v: &Vessel) -> MassProperties {
    let mass: f64 = v.parts.iter().map(|p| p.mass()).sum();
    let com = v.parts.iter().map(|p| p.com() * p.mass()).sum::<Vec3<VesselFrame>>() / mass;
    let inertia = v
        .parts
        .iter()
        .map(|p| p.inertia() + parallel_axis(p.mass(), (p.com() - com).raw()))
        .fold(DMat3::ZERO, |a, b| a + b);
    MassProperties { mass, com, inertia: InertiaTensor(inertia) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::PartLibrary;
    use crate::vessel::VesselBuilder;

    #[test]
    fn two_stacked_tanks_equal_one_long_cylinder() {
        let lib = PartLibrary::builtin();
        let tank = lib.part("tank");
        let mut b = VesselBuilder::new("Test", tank.clone());
        b.stack_below(0, tank.clone());
        let v = b.build();
        let mp = mass_properties(&v);

        let m = 2.0 * (tank.dry_mass + tank.fuel_capacity());
        let (r, l) = (0.625, 1.9);
        assert!((mp.mass - m).abs() < 1e-9);
        // Schwerpunkt genau an der Fuge
        assert!((mp.com - Vec3::new(0.0, 0.0, -l / 2.0)).length() < 1e-12);
        // Zylinder der Länge 2l
        let ixx = m * (3.0 * r * r + (2.0 * l) * (2.0 * l)) / 12.0;
        let izz = 0.5 * m * r * r;
        let i = mp.inertia.0;
        assert!((i.x_axis.x - ixx).abs() / ixx < 1e-12);
        assert!((i.y_axis.y - ixx).abs() / ixx < 1e-12);
        assert!((i.z_axis.z - izz).abs() / izz < 1e-12);
    }

    #[test]
    fn symmetric_fins_keep_com_on_axis() {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("tank"));
        b.radial_symmetric(0, lib.part("steuerflaeche"), 4, -0.5);
        let mp = mass_properties(&b.build());
        assert!(mp.com.x().abs() < 1e-12 && mp.com.y().abs() < 1e-12);
        // Symmetrie: keine Deviationsmomente
        let i = mp.inertia.0;
        assert!(i.x_axis.y.abs() < 1e-9 && i.x_axis.z.abs() < 1e-9 && i.y_axis.z.abs() < 1e-9);
    }
}
