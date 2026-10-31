//! Lasten in den Verbindungen, berechnet aus der Starrkörperbewegung.
//!
//! Idee: Jedes Teil i muss die Beschleunigung a_i des starren Körpers an
//! seinem Ort mitmachen. Dafür braucht es die Nettokraft m_i·a_i. Was davon
//! nicht durch äußere Kräfte am Teil selbst (Schub, Auftrieb, ...) gedeckt
//! ist, muss über die Verbindungen kommen. Summiert man das über den
//! Teilbaum unterhalb einer Verbindung, erhält man genau die Kraft und das
//! Moment, das der Elternteil auf diesen Teilbaum ausübt.
//!
//! Gravitation taucht nicht auf: In einem homogenen Schwerefeld wird jedes
//! Teil gleich beschleunigt, das erzeugt keine inneren Lasten. Deshalb wird
//! mit der *Eigenbeschleunigung* (ohne Gravitation) gerechnet, also dem,
//! was ein Beschleunigungsmesser an Bord anzeigen würde.
//!
//! Aufwand: ein Durchlauf über alle Teile, O(n).

use crate::mass::{InertiaTensor, MassProperties};
use crate::vessel::Vessel;
use sim_core::{Vec3, VesselFrame};

/// Äußere Kraft (ohne Gravitation) an einem Teil.
#[derive(Debug, Clone, Copy)]
pub struct ExternalLoad {
    pub part: usize,
    pub force: Vec3<VesselFrame>,
    /// Angriffspunkt der Kraft.
    pub point: Vec3<VesselFrame>,
    /// Reines Drehmoment (z. B. Reaktionsrad).
    pub torque: Vec3<VesselFrame>,
}

/// Bewegungszustand des starren Körpers, alles im Schiffssystem.
#[derive(Debug, Clone, Copy)]
pub struct Motion {
    /// Eigenbeschleunigung des Schwerpunkts (ohne Gravitation).
    pub accel: Vec3<VesselFrame>,
    pub omega: Vec3<VesselFrame>,
    pub alpha: Vec3<VesselFrame>,
}

/// Bewegung, die aus den äußeren Lasten folgt (Newton + Euler).
pub fn motion_from_loads(mp: &MassProperties, ext: &[ExternalLoad], omega: Vec3<VesselFrame>) -> Motion {
    let force: Vec3<VesselFrame> = ext.iter().map(|e| e.force).sum();
    let torque: Vec3<VesselFrame> = ext.iter().map(|e| (e.point - mp.com).cross(e.force) + e.torque).sum();
    let gyro = omega.cross(mp.inertia.apply(omega));
    Motion { accel: force / mp.mass, omega, alpha: mp.inertia.solve(torque - gyro) }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoadKind {
    Tension,
    Compression,
    Shear,
    Bending,
    Torsion,
}

impl LoadKind {
    pub fn label(self) -> &'static str {
        match self {
            LoadKind::Tension => "Zug",
            LoadKind::Compression => "Druck",
            LoadKind::Shear => "Scherung",
            LoadKind::Bending => "Biegung",
            LoadKind::Torsion => "Torsion",
        }
    }
}

/// Last in einer Verbindung (vom Eltern- auf den Kind-Teilbaum).
#[derive(Debug, Clone, Copy)]
pub struct JointLoad {
    pub child: usize,
    pub parent: usize,
    pub force: Vec3<VesselFrame>,
    /// Moment um den Verbindungspunkt.
    pub moment: Vec3<VesselFrame>,
    /// Axialkraft: positiv = Zug, negativ = Druck.
    pub axial: f64,
    pub shear: f64,
    pub bending: f64,
    pub torsion: f64,
    /// Größtes Verhältnis Last/Festigkeit; > 1 bedeutet Bruch.
    pub utilization: f64,
    pub critical: LoadKind,
}

#[derive(Debug, Clone)]
pub struct LoadAnalysis {
    pub joints: Vec<JointLoad>,
    /// Was am Wurzelteil übrig bleibt. Physikalisch muss das null sein;
    /// dient als Selbsttest der Rechnung.
    pub residual_force: Vec3<VesselFrame>,
    pub residual_moment: Vec3<VesselFrame>,
}

impl LoadAnalysis {
    pub fn max_utilization(&self) -> Option<&JointLoad> {
        self.joints.iter().max_by(|a, b| a.utilization.total_cmp(&b.utilization))
    }
    /// Verbindungen über ihrer Festigkeit (Index des Kindteils).
    pub fn broken(&self) -> impl Iterator<Item = &JointLoad> {
        self.joints.iter().filter(|j| j.utilization > 1.0)
    }
}

pub fn joint_loads(v: &Vessel, mp: &MassProperties, m: &Motion, ext: &[ExternalLoad]) -> LoadAnalysis {
    let n = v.parts.len();
    // Kraft und Moment (um den Ursprung des Schiffssystems), die jedes
    // Teil von außerhalb seines Teilbaums braucht.
    let mut f = vec![Vec3::<VesselFrame>::ZERO; n];
    let mut m0 = vec![Vec3::<VesselFrame>::ZERO; n];

    for (i, p) in v.parts.iter().enumerate() {
        let mi = p.mass();
        let c = p.com();
        let d = c - mp.com;
        let a_i = m.accel + m.alpha.cross(d) + m.omega.cross(m.omega.cross(d));
        let ii = InertiaTensor(p.inertia());
        let needed = a_i * mi;
        f[i] = needed;
        m0[i] = c.cross(needed) + ii.apply(m.alpha) + m.omega.cross(ii.apply(m.omega));
    }
    for e in ext {
        f[e.part] -= e.force;
        m0[e.part] -= e.point.cross(e.force) + e.torque;
    }

    // Von hinten nach vorn: Kinder vor Eltern (Index-Invariante).
    let mut joints = Vec::with_capacity(n.saturating_sub(1));
    for i in (1..n).rev() {
        let p = &v.parts[i];
        let parent = p.parent.expect("nur die Wurzel hat kein Elternteil");
        debug_assert!(parent < i, "Index-Invariante verletzt");
        let joint = p.joint.expect("Teil mit Elternteil braucht Verbindung");

        let force = f[i];
        let moment = m0[i] - joint.point.cross(force);
        let ax = joint.axis;
        // Kraft auf das Kind in Richtung Kind (+axis) = Druck, Richtung Eltern = Zug.
        let axial = -force.dot(ax);
        let shear = force.reject_from_unit(ax).length();
        let torsion = moment.dot(ax).abs();
        let bending = moment.reject_from_unit(ax).length();

        let s = joint.strength;
        let candidates = [
            if axial >= 0.0 { (axial / s.tension, LoadKind::Tension) } else { (-axial / s.compression, LoadKind::Compression) },
            (shear / s.shear, LoadKind::Shear),
            (bending / s.bending, LoadKind::Bending),
            (torsion / s.torsion, LoadKind::Torsion),
        ];
        let (utilization, critical) =
            candidates.into_iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();

        joints.push(JointLoad { child: i, parent, force, moment, axial, shear, bending, torsion, utilization, critical });

        let m_i = m0[i];
        f[parent] += force;
        m0[parent] += m_i;
    }
    joints.reverse();

    LoadAnalysis { joints, residual_force: f[0], residual_moment: m0[0] - mp.com.cross(f[0]) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mass::mass_properties;
    use crate::part::PartLibrary;
    use crate::vessel::VesselBuilder;

    /// Kapsel / Tank / Triebwerk, Triebwerk schiebt mit Schub T.
    fn simple_rocket(thrust: f64) -> (Vessel, Vec<ExternalLoad>) {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("kapsel"));
        let t = b.stack_below(0, lib.part("tank"));
        let e = b.stack_below(t, lib.part("triebwerk"));
        let v = b.build();
        let ext = vec![ExternalLoad {
            part: e,
            force: Vec3::new(0.0, 0.0, thrust),
            point: v.parts[e].origin,
            torque: Vec3::ZERO,
        }];
        (v, ext)
    }

    #[test]
    fn engine_joint_compression_matches_hand_calculation() {
        let thrust = 215_000.0;
        let (v, ext) = simple_rocket(thrust);
        let mp = mass_properties(&v);
        let m = motion_from_loads(&mp, &ext, Vec3::ZERO);
        let la = joint_loads(&v, &mp, &m, &ext);

        let a = thrust / mp.mass;
        let engine = la.joints.iter().find(|j| j.child == 2).unwrap();
        // Das Triebwerk schiebt den Rest (Kapsel + Tank) mit a an.
        let expected = (v.parts[0].mass() + v.parts[1].mass()) * a;
        assert_eq!(engine.critical, LoadKind::Compression);
        assert!((-engine.axial - expected).abs() < 1e-6);
        assert!(engine.shear < 1e-9 && engine.bending < 1e-6);

        // Verbindung Kapsel/Tank trägt nur die Kapsel
        let upper = la.joints.iter().find(|j| j.child == 1).unwrap();
        assert!((-upper.axial - v.parts[0].mass() * a).abs() < 1e-6);
    }

    #[test]
    fn residual_vanishes_for_arbitrary_loads() {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("kapsel"));
        let t = b.stack_below(0, lib.part("tank"));
        let fins = b.radial_symmetric(t, lib.part("steuerflaeche"), 3, -0.3);
        let e = b.stack_below(t, lib.part("triebwerk"));
        let v = b.build();
        let mp = mass_properties(&v);
        let ext = vec![
            ExternalLoad { part: e, force: Vec3::new(3_000.0, -1_000.0, 200_000.0), point: v.parts[e].origin, torque: Vec3::ZERO },
            ExternalLoad { part: fins[1], force: Vec3::new(500.0, 2_000.0, -300.0), point: v.parts[fins[1]].origin, torque: Vec3::ZERO },
            ExternalLoad { part: 0, force: Vec3::ZERO, point: Vec3::ZERO, torque: Vec3::new(1_000.0, 0.0, 4_000.0) },
        ];
        let m = motion_from_loads(&mp, &ext, Vec3::new(0.3, -0.2, 1.1));
        let la = joint_loads(&v, &mp, &m, &ext);
        assert!(la.residual_force.length() < 1e-6, "{:?}", la.residual_force);
        assert!(la.residual_moment.length() < 1e-5, "{:?}", la.residual_moment);
    }

    #[test]
    fn spinning_fins_feel_centripetal_tension() {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("tank"));
        b.radial_symmetric(0, lib.part("steuerflaeche"), 4, 0.0);
        let v = b.build();
        let mp = mass_properties(&v);
        let omega = 3.0; // rad/s um die Längsachse
        let m = motion_from_loads(&mp, &[], Vec3::new(0.0, 0.0, omega));
        let la = joint_loads(&v, &mp, &m, &[]);

        let fin = &v.parts[1];
        let r = fin.com().reject_from_unit(Vec3::Z).length();
        let expected = fin.mass() * omega * omega * r;
        for j in &la.joints {
            assert_eq!(j.critical, LoadKind::Tension);
            assert!((j.axial - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn excessive_thrust_breaks_the_engine_joint() {
        let (v, ext) = simple_rocket(5.0e6);
        let mp = mass_properties(&v);
        let m = motion_from_loads(&mp, &ext, Vec3::ZERO);
        let la = joint_loads(&v, &mp, &m, &ext);
        let broken: Vec<_> = la.broken().map(|j| j.child).collect();
        assert!(broken.contains(&2));
    }
}
