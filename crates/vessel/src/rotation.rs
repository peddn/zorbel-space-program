//! Kräftefreie Drehung eines starren Körpers mit beliebigem Trägheitstensor.
//!
//! Zusammengebaute Schiffe sind unsymmetrisch: Ihr Trägheitstensor ist im
//! Schiffssystem voll besetzt. Er wird auf Hauptachsen gebracht; dort ist
//! die Drehung eine Summe von drei Teilflüssen "Drehung um Hauptachse k mit
//! ω_k = L_k / I_k", die jeder für sich exakt lösbar sind. Sie werden
//! symmetrisch hintereinander ausgeführt (Splitting nach
//! Dullweber/Reich/Webster, 2. Ordnung):
//!
//! ```text
//! ½·X  ½·Y  Z  ½·Y  ½·X
//! ```
//!
//! Eigenschaften, unabhängig von der Schrittweite:
//! - Der Drehimpuls im Inertialsystem bleibt exakt erhalten (jeder
//!   Teilfluss dreht Lage und körperfesten Drehimpuls gegenläufig).
//! - Die Energie hat einen kleinen, beschränkten Fehler ohne Drift
//!   (symplektisches Verfahren).
//! - Ein langer Schritt wird in Teilschritte von höchstens
//!   [`MAX_SUBSTEP_ANGLE`] zerlegt. Darum trudelt ein Schiff im Zeitraffer
//!   genau so wie in Echtzeit, inklusive Dzhanibekov-Effekt.
//!
//! Optional bremsen die Reaktionsräder in jedem Teilschritt (Lagehalter im
//! Zeitraffer): Das Moment je Schiffsachse ist auf das Radmoment begrenzt
//! und wirkt als Drehimpulsstoß.

use crate::mass::InertiaTensor;
use glam::{DMat3, DQuat, DVec3};
use sim_core::{Inertial, Rotation, Vec3, VesselFrame};

/// Größter Drehwinkel pro Teilschritt (rad).
pub const MAX_SUBSTEP_ANGLE: f64 = 0.05;
/// Obergrenze der Teilschritte pro Aufruf (Schutz vor Endlosrechnung).
pub const MAX_SUBSTEPS: usize = 2_000_000;

/// Hauptachsen eines Trägheitstensors.
#[derive(Debug, Clone, Copy)]
pub struct PrincipalAxes {
    /// Hauptträgheitsmomente, aufsteigend.
    pub moments: DVec3,
    /// Spalten = Hauptachsen im Schiffssystem (Rechtssystem).
    pub axes: DMat3,
}

/// Eigenzerlegung des symmetrischen Tensors (zyklisches Jacobi-Verfahren).
pub fn principal_axes(inertia: DMat3) -> PrincipalAxes {
    let mut a = inertia.to_cols_array_2d(); // a[col][row], symmetrisch
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]; // v[col][row]
    let scale = (0..3).map(|i| a[i][i].abs()).fold(0.0, f64::max).max(f64::MIN_POSITIVE);
    for _ in 0..50 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off <= 1e-15 * scale {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q].abs() <= 1e-300 {
                continue;
            }
            // Drehwinkel, der a[p][q] zu null macht
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            // A ← Jᵀ A J
            for col in a.iter_mut() {
                let (akp, akq) = (col[p], col[q]);
                col[p] = c * akp - s * akq;
                col[q] = s * akp + c * akq;
            }
            let (ap, aq) = (a[p], a[q]);
            for k in 0..3 {
                a[p][k] = c * ap[k] - s * aq[k];
                a[q][k] = s * ap[k] + c * aq[k];
            }
            // V ← V J
            let (vp, vq) = (v[p], v[q]);
            for k in 0..3 {
                v[p][k] = c * vp[k] - s * vq[k];
                v[q][k] = s * vp[k] + c * vq[k];
            }
        }
    }
    let mut idx = [0, 1, 2];
    idx.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let col = |i: usize| DVec3::from_array(v[i]);
    let (x, y) = (col(idx[0]), col(idx[1]));
    // Rechtssystem erzwingen
    let z = x.cross(y);
    PrincipalAxes { moments: DVec3::new(a[idx[0]][idx[0]], a[idx[1]][idx[1]], a[idx[2]][idx[2]]), axes: DMat3::from_cols(x, y, z) }
}

/// Ergebnis eines Drehschritts.
#[derive(Debug, Clone, Copy)]
pub struct Spin {
    pub attitude: Rotation<VesselFrame, Inertial>,
    pub omega: Vec3<VesselFrame>,
}

/// Freie Drehung über `dt`. Mit `brake = Some(τ)` bremsen Reaktionsräder
/// mit höchstens τ (N·m) je Schiffsachse gegen den Drehimpuls.
pub fn drift(
    attitude: Rotation<VesselFrame, Inertial>,
    omega: Vec3<VesselFrame>,
    inertia: &InertiaTensor,
    dt: f64,
    brake: Option<f64>,
) -> Spin {
    drift_with(attitude, omega, inertia, dt, brake, MAX_SUBSTEP_ANGLE)
}

pub(crate) fn drift_with(
    attitude: Rotation<VesselFrame, Inertial>,
    omega: Vec3<VesselFrame>,
    inertia: &InertiaTensor,
    dt: f64,
    brake: Option<f64>,
    max_angle: f64,
) -> Spin {
    let brake = brake.filter(|&t| t > 0.0);
    if dt <= 0.0 || (omega == Vec3::ZERO && brake.is_none()) {
        return Spin { attitude, omega };
    }
    let pa = principal_axes(inertia.0);
    let r = pa.axes;
    let q_r = DQuat::from_mat3(&r);
    let m = pa.moments;
    // Lage der Hauptachsen im Inertialsystem, Drehimpuls in Hauptachsen
    let mut q = attitude.quat() * q_r;
    let mut l = r.transpose() * (inertia.0 * omega.raw());

    // Teilschrittweite aus einer Erhaltungsgröße: |L|/I_min ist die größte
    // mögliche Drehrate bei diesem Drehimpuls. Hinge sie vom momentanen ω
    // ab, wechselte die Schrittweite entlang der Bahn, und die Energie
    // würde langsam wegdriften.
    let w = l.length() / m.x;
    let n = ((w * dt / max_angle).ceil() as usize).clamp(1, MAX_SUBSTEPS);
    let h = dt / n as f64;
    let axes = [DVec3::X, DVec3::Y, DVec3::Z];
    let sub = |k: usize, tau: f64, q: &mut DQuat, l: &mut DVec3| {
        let theta = l[k] / m[k] * tau;
        if theta != 0.0 {
            let rot = DQuat::from_axis_angle(axes[k], theta);
            *q *= rot;
            *l = rot.inverse() * *l;
        }
    };
    for i in 0..n {
        sub(0, h / 2.0, &mut q, &mut l);
        sub(1, h / 2.0, &mut q, &mut l);
        sub(2, h, &mut q, &mut l);
        sub(1, h / 2.0, &mut q, &mut l);
        sub(0, h / 2.0, &mut q, &mut l);
        if let Some(tau) = brake {
            // Radmoment je Schiffsachse, höchstens bis zum Stillstand
            let lv = r * l;
            let limit = tau * h;
            let lv = DVec3::new(
                lv.x - lv.x.clamp(-limit, limit),
                lv.y - lv.y.clamp(-limit, limit),
                lv.z - lv.z.clamp(-limit, limit),
            );
            l = r.transpose() * lv;
            if lv == DVec3::ZERO {
                break;
            }
        }
        if i % 1024 == 1023 {
            q = q.normalize();
        }
    }
    let q = q.normalize();
    let attitude = Rotation::from_quat(q * q_r.inverse());
    let omega = inertia.solve(Vec3::from_raw(r * l));
    Spin { attitude, omega }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tensor(d: DVec3) -> InertiaTensor {
        InertiaTensor(DMat3::from_diagonal(d))
    }

    fn l_inertial(s: &Spin, i: &InertiaTensor) -> DVec3 {
        s.attitude.apply(i.apply(s.omega)).raw()
    }

    fn energy(s: &Spin, i: &InertiaTensor) -> f64 {
        0.5 * s.omega.dot(i.apply(s.omega))
    }

    fn start(omega: DVec3) -> Spin {
        Spin { attitude: Rotation::IDENTITY, omega: Vec3::from_raw(omega) }
    }

    fn run(mut s: Spin, i: &InertiaTensor, dt: f64, steps: usize) -> Spin {
        for _ in 0..steps {
            s = drift(s.attitude, s.omega, i, dt, None);
        }
        s
    }

    #[test]
    fn principal_axes_of_rotated_tensor() {
        let d = DVec3::new(3.0, 1.0, 2.0);
        let rot = DMat3::from_quat(DQuat::from_euler(glam::EulerRot::XYZ, 0.3, -0.7, 1.1));
        let i = rot * DMat3::from_diagonal(d) * rot.transpose();
        let pa = principal_axes(i);
        assert!((pa.moments - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-12, "{:?}", pa.moments);
        // I·a_k = λ_k·a_k und Rechtssystem
        for k in 0..3 {
            let a = pa.axes.col(k);
            assert!((i * a - a * pa.moments[k]).length() < 1e-12);
        }
        assert!((pa.axes.determinant() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dzhanibekov_intermediate_axis_flips() {
        // Drehung fast genau um die mittlere Achse: kippt periodisch um.
        let i = tensor(DVec3::new(1.0, 2.0, 3.0));
        let mut s = start(DVec3::new(1e-3, 1.0, 1e-3));
        let l0 = l_inertial(&s, &i);
        let e0 = energy(&s, &i);
        let mut sign = 1.0;
        let mut flips = 0;
        for _ in 0..2000 {
            s = drift(s.attitude, s.omega, &i, 0.1, None);
            if s.omega.y() * sign < -0.5 {
                sign = -sign;
                flips += 1;
            }
        }
        assert!(flips >= 3, "nur {flips} Umschläge in 200 s");
        assert!((l_inertial(&s, &i) - l0).length() < 1e-12 * l0.length());
        assert!((energy(&s, &i) - e0).abs() < 1e-3 * e0);
    }

    #[test]
    fn major_and_minor_axis_are_stable() {
        let i = tensor(DVec3::new(1.0, 2.0, 3.0));
        for (axis, spin) in [(2, DVec3::new(1e-3, 1e-3, 1.0)), (0, DVec3::new(1.0, 1e-3, 1e-3))] {
            let mut s = start(spin);
            for _ in 0..2000 {
                s = drift(s.attitude, s.omega, &i, 0.1, None);
                let w = s.omega.raw();
                assert!(w[axis].abs() > 0.99 * w.length(), "Achse {axis} instabil: {w:?}");
            }
        }
    }

    #[test]
    fn no_drift_over_a_million_seconds() {
        // Allgemeiner Fall, langer Zeitraffer: 1000 Schritte à 1000 s
        let i = tensor(DVec3::new(1.0, 2.0, 3.0));
        let s0 = start(DVec3::new(0.3, 0.2, 0.25));
        let (l0, e0) = (l_inertial(&s0, &i), energy(&s0, &i));
        let err = |s: &Spin| (energy(s, &i) - e0).abs() / e0;
        // Größte Abweichung in den ersten 1000 s (Schwingung des Fehlers)
        let mut s = s0;
        let mut early = 0.0f64;
        for _ in 0..1000 {
            s = drift(s.attitude, s.omega, &i, 1.0, None);
            early = early.max(err(&s));
        }
        // ... und über 10⁶ s in Schritten von 1000 s
        let mut s = s0;
        let mut late = 0.0f64;
        for _ in 0..1000 {
            s = drift(s.attitude, s.omega, &i, 1000.0, None);
            late = late.max(err(&s));
        }
        assert!((l_inertial(&s, &i) - l0).length() < 1e-9 * l0.length());
        // Energiefehler beschränkt, wächst nicht mit der Zeit
        assert!(late < 1e-5 && late < 1.5 * early, "{early:e} → {late:e}");
    }

    #[test]
    fn long_step_equals_many_short_steps() {
        // Zeitraffer: ein Schritt von 100 s ergibt dasselbe wie 1000 à 0,1 s
        let i = tensor(DVec3::new(1.0, 2.0, 3.0));
        let s0 = start(DVec3::new(0.3, 0.2, 0.25));
        let a = run(s0, &i, 100.0, 1);
        let b = run(s0, &i, 0.1, 1000);
        assert!(a.attitude.quat().angle_between(b.attitude.quat()) < 0.02);
        assert!((a.omega - b.omega).length() < 0.01 * s0.omega.length());
    }

    #[test]
    fn second_order_convergence() {
        // Unsymmetrischer, voll besetzter Tensor gegen sehr feine Referenz
        let rot = DMat3::from_quat(DQuat::from_euler(glam::EulerRot::XYZ, 0.4, 0.2, -0.5));
        let i = InertiaTensor(rot * DMat3::from_diagonal(DVec3::new(2.0, 5.0, 6.0)) * rot.transpose());
        let w = Vec3::new(0.4, -0.3, 0.5);
        let reference = drift_with(Rotation::IDENTITY, w, &i, 20.0, None, 1e-4);
        let err = |angle: f64| {
            let s = drift_with(Rotation::IDENTITY, w, &i, 20.0, None, angle);
            s.attitude.quat().angle_between(reference.attitude.quat())
        };
        let (e1, e2) = (err(0.08), err(0.04));
        let ratio = e1 / e2;
        assert!(ratio > 3.0 && ratio < 5.0, "Fehler {e1:e} / {e2:e} = {ratio}");
    }

    #[test]
    fn matches_euler_equations_with_full_tensor() {
        // Unabhängige Prüfung: RK4 auf L̇ = L × ω und q̇ = ½·q·ω im Schiffssystem
        let rot = DMat3::from_quat(DQuat::from_euler(glam::EulerRot::XYZ, -0.6, 0.9, 0.3));
        let i = InertiaTensor(rot * DMat3::from_diagonal(DVec3::new(1.5, 4.0, 4.5)) * rot.transpose());
        let w0 = DVec3::new(0.2, 0.5, -0.3);
        let inv = i.0.inverse();
        let deriv = |q: DQuat, l: DVec3| {
            let w = inv * l;
            (q * DQuat::from_xyzw(w.x, w.y, w.z, 0.0) * 0.5, l.cross(w))
        };
        let (mut q, mut l) = (DQuat::IDENTITY, i.0 * w0);
        let h = 1e-3;
        for _ in 0..30_000 {
            let (k1q, k1l) = deriv(q, l);
            let (k2q, k2l) = deriv(q + k1q * (h / 2.0), l + k1l * (h / 2.0));
            let (k3q, k3l) = deriv(q + k2q * (h / 2.0), l + k2l * (h / 2.0));
            let (k4q, k4l) = deriv(q + k3q * h, l + k3l * h);
            q = (q + (k1q + k2q * 2.0 + k3q * 2.0 + k4q) * (h / 6.0)).normalize();
            l += (k1l + k2l * 2.0 + k3l * 2.0 + k4l) * (h / 6.0);
        }
        let s = drift_with(Rotation::IDENTITY, Vec3::from_raw(w0), &i, 30.0, None, 0.005);
        assert!(s.attitude.quat().angle_between(q) < 1e-3, "{}", s.attitude.quat().angle_between(q));
        assert!((s.omega.raw() - inv * l).length() < 1e-3 * w0.length());
    }

    #[test]
    fn wheels_brake_at_their_torque_limit() {
        // |L| = 10 N·m·s um x, Rad 0,5 N·m: Stillstand nach 20 s, unabhängig vom Schritt
        let i = tensor(DVec3::new(10.0, 20.0, 30.0));
        for dt in [0.1, 5.0] {
            let mut s = start(DVec3::new(1.0, 0.0, 0.0));
            let mut t = 0.0;
            while s.omega.length() > 0.0 && t < 100.0 {
                s = drift(s.attitude, s.omega, &i, dt, Some(0.5));
                t += dt;
                if (t - 10.0).abs() < 1e-9 {
                    assert!((s.omega.length() - 0.5).abs() < 1e-9, "dt {dt}: {}", s.omega.length());
                }
            }
            assert!((t - 20.0).abs() <= dt + 1e-9, "dt {dt}: erst nach {t} s gebremst");
        }
        // Ein einziger langer Schritt endet in Ruhe
        let s = drift(Rotation::IDENTITY, Vec3::new(1.0, 0.0, 0.0), &i, 1000.0, Some(0.5));
        assert_eq!(s.omega, Vec3::ZERO);
    }
}
