//! Typisierte Bezugssysteme.
//!
//! Ein großer Teil der "Kraken"-Fehler in KSP entsteht, weil Vektoren aus
//! unterschiedlichen Bezugssystemen (Weltursprung, verschobener Ursprung,
//! mitbewegtes Krakensbane-System, Schiffssystem) versehentlich verrechnet
//! werden. Hier ist das ein Compile-Fehler:
//!
//! ```compile_fail
//! use sim_core::{Vec3, Inertial, VesselFrame};
//! let a: Vec3<Inertial> = Vec3::new(1.0, 0.0, 0.0);
//! let b: Vec3<VesselFrame> = Vec3::new(0.0, 1.0, 0.0);
//! let _ = a + b; // Fehler: verschiedene Bezugssysteme
//! ```
//!
//! Umrechnen geht nur explizit:
//!
//! ```
//! use sim_core::{Vec3, Inertial, VesselFrame, Rotation};
//! let att: Rotation<VesselFrame, Inertial> = Rotation::IDENTITY;
//! let a: Vec3<Inertial> = Vec3::new(1.0, 0.0, 0.0);
//! let b: Vec3<VesselFrame> = Vec3::new(0.0, 1.0, 0.0);
//! let _ = a + att.apply(b);
//! ```

use glam::{DQuat, DVec3};
use std::fmt;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// Marker für ein Bezugssystem.
pub trait Frame: Copy + Clone + fmt::Debug + PartialEq + 'static {
    const NAME: &'static str;
}

/// Inertialsystem mit Ursprung im Mittelpunkt des aktuellen Zentralkörpers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Inertial;
impl Frame for Inertial {
    const NAME: &'static str = "Inertial";
}

/// Mit dem Zentralkörper mitrotierendes System, Ursprung im Mittelpunkt,
/// +Z = Drehachse. Startplatz, Luft und Boden ruhen hier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyFixed;
impl Frame for BodyFixed {
    const NAME: &'static str = "Körperfest";
}

/// Schiffsfestes System. Ursprung: Mitte des Wurzelteils, +Z = Spitze.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VesselFrame;
impl Frame for VesselFrame {
    const NAME: &'static str = "Schiff";
}

/// Lokales System eines Bauteils (aus der Teiledefinition).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartLocal;
impl Frame for PartLocal {
    const NAME: &'static str = "Teil";
}

/// 3D-Vektor (f64) in einem bestimmten Bezugssystem.
#[derive(Clone, Copy, PartialEq)]
pub struct Vec3<F: Frame> {
    raw: DVec3,
    _frame: PhantomData<F>,
}

impl<F: Frame> fmt::Debug for Vec3<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}({:.6}, {:.6}, {:.6})", F::NAME, self.raw.x, self.raw.y, self.raw.z)
    }
}

impl<F: Frame> Default for Vec3<F> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<F: Frame> Vec3<F> {
    pub const ZERO: Self = Self::from_raw(DVec3::ZERO);
    pub const X: Self = Self::from_raw(DVec3::X);
    pub const Y: Self = Self::from_raw(DVec3::Y);
    pub const Z: Self = Self::from_raw(DVec3::Z);

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self::from_raw(DVec3::new(x, y, z))
    }

    /// Ausstieg ins untypisierte glam. Bewusst explizit benannt:
    /// wer das benutzt, übernimmt die Verantwortung für das Bezugssystem.
    pub const fn from_raw(raw: DVec3) -> Self {
        Self { raw, _frame: PhantomData }
    }

    pub const fn raw(self) -> DVec3 {
        self.raw
    }

    pub fn x(self) -> f64 {
        self.raw.x
    }
    pub fn y(self) -> f64 {
        self.raw.y
    }
    pub fn z(self) -> f64 {
        self.raw.z
    }

    pub fn dot(self, o: Self) -> f64 {
        self.raw.dot(o.raw)
    }
    pub fn cross(self, o: Self) -> Self {
        Self::from_raw(self.raw.cross(o.raw))
    }
    pub fn length(self) -> f64 {
        self.raw.length()
    }
    pub fn length_squared(self) -> f64 {
        self.raw.length_squared()
    }
    /// Einheitsvektor; Nullvektor bleibt Nullvektor.
    pub fn normalize_or_zero(self) -> Self {
        Self::from_raw(self.raw.normalize_or_zero())
    }
    /// Komponente senkrecht zu `axis` (axis muss normiert sein).
    pub fn reject_from_unit(self, axis: Self) -> Self {
        self - axis * self.dot(axis)
    }
    pub fn is_finite(self) -> bool {
        self.raw.is_finite()
    }
}

impl<F: Frame> Add for Vec3<F> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::from_raw(self.raw + o.raw)
    }
}
impl<F: Frame> Sub for Vec3<F> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::from_raw(self.raw - o.raw)
    }
}
impl<F: Frame> Neg for Vec3<F> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::from_raw(-self.raw)
    }
}
impl<F: Frame> Mul<f64> for Vec3<F> {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::from_raw(self.raw * s)
    }
}
impl<F: Frame> Mul<Vec3<F>> for f64 {
    type Output = Vec3<F>;
    fn mul(self, v: Vec3<F>) -> Vec3<F> {
        v * self
    }
}
impl<F: Frame> Div<f64> for Vec3<F> {
    type Output = Self;
    fn div(self, s: f64) -> Self {
        Self::from_raw(self.raw / s)
    }
}
impl<F: Frame> AddAssign for Vec3<F> {
    fn add_assign(&mut self, o: Self) {
        self.raw += o.raw;
    }
}
impl<F: Frame> SubAssign for Vec3<F> {
    fn sub_assign(&mut self, o: Self) {
        self.raw -= o.raw;
    }
}
impl<F: Frame> std::iter::Sum for Vec3<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |a, b| a + b)
    }
}

/// Rotation, die Vektoren aus `From` nach `To` abbildet.
#[derive(Clone, Copy, PartialEq)]
pub struct Rotation<From: Frame, To: Frame> {
    q: DQuat,
    _f: PhantomData<(From, To)>,
}

impl<A: Frame, B: Frame> fmt::Debug for Rotation<A, B> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Rotation<{}→{}>({:?})", A::NAME, B::NAME, self.q)
    }
}

impl<A: Frame, B: Frame> Rotation<A, B> {
    pub const IDENTITY: Self = Self::from_quat(DQuat::IDENTITY);

    /// Untypisierter Einstieg, siehe [`Vec3::from_raw`].
    pub const fn from_quat(q: DQuat) -> Self {
        Self { q, _f: PhantomData }
    }
    pub fn quat(self) -> DQuat {
        self.q
    }
    pub fn apply(self, v: Vec3<A>) -> Vec3<B> {
        Vec3::from_raw(self.q * v.raw())
    }
    pub fn inverse(self) -> Rotation<B, A> {
        Rotation::from_quat(self.q.inverse())
    }
    /// Gegen numerische Drift nach vielen Integrationsschritten.
    pub fn renormalize(self) -> Self {
        Self::from_quat(self.q.normalize())
    }
    /// Rotationsmatrix (Spalten = Basisvektoren von `A` in `B`).
    pub fn matrix(self) -> glam::DMat3 {
        glam::DMat3::from_quat(self.q)
    }
}

impl<F: Frame> Rotation<F, F> {
    /// Drehung um `axis * angle` (Drehvektor) innerhalb eines Systems.
    pub fn from_scaled_axis(v: Vec3<F>) -> Self {
        Self::from_quat(DQuat::from_scaled_axis(v.raw()))
    }
    pub fn from_axis_angle(axis: Vec3<F>, angle: f64) -> Self {
        Self::from_quat(DQuat::from_axis_angle(axis.normalize_or_zero().raw(), angle))
    }
}

/// Verkettung: (B→C) * (A→B) = (A→C). Falsche Reihenfolge kompiliert nicht.
impl<A: Frame, B: Frame, C: Frame> Mul<Rotation<A, B>> for Rotation<B, C> {
    type Output = Rotation<A, C>;
    fn mul(self, rhs: Rotation<A, B>) -> Rotation<A, C> {
        Rotation::from_quat(self.q * rhs.q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn rotation_roundtrip_and_composition() {
        let r1: Rotation<PartLocal, VesselFrame> =
            Rotation::from_quat(DQuat::from_rotation_z(FRAC_PI_2));
        let r2: Rotation<VesselFrame, Inertial> =
            Rotation::from_quat(DQuat::from_rotation_x(FRAC_PI_2));
        let v: Vec3<PartLocal> = Vec3::new(1.0, 0.0, 0.0);

        let via_steps = r2.apply(r1.apply(v));
        let composed = (r2 * r1).apply(v);
        assert!((via_steps - composed).length() < 1e-12);

        let back = r1.inverse().apply(r1.apply(v));
        assert!((back - v).length() < 1e-12);
        // x -> y (um z) -> z (um x)
        assert!((composed - Vec3::new(0.0, 0.0, 1.0)).length() < 1e-12);
    }
}
