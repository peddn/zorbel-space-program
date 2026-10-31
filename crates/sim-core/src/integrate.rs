//! Numerische Integration.
//!
//! Das Hauptverfahren ist [`kepler_kdk`]: die Bewegung im Schwerefeld
//! ("Drift") wird exakt über die Keplerbahn gelöst, alle anderen Kräfte
//! (Schub, Luft) wirken als Stöße ("Kicks") vor und nach der Drift.
//! Sind die anderen Kräfte null, ist das Ergebnis exakt die Keplerbahn.
//! Es gibt also keine getrennten Modi "auf Schienen" und "mit Physik",
//! sondern ein Verfahren, dessen Fehler mit den Störkräften verschwindet.

use crate::frame::{Frame, Vec3};
use crate::kepler::{self, StateVector};

/// Ein Kick-Drift-Kick-Schritt (2. Ordnung). `perturbation` liefert die
/// Beschleunigung ohne Zentralgravitation.
pub fn kepler_kdk<F: Frame>(
    s: StateVector<F>,
    mu: f64,
    dt: f64,
    perturbation: impl Fn(&StateVector<F>) -> Vec3<F>,
) -> StateVector<F> {
    let mut s = s;
    s.v += perturbation(&s) * (dt / 2.0);
    s = kepler::propagate(s, mu, dt);
    s.v += perturbation(&s) * (dt / 2.0);
    s
}

/// Gravitation eines Punktmassen-Körpers am Ort `r` (relativ zum Mittelpunkt).
pub fn point_mass_gravity<F: Frame>(r: Vec3<F>, mu: f64) -> Vec3<F> {
    let d2 = r.length_squared();
    -r * (mu / (d2 * d2.sqrt()))
}

/// Ein klassischer Runge-Kutta-Schritt 4. Ordnung.
/// `accel` liefert die Gesamtbeschleunigung für einen Zustand.
pub fn rk4_step<F: Frame>(
    s: StateVector<F>,
    dt: f64,
    accel: impl Fn(&StateVector<F>) -> Vec3<F>,
) -> StateVector<F> {
    let k1v = accel(&s);
    let k1r = s.v;

    let s2 = StateVector::new(s.r + k1r * (dt / 2.0), s.v + k1v * (dt / 2.0));
    let k2v = accel(&s2);
    let k2r = s2.v;

    let s3 = StateVector::new(s.r + k2r * (dt / 2.0), s.v + k2v * (dt / 2.0));
    let k3v = accel(&s3);
    let k3r = s3.v;

    let s4 = StateVector::new(s.r + k3r * dt, s.v + k3v * dt);
    let k4v = accel(&s4);
    let k4r = s4.v;

    StateVector::new(
        s.r + (k1r + k2r * 2.0 + k3r * 2.0 + k4r) * (dt / 6.0),
        s.v + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (dt / 6.0),
    )
}
