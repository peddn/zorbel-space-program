//! Analytische Zweikörper-Bahnen ("on rails").
//!
//! Propagation über universelle Variablen (Stumpff-Funktionen), damit
//! elliptische, parabolische und hyperbolische Bahnen denselben Codepfad
//! nehmen. Das Ergebnis ist exakt bis auf Rundung, unabhängig von der
//! Schrittweite; deshalb ist beliebiger Zeitraffer ohne Drift möglich.

use crate::frame::{Frame, Vec3};
use std::f64::consts::{PI, TAU};

/// Position und Geschwindigkeit relativ zum Zentralkörper.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateVector<F: Frame> {
    pub r: Vec3<F>,
    pub v: Vec3<F>,
}

impl<F: Frame> StateVector<F> {
    pub fn new(r: Vec3<F>, v: Vec3<F>) -> Self {
        Self { r, v }
    }
    /// Spezifische Bahnenergie ε = v²/2 − μ/r (Erhaltungsgröße).
    pub fn specific_energy(&self, mu: f64) -> f64 {
        0.5 * self.v.length_squared() - mu / self.r.length()
    }
    /// Spezifischer Drehimpuls h = r × v (Erhaltungsgröße).
    pub fn angular_momentum(&self) -> Vec3<F> {
        self.r.cross(self.v)
    }
}

/// Stumpff-Funktion C(z).
fn stumpff_c(z: f64) -> f64 {
    if z > 1e-6 {
        let s = z.sqrt();
        (1.0 - s.cos()) / z
    } else if z < -1e-6 {
        let s = (-z).sqrt();
        (s.cosh() - 1.0) / (-z)
    } else {
        1.0 / 2.0 - z / 24.0 + z * z / 720.0 - z * z * z / 40_320.0
    }
}

/// Stumpff-Funktion S(z).
fn stumpff_s(z: f64) -> f64 {
    if z > 1e-6 {
        let s = z.sqrt();
        (s - s.sin()) / (s * s * s)
    } else if z < -1e-6 {
        let s = (-z).sqrt();
        (s.sinh() - s) / (s * s * s)
    } else {
        1.0 / 6.0 - z / 120.0 + z * z / 5_040.0 - z * z * z / 362_880.0
    }
}

/// Zustand nach `dt` Sekunden auf der Keplerbahn (dt darf negativ sein).
pub fn propagate<F: Frame>(s: StateVector<F>, mu: f64, dt: f64) -> StateVector<F> {
    if dt == 0.0 {
        return s;
    }
    let r0v = s.r;
    let v0v = s.v;
    let r0 = r0v.length();
    let v0 = v0v.length();
    let sqrt_mu = mu.sqrt();
    let rdotv = r0v.dot(v0v);
    let alpha = 2.0 / r0 - v0 * v0 / mu; // 1/a

    // Auf Ellipsen ganze Umläufe abziehen: hält χ klein und präzise.
    let mut dt = dt;
    if alpha > 1e-12 {
        let period = TAU / (alpha.powi(3) * mu).sqrt();
        dt %= period;
    }

    // Startwert nach Vallado.
    let mut chi = if alpha > 1e-12 {
        sqrt_mu * dt * alpha
    } else if alpha < -1e-12 && dt.abs() * v0 > r0 {
        // Schätzung nach Vallado; gilt für lange Schritte. Bei kurzen Schritten
        // (Strecke kleiner als der Abstand) führt sie Newton zur falschen
        // Lösung, dort ist die lineare Schätzung unten richtig.
        let a = 1.0 / alpha;
        let sign = dt.signum();
        let num = -2.0 * mu * alpha * dt;
        let den = rdotv + sign * (-mu * a).sqrt() * (1.0 - r0 * alpha);
        let guess = sign * (-a).sqrt() * (num / den).ln();
        if guess.is_finite() { guess } else { sqrt_mu * dt / r0 }
    } else {
        sqrt_mu * dt / r0
    };

    // Newton-Verfahren auf der universellen Kepler-Gleichung.
    let c1 = rdotv / sqrt_mu;
    for _ in 0..100 {
        let z = alpha * chi * chi;
        let c = stumpff_c(z);
        let ss = stumpff_s(z);
        let chi2 = chi * chi;
        let f = c1 * chi2 * c + (1.0 - alpha * r0) * chi2 * chi * ss + r0 * chi - sqrt_mu * dt;
        let df = c1 * chi * (1.0 - z * ss) + (1.0 - alpha * r0) * chi2 * c + r0;
        let step = f / df;
        chi -= step;
        if step.abs() < 1e-12 * chi.abs().max(1.0) {
            break;
        }
    }

    let z = alpha * chi * chi;
    let c = stumpff_c(z);
    let ss = stumpff_s(z);
    let chi2 = chi * chi;

    let f = 1.0 - chi2 / r0 * c;
    let g = dt - chi2 * chi / sqrt_mu * ss;
    let r = r0v * f + v0v * g;
    let rn = r.length();
    let fdot = sqrt_mu / (rn * r0) * (alpha * chi2 * chi * ss - chi);
    let gdot = 1.0 - chi2 / rn * c;
    let v = r0v * fdot + v0v * gdot;
    StateVector { r, v }
}

/// Frühester Zeitpunkt in `[0, t_max]`, zu dem die Bahn den Radius
/// `r_limit` unterschreitet. Dient der Schrittweitensteuerung: ein großer
/// Schritt soll nicht über den Eintritt in die Atmosphäre hinwegspringen.
pub fn time_until_below<F: Frame>(s: StateVector<F>, mu: f64, r_limit: f64, t_max: f64) -> Option<f64> {
    if s.r.length() <= r_limit {
        return Some(0.0);
    }
    let el = elements(s, mu);
    if el.periapsis >= r_limit {
        return None;
    }
    // Offene Bahn, die sich schon entfernt: kommt nie zurück.
    let bound = s.specific_energy(mu) < 0.0;
    if !bound && s.r.dot(s.v) >= 0.0 {
        return None;
    }
    // Spätestens nach einem Umlauf ist ein Unterschreiten passiert.
    let t_end = match el.period {
        Some(p) if bound => t_max.min(p),
        _ => t_max,
    };

    // Vorwärts tasten in Schritten, die an die lokale Zeitskala r/v
    // angepasst sind; Periapsiszeit als zusätzlicher Prüfpunkt, damit kein
    // flacher Durchgang zwischen zwei Proben verloren geht. Robust auch für
    // entartete (radiale) Bahnen.
    let r_at = |t: f64| propagate(s, mu, t).r.length();
    let t_peri = el.time_to_periapsis;
    let mut t = 0.0;
    while t < t_end {
        let st = propagate(s, mu, t);
        let r = st.r.length();
        let scale = r / st.v.length().max((mu / r).sqrt());
        let mut next = (t + 0.05 * scale).min(t_end);
        if t_peri.is_finite() && t_peri > t && t_peri < next {
            next = t_peri;
        }
        if r_at(next) <= r_limit {
            let (mut lo, mut hi) = (t, next);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if r_at(mid) > r_limit { lo = mid } else { hi = mid }
            }
            return Some(hi);
        }
        t = next;
    }
    None
}

/// Bahnelemente zur Anzeige und für Manöverlogik.
#[derive(Debug, Clone, Copy)]
pub struct Elements {
    /// Große Halbachse (negativ bei Hyperbeln, ∞ bei Parabeln).
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    /// Bahnneigung in rad.
    pub inclination: f64,
    /// Periapsis-Radius (vom Mittelpunkt).
    pub periapsis: f64,
    /// Apoapsis-Radius; `None` bei offenen Bahnen.
    pub apoapsis: Option<f64>,
    /// Umlaufzeit; `None` bei offenen Bahnen.
    pub period: Option<f64>,
    /// Wahre Anomalie in rad (0..2π).
    pub true_anomaly: f64,
    /// Zeit bis zur nächsten Apoapsis; `None` bei offenen Bahnen.
    pub time_to_apoapsis: Option<f64>,
    /// Zeit bis zur nächsten Periapsis (bei Hyperbeln negativ, falls schon vorbei).
    pub time_to_periapsis: f64,
}

pub fn elements<F: Frame>(s: StateVector<F>, mu: f64) -> Elements {
    let r = s.r.length();
    let v2 = s.v.length_squared();
    let h = s.angular_momentum();
    let hn = h.length();
    let rdotv = s.r.dot(s.v);
    let e_vec = (s.r * (v2 - mu / r) - s.v * rdotv) / mu;
    let e = e_vec.length();
    let energy = 0.5 * v2 - mu / r;
    let a = if energy.abs() > 1e-12 { -mu / (2.0 * energy) } else { f64::INFINITY };
    let inclination = (h.z() / hn).clamp(-1.0, 1.0).acos();
    let periapsis = hn * hn / mu / (1.0 + e);

    // Wahre Anomalie; bei (fast) Kreisbahnen ist die Periapsis-Richtung
    // undefiniert, dann zählen wir von der aktuellen Position.
    let mut nu = if e > 1e-10 {
        let cos_nu = (e_vec.dot(s.r) / (e * r)).clamp(-1.0, 1.0);
        let nu = cos_nu.acos();
        if rdotv < 0.0 { TAU - nu } else { nu }
    } else {
        0.0
    };
    if !nu.is_finite() {
        nu = 0.0;
    }

    // Gebunden oder nicht entscheidet die Energie, nicht e: eine radiale
    // Bahn (z. B. Start aus der Ruhe) hat e = 1, ist aber gebunden.
    let (apoapsis, period, time_to_apoapsis, time_to_periapsis) = if energy < 0.0 {
        let n = (mu / a.powi(3)).sqrt();
        let e = e.min(1.0 - 1e-12);
        // Exzentrische und mittlere Anomalie.
        let ecc_anom = 2.0 * (((1.0 - e) / (1.0 + e)).sqrt() * (nu / 2.0).tan()).atan();
        let m = (ecc_anom - e * ecc_anom.sin()).rem_euclid(TAU);
        let to_apo = (PI - m).rem_euclid(TAU) / n;
        let to_peri = (TAU - m).rem_euclid(TAU) / n;
        (Some(a * (1.0 + e)), Some(TAU / n), Some(to_apo), to_peri)
    } else if e > 1.0 && energy > 0.0 {
        let n = (mu / (-a).powi(3)).sqrt();
        let nu_s = if nu > PI { nu - TAU } else { nu };
        let hyp_anom = 2.0 * (((e - 1.0) / (e + 1.0)).sqrt() * (nu_s / 2.0).tan()).atanh();
        let m = e * hyp_anom.sinh() - hyp_anom;
        (None, None, None, -m / n)
    } else {
        (None, None, None, f64::NAN)
    };

    Elements {
        semi_major_axis: a,
        eccentricity: e,
        inclination,
        periapsis,
        apoapsis,
        period,
        true_anomaly: nu,
        time_to_apoapsis,
        time_to_periapsis,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Inertial;
    use crate::integrate::{point_mass_gravity, rk4_step};

    const MU: f64 = 3.5316e12; // Testplanet, 600 km Radius

    fn ellipse() -> StateVector<Inertial> {
        // Periapsis 700 km, deutlich exzentrisch
        let r = 700_000.0;
        let v = (MU / r).sqrt() * 1.25;
        StateVector::new(Vec3::new(r, 0.0, 0.0), Vec3::new(0.0, v * 0.96, v * 0.28))
    }

    fn rel_err(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs().max(1e-30)
    }

    #[test]
    fn full_period_returns_to_start() {
        for s in [
            StateVector::new(
                Vec3::<Inertial>::new(700_000.0, 0.0, 0.0),
                Vec3::new(0.0, (MU / 700_000.0).sqrt(), 0.0),
            ),
            ellipse(),
        ] {
            let el = elements(s, MU);
            let p = el.period.unwrap();
            let back = propagate(s, MU, p);
            assert!((back.r - s.r).length() < 1e-3, "Positionsfehler {:?}", back.r - s.r);
            assert!((back.v - s.v).length() < 1e-6);
        }
    }

    #[test]
    fn energy_and_angular_momentum_conserved() {
        let s0 = ellipse();
        let e0 = s0.specific_energy(MU);
        let h0 = s0.angular_momentum();
        // 10 000 Umläufe in krummen Schritten, ohne Fehlerakkumulation
        let p = elements(s0, MU).period.unwrap();
        let mut s = s0;
        for i in 0..2_000 {
            s = propagate(s, MU, p * 5.0 + 1.37 * i as f64);
            assert!(rel_err(s.specific_energy(MU), e0) < 1e-9);
            assert!((s.angular_momentum() - h0).length() / h0.length() < 1e-9);
        }
    }

    #[test]
    fn hyperbola_short_steps_far_out() {
        // Weit draußen auf einer Fluchtbahn: kurze und lange Schritte müssen
        // mit der Geschwindigkeit zusammenpassen (früher Fehler bei 1 s).
        let r0 = 700_000.0;
        let s0 = StateVector::<Inertial>::new(Vec3::new(r0, 0.0, 0.0), Vec3::new(0.0, (MU / r0).sqrt() * 2f64.sqrt() * 1.3, 0.0));
        let far = propagate(s0, MU, 30_000.0);
        for dt in [0.01, 1.0, 10.0, 100.0, 1_000.0, 1.0e5] {
            let a = propagate(far, MU, dt);
            let b = (0..10).fold(far, |s, _| propagate(s, MU, dt / 10.0));
            assert!((a.r - b.r).length() < 1e-6 * a.r.length(), "dt {dt}: {:?}", a.r - b.r);
            if dt <= 100.0 {
                let moved = (a.r - far.r).length();
                assert!((moved - far.v.length() * dt).abs() < 1e-3 * moved, "dt {dt}: {moved} m");
            }
        }
    }

    #[test]
    fn hyperbola_forward_and_back() {
        let r = 700_000.0;
        let v_esc = (2.0 * MU / r).sqrt();
        let s0 = StateVector::<Inertial>::new(Vec3::new(r, 0.0, 0.0), Vec3::new(0.0, v_esc * 1.3, 0.0));
        assert!(elements(s0, MU).eccentricity > 1.0);
        let fwd = propagate(s0, MU, 20_000.0);
        assert!(fwd.r.length() > 10.0 * r);
        assert!(rel_err(fwd.specific_energy(MU), s0.specific_energy(MU)) < 1e-9);
        let back = propagate(fwd, MU, -20_000.0);
        assert!((back.r - s0.r).length() < 1e-2);
    }

    #[test]
    fn analytic_matches_numeric_integration() {
        let s0 = ellipse();
        let mut num = s0;
        let dt = 0.5;
        let steps = 4_000;
        for _ in 0..steps {
            num = rk4_step(num, dt, |st| point_mass_gravity(st.r, MU));
        }
        let ana = propagate(s0, MU, dt * steps as f64);
        assert!((num.r - ana.r).length() < 1.0, "Abweichung {} m", (num.r - ana.r).length());
    }

    #[test]
    fn kdk_without_perturbation_is_exactly_kepler() {
        use crate::integrate::kepler_kdk;
        let s0 = ellipse();
        let a = kepler_kdk(s0, MU, 123.4, |_| Vec3::ZERO);
        let b = propagate(s0, MU, 123.4);
        assert_eq!(a, b); // bitgleich
    }

    #[test]
    fn kdk_with_constant_thrust_matches_rk4() {
        use crate::integrate::kepler_kdk;
        let s0 = ellipse();
        let thrust = Vec3::<Inertial>::new(0.0, 15.0, 3.0);
        let (mut a, mut b) = (s0, s0);
        for _ in 0..3_000 {
            a = kepler_kdk(a, MU, 0.02, |_| thrust);
            b = rk4_step(b, 0.02, |st| point_mass_gravity(st.r, MU) + thrust);
        }
        assert!((a.r - b.r).length() < 0.05, "{} m", (a.r - b.r).length());
    }

    #[test]
    fn finds_atmosphere_entry() {
        // Suborbitaler Bogen: steigt erst, fällt dann unter 670 km
        let r = 680_000.0;
        let s0 = StateVector::<Inertial>::new(Vec3::new(r, 0.0, 0.0), Vec3::new(800.0, 1_500.0, 0.0));
        let t = time_until_below(s0, MU, 670_000.0, 1.0e5).unwrap();
        assert!((propagate(s0, MU, t).r.length() - 670_000.0).abs() < 1e-3);
        // davor nie darunter
        for k in 0..100 {
            assert!(propagate(s0, MU, t * k as f64 / 100.0).r.length() >= 670_000.0 - 1e-3);
        }
        // Zu kurzes Zeitfenster: kein Ereignis
        assert!(time_until_below(s0, MU, 670_000.0, t * 0.5).is_none());
        // Radialer Fall aus der Ruhe (entartete Bahn)
        let drop = StateVector::<Inertial>::new(Vec3::new(800_000.0, 0.0, 0.0), Vec3::ZERO);
        let t = time_until_below(drop, MU, 670_000.0, 1.0e6).unwrap();
        assert!((propagate(drop, MU, t).r.length() - 670_000.0).abs() < 1e-3);
        // Stabile Bahn darüber: nie
        let circ = StateVector::<Inertial>::new(Vec3::new(r, 0.0, 0.0), Vec3::new(0.0, (MU / r).sqrt(), 0.0));
        assert!(time_until_below(circ, MU, 670_000.0, 1.0e7).is_none());
    }

    #[test]
    fn radial_bound_trajectory_has_apoapsis() {
        // Ruhe über der Oberfläche: e = 1, aber gebunden
        let s0 = StateVector::<Inertial>::new(Vec3::new(600_010.0, 0.0, 0.0), Vec3::ZERO);
        let el = elements(s0, MU);
        assert!((el.apoapsis.unwrap() - 600_010.0).abs() < 1e-3);
        assert!(el.time_to_apoapsis.unwrap().abs() < 1e-6);
    }

    #[test]
    fn apoapsis_timing() {
        let s0 = ellipse();
        let el = elements(s0, MU);
        let at_apo = propagate(s0, MU, el.time_to_apoapsis.unwrap());
        assert!(rel_err(at_apo.r.length(), el.apoapsis.unwrap()) < 1e-9);
        // am Apoapsis: radiale Geschwindigkeit null
        assert!(at_apo.r.normalize_or_zero().dot(at_apo.v).abs() < 1e-6);
    }
}
