//! Standard-Szenario: Testplanet, Erstflug-Rakete, Startplatz.
//! Gemeinsam genutzt von Testflug, Spiel und Screenshot-Modus.

use crate::dynamics::RigidState;
use crate::mass::mass_properties;
use crate::part::PartLibrary;
use crate::vessel::{Vessel, VesselBuilder};
use sim_core::{Atmosphere, CelestialBody, GeoPosition};

/// Planet mit KSP-artigen Kenndaten: 600 km Radius, μ = 3,5316·10¹² m³/s²,
/// Atmosphäre bis 70 km, siderischer Tag 21 549 s.
pub fn test_planet() -> CelestialBody {
    CelestialBody {
        name: "Testplanet".into(),
        mu: 3.5316e12,
        radius: 600_000.0,
        soi_radius: 8.4e7,
        atmosphere: Atmosphere::earthlike(),
        rotation_rate: CelestialBody::rate_from_period(21_549.425),
        rotation_offset: 0.0,
    }
}

/// Zweistufige Rakete aus den fünf Basisteilen.
pub fn first_rocket(lib: &PartLibrary) -> Vessel {
    let mut b = VesselBuilder::new("Erstflug", lib.part("kapsel"));
    let t2 = b.stack_below(0, lib.part("tank"));
    let e2 = b.stack_below(t2, lib.part("triebwerk"));
    let dec = b.stack_below(e2, lib.part("entkoppler"));
    let t1a = b.stack_below(dec, lib.part("tank"));
    let t1b = b.stack_below(t1a, lib.part("tank"));
    let t1c = b.stack_below(t1b, lib.part("tank"));
    let e1 = b.stack_below(t1c, lib.part("triebwerk"));
    b.radial_symmetric(t1c, lib.part("steuerflaeche"), 4, -0.5);
    // Höchste Stufe zuerst: 1 = Zündung unten, 0 = Trennen + Zündung oben
    b.set_stage(e1, 1).set_stage(dec, 0).set_stage(e2, 0);
    b.build()
}

/// Startplatz am Äquator, Länge 0.
pub const LAUNCH_SITE: GeoPosition = GeoPosition { latitude: 0.0, longitude: 0.0, altitude: 0.0 };

/// Höhe des Schwerpunkts über der Unterkante des Schiffs.
pub fn com_height(v: &Vessel) -> f64 {
    let mp = mass_properties(v);
    let bottom = v.parts.iter().map(|p| p.origin.z() - p.def.hull.stack_length() / 2.0).fold(f64::INFINITY, f64::min);
    mp.com.z() - bottom
}

/// Schiff steht auf dem Startplatz (dreht mit dem Planeten mit).
pub fn on_launch_site(body: &CelestialBody, v: &Vessel, time: f64) -> RigidState {
    RigidState::at_rest(body, GeoPosition { altitude: com_height(v) + 0.1, ..LAUNCH_SITE }, time)
}
