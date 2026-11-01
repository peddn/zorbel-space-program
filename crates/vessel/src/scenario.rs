//! Standard-Szenario: Universum, Heimatplanet Zorb, Erstflug-Rakete, Startplatz.
//! Gemeinsam genutzt von Testflug, Spiel und Screenshot-Modus.

use crate::dynamics::RigidState;
use crate::mass::mass_properties;
use crate::part::PartLibrary;
use crate::vessel::{Vessel, VesselBuilder};
use sim_core::{CelestialBody, GeoPosition, Universe};

/// Das mitgelieferte Universum (Zol mit Zun, Zorb, Muz; Nachbarsystem).
pub fn universe() -> Universe {
    Universe::builtin()
}

/// Der Heimatplanet Zorb (600 km Radius, Atmosphäre bis 70 km).
pub fn home_planet(u: &Universe) -> &CelestialBody {
    u.body(u.find("Zorb").expect("Zorb fehlt im Universum"))
}

/// Zorb als eigenständiger Körper, für Tests der Physik auf einem Körper.
pub fn test_planet() -> CelestialBody {
    home_planet(&universe()).clone()
}

/// Zweistufige Rakete aus den fünf Basisteilen.
pub fn first_rocket(lib: &PartLibrary) -> Vessel {
    first_rocket_with(lib, true)
}

/// Dieselbe Rakete, wahlweise ohne Finnen (zum Ausprobieren der Stabilität).
pub fn first_rocket_with(lib: &PartLibrary, fins: bool) -> Vessel {
    let name = if fins { "Erstflug" } else { "Erstflug ohne Finnen" };
    let mut b = VesselBuilder::new(name, lib.part("kapsel"));
    let t2 = b.stack_below(0, lib.part("tank"));
    let e2 = b.stack_below(t2, lib.part("triebwerk"));
    let dec = b.stack_below(e2, lib.part("entkoppler"));
    let t1a = b.stack_below(dec, lib.part("tank"));
    let t1b = b.stack_below(t1a, lib.part("tank"));
    let t1c = b.stack_below(t1b, lib.part("tank"));
    let e1 = b.stack_below(t1c, lib.part("triebwerk"));
    if fins {
        b.radial_symmetric(t1c, lib.part("steuerflaeche"), 4, -0.5);
    }
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
