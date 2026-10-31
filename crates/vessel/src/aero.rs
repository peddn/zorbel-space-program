//! Aerodynamik pro Bauteil.
//!
//! Jedes Teil bekommt seine eigene Luftkraft am eigenen Ort. Damit ergeben
//! sich Stabilität (Finnen hinten), Dämpfung (Drehung erzeugt Anströmung)
//! und Strukturlasten (Biegung bei Anstellwinkel) ohne Sonderregeln: die
//! Kräfte gehen in die Starrkörperbewegung und in die Lastberechnung ein.
//!
//! Alle Kräfte sind proportional zur Luftdichte. Im Vakuum (Dichte 0) sind
//! sie exakt null, ohne dass irgendwo ein Fall unterschieden wird.
//!
//! Modell (bewusst einfach, aber stetig in Richtung und Machzahl):
//! - Stirnflächen: angeströmt mit `CD_FLAT`, im Lee mit `CD_BASE`
//!   (Heckwiderstand). Von Nachbarteilen verdeckte Flächen zählen nicht.
//! - Kegelmantel: Newtonscher Widerstand 2·sin²θ auf die Ringfläche.
//! - Mantel quer angeströmt: `CD_SIDE` auf die Projektionsfläche.
//! - Platten (Finnen, Steuerflächen): Normalkraft aus Anstellwinkel
//!   (lineare Steigung + Querströmungsanteil), Reibung tangential.
//! - Machzahl: glatter Anstieg des Widerstands im Transschall.

use crate::part::Shape;
use crate::vessel::{PartInstance, Vessel};
use sim_core::{AtmosphereSample, PartLocal, Vec3};
use std::f64::consts::PI;

pub const CD_FLAT: f64 = 0.8;
pub const CD_BASE: f64 = 0.12;
pub const CD_SIDE: f64 = 1.0;
pub const CD_FRICTION: f64 = 0.01;
/// Normalkraftbeiwert einer quer angeströmten Platte.
pub const CN_CROSSFLOW: f64 = 1.8;

/// Aerodynamisch wirksame Flächen eines Teils, inklusive Verdeckung.
#[derive(Debug, Clone, Copy, Default)]
pub struct PartAero {
    /// Freie ebene Fläche der +Z-Stirnseite.
    pub top_flat: f64,
    /// Freie ebene Fläche der −Z-Stirnseite.
    pub bottom_flat: f64,
    /// Widerstandsfläche des Kegelmantels bei Anströmung von +Z.
    pub top_slope_cda: f64,
    /// … von −Z (Kegel, der unten schmaler ist).
    pub bottom_slope_cda: f64,
    /// Widerstandsfläche bei Queranströmung.
    pub side_cda: f64,
}

fn disc(r: f64) -> f64 {
    PI * r * r
}

/// Ringfläche zwischen `r_inner` und `r_outer` (0, wenn verdeckt).
fn annulus(r_outer: f64, r_inner: f64) -> f64 {
    (disc(r_outer) - disc(r_inner)).max(0.0)
}

/// Flächen aller Teile. Die Verdeckung folgt aus dem Verbindungsgraphen,
/// deshalb ändert sie sich automatisch, wenn eine Stufe abgetrennt wird.
pub fn aero_geometry(v: &Vessel) -> Vec<PartAero> {
    let n = v.parts.len();
    // Radius des Nachbarn an der oberen bzw. unteren Stirnseite.
    let mut cover_top = vec![0.0_f64; n];
    let mut cover_bottom = vec![0.0_f64; n];
    for (i, p) in v.parts.iter().enumerate() {
        let (Some(par), Some(j)) = (p.parent, p.joint) else { continue };
        // Nur Stapelverbindungen entlang der Längsachse verdecken.
        if j.axis.z() > -0.999 {
            continue;
        }
        if let (Some((_, my_top)), Some((par_bottom, _))) =
            (p.def.shape.end_radii(), v.parts[par].def.shape.end_radii())
        {
            cover_top[i] = cover_top[i].max(par_bottom);
            cover_bottom[par] = cover_bottom[par].max(my_top);
        }
    }

    v.parts
        .iter()
        .enumerate()
        .map(|(i, p)| match p.def.shape {
            Shape::Cylinder { radius, length } => PartAero {
                top_flat: annulus(radius, cover_top[i]),
                bottom_flat: annulus(radius, cover_bottom[i]),
                side_cda: CD_SIDE * 2.0 * radius * length,
                ..Default::default()
            },
            Shape::Cone { radius_bottom: rb, radius_top: rt, length } => {
                let sin_t = (rb - rt).abs() / ((rb - rt).powi(2) + length * length).sqrt();
                let newton = 2.0 * sin_t * sin_t;
                let (top_slope_cda, bottom_slope_cda) = if rb >= rt {
                    (newton * annulus(rb, rt.max(cover_top[i])), 0.0)
                } else {
                    (0.0, newton * annulus(rt, rb.max(cover_bottom[i])))
                };
                PartAero {
                    top_flat: annulus(rt, cover_top[i]),
                    bottom_flat: annulus(rb, cover_bottom[i]),
                    top_slope_cda,
                    bottom_slope_cda,
                    side_cda: CD_SIDE * (rb + rt) * length,
                }
            }
            Shape::Plate { .. } => PartAero::default(),
        })
        .collect()
}

/// Widerstandsanstieg im Transschall, glatt in der Machzahl.
pub fn mach_factor(mach: f64) -> f64 {
    let bump = 0.7 * (-((mach - 1.1) / 0.3).powi(2)).exp();
    let t = ((mach - 0.8) / 0.6).clamp(0.0, 1.0);
    let supersonic = 0.3 * t * t * (3.0 - 2.0 * t);
    1.0 + bump + supersonic
}

/// Luftkraft auf ein Teil, im lokalen System des Teils.
/// `v_local`: Geschwindigkeit des Teils relativ zur Luft.
/// `deflection`: Ausschlag einer Steuerfläche in rad (sonst 0).
pub fn part_force(
    p: &PartInstance,
    geo: &PartAero,
    v_local: Vec3<PartLocal>,
    air: &AtmosphereSample,
    deflection: f64,
) -> Vec3<PartLocal> {
    let u2 = v_local.length_squared();
    if u2 == 0.0 || air.density == 0.0 {
        return Vec3::ZERO;
    }
    let u = u2.sqrt();
    let d = v_local / u;
    let q = 0.5 * air.density * u2;
    let mach = mach_factor(u / air.speed_of_sound);

    match p.def.shape {
        Shape::Plate { span, chord, thickness } => {
            let area = span * chord;
            let (sd, cd) = deflection.sin_cos();
            // Normale (Dickenrichtung), um die Spannweite (X) gedreht.
            let n = Vec3::<PartLocal>::new(0.0, cd, sd);
            let lift_slope = p.def.control_surface().map_or(2.0 * PI * 0.8, |c| c.lift_slope);
            let s = d.dot(n);
            let c = (1.0 - s * s).max(0.0).sqrt();
            let cn = lift_slope * s * c + CN_CROSSFLOW * s * s.abs();
            let normal = n * (-q * area * cn);
            let tangential = d.reject_from_unit(n);
            let edge = span * thickness * CD_FLAT * d.z().abs();
            let friction = tangential * (-q * (2.0 * area * CD_FRICTION + edge));
            (normal + friction) * mach
        }
        _ => {
            let dz = d.z();
            let axial_top = if dz > 0.0 { CD_FLAT } else { CD_BASE } * geo.top_flat;
            let axial_bottom = if dz < 0.0 { CD_FLAT } else { CD_BASE } * geo.bottom_flat;
            let cda = (axial_top + axial_bottom) * dz.abs()
                + geo.top_slope_cda * dz.max(0.0)
                + geo.bottom_slope_cda * (-dz).max(0.0)
                + geo.side_cda * (d.x() * d.x() + d.y() * d.y()).sqrt();
            d * (-q * cda * mach)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::PartLibrary;
    use crate::vessel::VesselBuilder;
    use sim_core::Atmosphere;

    #[test]
    fn stacking_hides_the_joined_faces() {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("tank"));
        b.stack_below(0, lib.part("tank"));
        let geo = aero_geometry(&b.build());
        let full = disc(0.625);
        assert!((geo[0].top_flat - full).abs() < 1e-12);
        assert_eq!(geo[0].bottom_flat, 0.0);
        assert_eq!(geo[1].top_flat, 0.0);
        assert!((geo[1].bottom_flat - full).abs() < 1e-12);
    }

    #[test]
    fn capsule_nose_first_has_less_drag_than_tail_first() {
        let lib = PartLibrary::builtin();
        let v = VesselBuilder::new("Test", lib.part("kapsel")).build();
        let geo = aero_geometry(&v);
        let air = Atmosphere::earthlike().sample(0.0);
        let p = &v.parts[0];
        let nose = part_force(p, &geo[0], Vec3::new(0.0, 0.0, 100.0), &air, 0.0).length();
        let tail = part_force(p, &geo[0], Vec3::new(0.0, 0.0, -100.0), &air, 0.0).length();
        assert!(nose < tail, "Spitze voran {nose} N, Boden voran {tail} N");
    }

    #[test]
    fn forces_vanish_in_vacuum_and_scale_with_density() {
        let lib = PartLibrary::builtin();
        let v = VesselBuilder::new("Test", lib.part("steuerflaeche")).build();
        let geo = aero_geometry(&v);
        let vel = Vec3::new(10.0, 20.0, 300.0);
        let mut air = Atmosphere::earthlike().sample(10_000.0);
        let f1 = part_force(&v.parts[0], &geo[0], vel, &air, 0.1);
        air.density *= 1e-9;
        let f2 = part_force(&v.parts[0], &geo[0], vel, &air, 0.1);
        assert!((f2 * 1e9 - f1).length() < 1e-9 * f1.length());
        air.density = 0.0;
        assert_eq!(part_force(&v.parts[0], &geo[0], vel, &air, 0.1), Vec3::ZERO);
    }

    #[test]
    fn plate_force_is_continuous_through_zero_angle() {
        let lib = PartLibrary::builtin();
        let v = VesselBuilder::new("Test", lib.part("steuerflaeche")).build();
        let geo = aero_geometry(&v);
        let air = Atmosphere::earthlike().sample(0.0);
        let at = |a: f64| part_force(&v.parts[0], &geo[0], Vec3::new(0.0, 100.0 * a.sin(), 100.0 * a.cos()), &air, 0.0);
        let eps = 1e-7;
        let scale = at(0.1).length();
        assert!((at(eps) - at(-eps)).length() < 1e-5 * scale);
        // Auftrieb entgegen der Anströmseite
        assert!(at(0.1).y() < 0.0 && at(-0.1).y() > 0.0);
    }

    #[test]
    fn mach_factor_is_smooth() {
        let mut last = mach_factor(0.0);
        for k in 1..=5_000 {
            let m = k as f64 * 0.001;
            let f = mach_factor(m);
            assert!((f - last).abs() < 0.01, "Sprung bei M {m}");
            last = f;
        }
    }
}
