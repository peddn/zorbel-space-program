//! Physikgeometrie eines Schiffs: das einfache Polygonmodell aller Teile.
//!
//! Jedes Teil beschreibt in seiner Teiledatei eine einfache Hülle
//! ([`crate::part::Shape`]). Daraus entstehen Flächenelemente (Polygone mit
//! Normale und Fläche). Alles, was von der Anordnung der Teile abhängt,
//! wird hier einmal pro Schiff berechnet und zwischengespeichert:
//!
//! - **Verdeckung**: Für 18 Blickrichtungen (längs vorn/hinten, 16 quer)
//!   wird das Schiff mit einem Tiefenpuffer gerastert. Jedes Flächenelement
//!   weiß dann, welcher Anteil von ihm aus dieser Richtung frei liegt.
//!   Verbindungsflächen zwischen gestapelten Teilen sind so automatisch
//!   verdeckt, ebenso alles hinter der Spitze beim senkrechten Flug.
//! - **Umströmte Fläche**: Anteil, der aus irgendeiner Richtung frei liegt
//!   (für Reibung; innere Flächen zählen nicht).
//! - **Spannweitenverteilung** längs der Schiffsachse für die
//!   Schlankkörpertheorie (Auftrieb dort, wo das Schiff breiter wird).
//! - **Kontaktpunkte** für den Boden: die Eckpunkte der Polygone.
//!
//! Die Darstellung benutzt eigene, detaillierte Modelle; dieses Polygonmodell
//! dient nur der Physik.

use crate::part::Shape;
use crate::vessel::Vessel;
use glam::DVec3;
use sim_core::{PartLocal, Vec3, VesselFrame};
use std::f64::consts::TAU;

/// Seitenzahl, mit der runde Hüllen als Vieleck angenähert werden.
pub const HULL_SIDES: usize = 12;
/// Anzahl der Querrichtungen für die Verdeckung.
pub const CROSS_DIRECTIONS: usize = 16;
/// Auflösung des Tiefenpuffers (Zellen über die größte Ausdehnung).
const GRID: usize = 160;

impl Shape {
    /// Polygone der Hülle im Teilsystem, Eckpunkte gegen den Uhrzeigersinn
    /// von außen gesehen (Normale zeigt nach außen).
    pub fn polygons(&self) -> Vec<Vec<Vec3<PartLocal>>> {
        match *self {
            Shape::Cylinder { radius, length } => frustum(radius, radius, length),
            Shape::Cone { radius_bottom, radius_top, length } => frustum(radius_bottom, radius_top, length),
            Shape::Plate { span, chord, thickness } => cuboid(span, thickness, chord),
        }
    }
}

fn ring(r: f64, z: f64) -> Vec<Vec3<PartLocal>> {
    (0..HULL_SIDES)
        .map(|k| {
            let a = TAU * k as f64 / HULL_SIDES as f64;
            Vec3::new(r * a.cos(), r * a.sin(), z)
        })
        .collect()
}

fn frustum(rb: f64, rt: f64, length: f64) -> Vec<Vec<Vec3<PartLocal>>> {
    let (zb, zt) = (-length / 2.0, length / 2.0);
    let (bottom, top) = (ring(rb, zb), ring(rt, zt));
    let mut polys = Vec::new();
    for k in 0..HULL_SIDES {
        let j = (k + 1) % HULL_SIDES;
        if rt > 1e-9 {
            polys.push(vec![bottom[k], bottom[j], top[j], top[k]]);
        } else {
            polys.push(vec![bottom[k], bottom[j], Vec3::new(0.0, 0.0, zt)]);
        }
    }
    if rt > 1e-9 {
        polys.push(top.clone());
    }
    if rb > 1e-9 {
        polys.push(bottom.into_iter().rev().collect());
    }
    polys
}

fn cuboid(x: f64, y: f64, z: f64) -> Vec<Vec<Vec3<PartLocal>>> {
    let (hx, hy, hz) = (x / 2.0, y / 2.0, z / 2.0);
    let p = |sx: f64, sy: f64, sz: f64| Vec3::new(sx * hx, sy * hy, sz * hz);
    vec![
        vec![p(1., -1., -1.), p(1., 1., -1.), p(1., 1., 1.), p(1., -1., 1.)],
        vec![p(-1., -1., -1.), p(-1., -1., 1.), p(-1., 1., 1.), p(-1., 1., -1.)],
        vec![p(-1., 1., -1.), p(-1., 1., 1.), p(1., 1., 1.), p(1., 1., -1.)],
        vec![p(-1., -1., -1.), p(1., -1., -1.), p(1., -1., 1.), p(-1., -1., 1.)],
        vec![p(-1., -1., 1.), p(1., -1., 1.), p(1., 1., 1.), p(-1., 1., 1.)],
        vec![p(-1., -1., -1.), p(-1., 1., -1.), p(1., 1., -1.), p(1., -1., -1.)],
    ]
}

/// Ein ebenes Flächenelement im Schiffssystem.
#[derive(Debug, Clone)]
pub struct Panel {
    pub part: usize,
    pub verts: Vec<DVec3>,
    /// Nach außen zeigende Einheitsnormale.
    pub normal: DVec3,
    pub area: f64,
    pub centroid: DVec3,
}

/// Auftriebsbeitrag eines Teils nach der Schlankkörpertheorie für eine
/// Anströmseite und eine Querrichtung.
#[derive(Debug, Clone, Copy, Default)]
pub struct LiftShare {
    /// Zuwachs der "scheinbaren Fläche" π(s/2)², den das Teil verursacht (m²).
    pub delta_area: f64,
    /// Lage des Beitrags auf der Schiffsachse (z im Schiffssystem).
    pub center_z: f64,
}

#[derive(Debug, Clone)]
pub struct HullGeometry {
    pub panels: Vec<Panel>,
    /// Blickrichtungen (Einheitsvektoren, Schiffssystem): 0 = von +Z,
    /// 1 = von −Z, danach die Querrichtungen.
    pub directions: Vec<DVec3>,
    /// `visible[d][j]`: frei liegender Anteil von Panel j aus Richtung d.
    pub visible: Vec<Vec<f64>>,
    /// Umströmter Anteil jedes Panels.
    pub wetted: Vec<f64>,
    /// `lift[e][a][part]`: e = 0 Spitze voran (+Z angeströmt), 1 Heck voran;
    /// a = 0 Querströmung entlang X, 1 entlang Y.
    pub lift: [[Vec<LiftShare>; 2]; 2],
    /// Seitenverhältnis-Korrektur pro Teil (1 für schlanke Körper).
    pub aspect_correction: Vec<f64>,
    /// Kontaktpunkte für den Boden: (Teil, Punkt im Schiffssystem).
    pub contact_points: Vec<(usize, Vec3<VesselFrame>)>,
}

impl HullGeometry {
    pub fn build(v: &Vessel) -> Self {
        let mut panels = Vec::new();
        let mut contact_points = Vec::new();
        for (i, p) in v.parts.iter().enumerate() {
            let to_vessel = |q: Vec3<PartLocal>| (p.origin + p.orientation.apply(q)).raw();
            let polys = p.def.hull.polygons();
            let mut seen: Vec<DVec3> = Vec::new();
            for poly in polys {
                let verts: Vec<DVec3> = poly.into_iter().map(to_vessel).collect();
                for &q in &verts {
                    if !seen.iter().any(|s| s.distance_squared(q) < 1e-12) {
                        seen.push(q);
                    }
                }
                if let Some(panel) = make_panel(i, verts) {
                    panels.push(panel);
                }
            }
            contact_points.extend(seen.into_iter().map(|q| (i, Vec3::from_raw(q))));
        }

        let mut directions = vec![DVec3::Z, -DVec3::Z];
        for k in 0..CROSS_DIRECTIONS {
            let a = TAU * k as f64 / CROSS_DIRECTIONS as f64;
            directions.push(DVec3::new(a.cos(), a.sin(), 0.0));
        }
        let visible: Vec<Vec<f64>> = directions.iter().map(|&d| visibility(&panels, d)).collect();
        let wetted = (0..panels.len()).map(|j| visible.iter().map(|v| v[j]).fold(0.0, f64::max)).collect();
        let (lift, aspect_correction) = slender_shares(v, &panels);

        Self { panels, directions, visible, wetted, lift, aspect_correction, contact_points }
    }

    /// Index der Querrichtung, die `dir` (senkrecht zur Achse) am nächsten liegt.
    pub fn nearest_cross_direction(&self, dir: DVec3) -> usize {
        let a = dir.y.atan2(dir.x).rem_euclid(TAU);
        2 + ((a / TAU * CROSS_DIRECTIONS as f64).round() as usize) % CROSS_DIRECTIONS
    }

    /// Größter Abstand eines Kontaktpunkts von `center`.
    pub fn extent_from(&self, center: Vec3<VesselFrame>) -> f64 {
        self.contact_points.iter().map(|(_, p)| (*p - center).length()).fold(0.0, f64::max)
    }
}

fn make_panel(part: usize, verts: Vec<DVec3>) -> Option<Panel> {
    // Newell-Verfahren: robuste Normale und Fläche für ebene Vielecke
    let mut n = DVec3::ZERO;
    for k in 0..verts.len() {
        let (a, b) = (verts[k], verts[(k + 1) % verts.len()]);
        n += DVec3::new((a.y - b.y) * (a.z + b.z), (a.z - b.z) * (a.x + b.x), (a.x - b.x) * (a.y + b.y));
    }
    let len = n.length();
    if len < 1e-12 {
        return None;
    }
    let centroid = verts.iter().copied().sum::<DVec3>() / verts.len() as f64;
    Some(Panel { part, normal: n / len, area: len / 2.0, centroid, verts })
}

/// Frei liegender Anteil jedes Panels, gesehen aus Richtung `d` (Kamera
/// weit draußen bei +d, Blick nach −d).
fn visibility(panels: &[Panel], d: DVec3) -> Vec<f64> {
    let e1 = d.any_orthonormal_vector();
    let e2 = d.cross(e1);
    let proj = |p: DVec3| (p.dot(e1), p.dot(e2), -p.dot(d));

    let (mut umin, mut umax, mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for p in panels.iter().flat_map(|p| &p.verts) {
        let (u, v, _) = proj(*p);
        umin = umin.min(u);
        umax = umax.max(u);
        vmin = vmin.min(v);
        vmax = vmax.max(v);
    }
    let size = (umax - umin).max(vmax - vmin).max(1e-6);
    let cell = size / GRID as f64;
    let (nu, nv) = (((umax - umin) / cell) as usize + 2, ((vmax - vmin) / cell) as usize + 2);
    let index = |u: f64, v: f64| {
        let i = (((u - umin) / cell) as isize).clamp(0, nu as isize - 1) as usize;
        let j = (((v - vmin) / cell) as isize).clamp(0, nv as isize - 1) as usize;
        j * nu + i
    };
    let mut depth = vec![f64::INFINITY; nu * nv];
    let facing = |p: &Panel| p.normal.dot(d) > 1e-9;

    // Tiefenpuffer: nur zur Kamera gewandte Flächen (geschlossene Hüllen)
    for p in panels.iter().filter(|p| facing(p)) {
        let pts: Vec<(f64, f64, f64)> = p.verts.iter().map(|q| proj(*q)).collect();
        for k in 1..pts.len() - 1 {
            raster_triangle(pts[0], pts[k], pts[k + 1], umin, vmin, cell, nu, nv, &mut depth);
            // Sehr kleine Dreiecke treffen evtl. keinen Zellmittelpunkt
            let c = ((pts[0].0 + pts[k].0 + pts[k + 1].0) / 3.0, (pts[0].1 + pts[k].1 + pts[k + 1].1) / 3.0);
            let cz = (pts[0].2 + pts[k].2 + pts[k + 1].2) / 3.0;
            let idx = index(c.0, c.1);
            depth[idx] = depth[idx].min(cz);
        }
    }

    let bias = 1e-3 + size * 2e-3;
    panels
        .iter()
        .map(|p| {
            if !facing(p) {
                return 0.0;
            }
            // Stichproben: Mittelpunkt und halbe Wege zu den Ecken
            let samples = std::iter::once(p.centroid).chain(p.verts.iter().map(|q| (p.centroid + *q) / 2.0));
            // Schräge Flächen ändern ihre Tiefe innerhalb einer Zelle stark:
            // Toleranz mit der Neigung skalieren (wie bei Schattenkarten).
            let c = p.normal.dot(d);
            let slope = ((1.0 - c * c).max(0.0).sqrt() / c).min(10.0);
            let tol = bias + 1.5 * cell * slope;
            let (mut hit, mut total) = (0, 0);
            for s in samples {
                let (u, v, z) = proj(s);
                total += 1;
                if z <= depth[index(u, v)] + tol {
                    hit += 1;
                }
            }
            hit as f64 / total as f64
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn raster_triangle(
    a: (f64, f64, f64),
    b: (f64, f64, f64),
    c: (f64, f64, f64),
    umin: f64,
    vmin: f64,
    cell: f64,
    nu: usize,
    nv: usize,
    depth: &mut [f64],
) {
    let area = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    if area.abs() < 1e-18 {
        return;
    }
    let lo = |x: f64, m: f64| (((x - m) / cell - 0.5).floor().max(0.0)) as usize;
    let hi = |x: f64, m: f64, n: usize| ((((x - m) / cell - 0.5).ceil()).max(0.0) as usize).min(n - 1);
    let (i0, i1) = (lo(a.0.min(b.0).min(c.0), umin), hi(a.0.max(b.0).max(c.0), umin, nu));
    let (j0, j1) = (lo(a.1.min(b.1).min(c.1), vmin), hi(a.1.max(b.1).max(c.1), vmin, nv));
    for j in j0..=j1 {
        for i in i0..=i1 {
            let (u, v) = (umin + (i as f64 + 0.5) * cell, vmin + (j as f64 + 0.5) * cell);
            let w0 = ((b.0 - u) * (c.1 - v) - (b.1 - v) * (c.0 - u)) / area;
            let w1 = ((c.0 - u) * (a.1 - v) - (c.1 - v) * (a.0 - u)) / area;
            let w2 = 1.0 - w0 - w1;
            if w0 >= -1e-9 && w1 >= -1e-9 && w2 >= -1e-9 {
                let z = w0 * a.2 + w1 * b.2 + w2 * c.2;
                let k = j * nu + i;
                if z < depth[k] {
                    depth[k] = z;
                }
            }
        }
    }
}

/// Seitenverhältnis-Korrektur nach Helmbold relativ zur Schlankkörpertheorie:
/// 1 für sehr schlanke Formen, kleiner für breite, kurze Flächen.
pub fn aspect_correction(aspect: f64) -> f64 {
    4.0 / (2.0 + (4.0 + aspect * aspect).sqrt())
}

/// Scheinbare Fläche eines Querschnitts aus Rumpf (Radius `a`) und dünnen
/// Flächen bis zur Halbspannweite `s` (exakt nach der Schlankkörpertheorie).
/// Reiner Rumpf (s = a): π a². Reine Fläche (a = 0): π s².
pub fn apparent_area(s: f64, a: f64) -> f64 {
    let s = s.max(a);
    if s <= 0.0 {
        return 0.0;
    }
    let k = (a / s).powi(2);
    std::f64::consts::PI * s * s * (1.0 - k + k * k)
}

/// Spannweitenverteilung längs der Achse → Auftriebsbeiträge pro Teil.
#[allow(clippy::type_complexity)]
fn slender_shares(v: &Vessel, panels: &[Panel]) -> ([[Vec<LiftShare>; 2]; 2], Vec<f64>) {
    let n_parts = v.parts.len();
    // Punkte auf allen Kanten (Abstand ≤ 2 cm)
    let mut pts: Vec<(usize, DVec3)> = Vec::new();
    for p in panels {
        for k in 0..p.verts.len() {
            let (a, b) = (p.verts[k], p.verts[(k + 1) % p.verts.len()]);
            let n = ((a.distance(b) / 0.02).ceil() as usize).max(1);
            pts.extend((0..n).map(|s| (p.part, a.lerp(b, s as f64 / n as f64))));
        }
    }
    let empty = || [[vec![LiftShare::default(); n_parts], vec![LiftShare::default(); n_parts]], [vec![LiftShare::default(); n_parts], vec![LiftShare::default(); n_parts]]];
    if pts.is_empty() {
        return (empty(), vec![1.0; n_parts]);
    }
    // Punkte in Strömungsrichtung durchlaufen. Wer die Breite (quer zur
    // Querströmung) über das bisherige Maximum hinaus vergrößert, bekommt den
    // Zuwachs der scheinbaren Fläche π(s/2)² gutgeschrieben, an seiner Stelle.
    let on_axis: Vec<bool> = v.parts.iter().map(|p| p.origin.x().hypot(p.origin.y()) < 1e-6).collect();
    let mut lift = empty();
    for (end, per_end) in lift.iter_mut().enumerate() {
        let mut order: Vec<usize> = (0..pts.len()).collect();
        // Auf µm gerundet sortieren; bei gleicher Lage gewinnt das Teil mit
        // dem kleineren Index (näher an der Wurzel), damit Stoßstellen nicht
        // durch Rundung dem falschen Teil zugeschlagen werden.
        let key = |k: usize| (pts[k].1.z * 1e6).round() as i64;
        order.sort_by(|&x, &y| {
            let by_z = if end == 0 { key(y).cmp(&key(x)) } else { key(x).cmp(&key(y)) };
            by_z.then(pts[x].0.cmp(&pts[y].0))
        });
        for (axis, slot) in per_end.iter_mut().enumerate() {
            // Strömung entlang X → Breite entlang Y und umgekehrt
            let coord = |q: DVec3| if axis == 0 { q.y } else { q.x };
            let (mut lo, mut hi) = (coord(pts[order[0]].1), coord(pts[order[0]].1));
            let mut body: f64 = 0.0;
            let mut area = 0.0;
            let mut acc = vec![(0.0, 0.0); n_parts];
            for &k in &order {
                let (part, q) = pts[k];
                let c = coord(q);
                // Rumpfradius: Teile auf der Achse
                let r = if on_axis[part] { q.x.hypot(q.y) } else { 0.0 };
                if c < lo || c > hi || r > body {
                    lo = lo.min(c);
                    hi = hi.max(c);
                    body = body.max(r);
                    let new_area = apparent_area((hi - lo) / 2.0, body).max(area);
                    let gain = new_area - area;
                    if gain > 0.0 {
                        area = new_area;
                        acc[part].0 += gain;
                        acc[part].1 += gain * q.z;
                    }
                }
            }
            *slot = acc
                .iter()
                .map(|&(ds, mz)| LiftShare { delta_area: ds, center_z: if ds > 0.0 { mz / ds } else { 0.0 } })
                .collect();
        }
    }

    // Seitenverhältnis pro Teil: 2 · seitliche Ausdehnung / Länge
    let correction = v
        .parts
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let mine: Vec<DVec3> = pts.iter().filter(|q| q.0 == i).map(|q| q.1).collect();
            if mine.is_empty() {
                return 1.0;
            }
            let r = |q: &DVec3| q.x.hypot(q.y);
            let rmax = mine.iter().map(r).fold(0.0, f64::max);
            let rmin = if on_axis[i] { 0.0 } else { mine.iter().map(r).fold(f64::INFINITY, f64::min) };
            let len = mine.iter().map(|q| q.z).fold(f64::NEG_INFINITY, f64::max) - mine.iter().map(|q| q.z).fold(f64::INFINITY, f64::min);
            aspect_correction(2.0 * (rmax - rmin) / len.max(1e-6))
        })
        .collect();
    (lift, correction)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::PartLibrary;
    use crate::scenario::first_rocket;
    use crate::vessel::VesselBuilder;

    #[test]
    fn hull_polygons_are_closed_and_outward() {
        for shape in [
            Shape::Cylinder { radius: 0.6, length: 2.0 },
            Shape::Cone { radius_bottom: 0.6, radius_top: 0.25, length: 1.2 },
            Shape::Cone { radius_bottom: 0.6, radius_top: 0.0, length: 1.2 },
            Shape::Plate { span: 0.8, chord: 1.0, thickness: 0.05 },
        ] {
            let panels: Vec<Panel> = shape.polygons().into_iter().filter_map(|p| make_panel(0, p.into_iter().map(|q| q.raw()).collect())).collect();
            // Geschlossen: Summe der Flächenvektoren ist null
            let sum: DVec3 = panels.iter().map(|p| p.normal * p.area).sum();
            assert!(sum.length() < 1e-9, "{shape:?}: {sum:?}");
            // Nach außen: Normale zeigt vom Mittelpunkt weg
            for p in &panels {
                assert!(p.normal.dot(p.centroid) > 0.0, "{shape:?}");
            }
        }
    }

    #[test]
    fn stacked_tanks_hide_the_joint() {
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("tank"));
        b.stack_below(0, lib.part("tank"));
        let v = b.build();
        let g = v.geometry();
        let caps: Vec<(usize, &Panel)> = g.panels.iter().enumerate().filter(|(_, p)| p.normal.z.abs() > 0.99).collect();
        for (j, p) in caps {
            let outer = (p.normal.z > 0.0 && p.part == 0) || (p.normal.z < 0.0 && p.part == 1);
            let dir = if p.normal.z > 0.0 { 0 } else { 1 };
            if outer {
                assert!(g.visible[dir][j] > 0.99, "Stirnfläche außen muss frei sein");
                assert!(g.wetted[j] > 0.99);
            } else {
                assert_eq!(g.visible[dir][j], 0.0, "Fuge muss verdeckt sein");
                assert_eq!(g.wetted[j], 0.0, "Fuge ist nicht umströmt");
            }
        }
    }

    #[test]
    fn vertical_flight_sees_only_nose_and_fin_edges() {
        // Dein Beispiel: senkrecht nach oben angeströmt sind nur die Spitze
        // (Kapsel) und die Vorderkanten der Finnen frei.
        let v = first_rocket(&PartLibrary::builtin());
        let g = v.geometry();
        let mut by_part = vec![0.0; v.parts.len()];
        for (j, p) in g.panels.iter().enumerate() {
            by_part[p.part] += g.visible[0][j] * p.area * p.normal.z.max(0.0);
        }
        let capsule = by_part[0];
        let r: f64 = 0.625;
        let ring_area = std::f64::consts::PI * r * r * (HULL_SIDES as f64 / TAU * (TAU / HULL_SIDES as f64).sin());
        assert!((capsule - ring_area).abs() / ring_area < 0.02, "Kapsel {capsule} m², erwartet {ring_area} m²");
        for (i, p) in v.parts.iter().enumerate().skip(1) {
            if p.def.control_surface().is_some() {
                let edge = 0.8 * 0.05;
                assert!((by_part[i] - edge).abs() / edge < 0.05, "Finne {i}: {} m²", by_part[i]);
            } else {
                // Tank-/Triebwerksringe, die über die Engstelle am Triebwerk
                // hinausragen, liegen frei (Ringflächen), sonst nichts
                assert!(by_part[i] < 0.75, "Teil {} ({}) frei: {} m²", i, p.def.name, by_part[i]);
            }
        }
    }

    #[test]
    fn body_of_revolution_gives_munk_area() {
        // Rumpf ohne Finnen: Spannweitenzuwachs von vorn = Querschnitt π r²
        let lib = PartLibrary::builtin();
        let mut b = VesselBuilder::new("Test", lib.part("kapsel"));
        b.stack_below(0, lib.part("tank"));
        let v = b.build();
        let g = v.geometry();
        let total: f64 = g.lift[0][0].iter().map(|s| s.delta_area).sum();
        let r = 0.625_f64;
        let expect = std::f64::consts::PI * r * r;
        assert!((total - expect).abs() / expect < 0.02, "{total} vs {expect}");
        // Der Zuwachs entsteht an der Kapsel, nicht am Tank
        assert!(g.lift[0][0][1].delta_area < 0.01 * total, "Tank: {}", g.lift[0][0][1].delta_area);
    }
}
