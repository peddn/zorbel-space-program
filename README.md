# Zorbel Space Program

Moderne Umsetzung der Kernmechaniken von Kerbal Space Program in Rust.
Stand: headless Simulationskern mit Atmosphäre, Teile, Schiffsgraph,
Aerodynamik pro Teil, Lastberechnung und Staging.

## Starten

```sh
cargo test --workspace                          # alle Tests
cargo run --bin zorbel                          # Headless-Testflug in eine 100-km-Kreisbahn
cargo run -p zorbel-game                        # Spiel im Fenster
cargo run -p zorbel-game -- --screenshots       # fester Flug ohne Fenster, Bilder nach screenshots/
```

Linux braucht für das Spiel die üblichen Bevy-Abhängigkeiten (Vulkan- oder
GL-Treiber, X11 oder Wayland). Ohne GPU reicht für den Screenshot-Modus
Mesas Software-Treiber (`mesa-vulkan-drivers`, Adapter "llvmpipe").

### Steuerung

| Taste | Wirkung |
|---|---|
| Leertaste | Nächste Stufe (die erste löst die Startklammern; ohne Gas bleibt die Rakete stehen) |
| Shift / Strg | Gas auf / ab, Z voll, X aus |
| W S / A D / Q E | Nicken / Gieren / Rollen (Drehrate bis 15°/s bzw. 30°/s; Eingriff schaltet den Autopiloten ab) |
| T | Lagehalter an/aus (bremst die Drehung und hält dann die Lage) |
| P | Autopilot: Aufstieg in eine 100-km-Kreisbahn |
| , / . | Zeitraffer langsamer / schneller (1× bis 100 000×) |
| M | Karte (weit herauszoomen) |
| Rechte Maus, Mausrad | Kamera drehen, zoomen |

Befehlstasten (Z, X, P, T, M, R, Komma, Punkt) gelten nach ihrer Beschriftung,
funktionieren also auf QWERTZ wie auf QWERTY. Gehaltene Tasten (WASD, QE,
Shift, Strg) gelten nach ihrer Position.
| R | Neustart |

## Aufbau

| Crate | Inhalt |
|---|---|
| `sim-core` | f64, typisierte Bezugssysteme (`Vec3<Inertial>`, `Vec3<VesselFrame>`, …), Keplerbahnen über universelle Variablen, RK4. Keine Engine-Abhängigkeit. |
| `vessel` | Teile aus RON-Dateien, Schiff als Baum, Masse/Schwerpunkt/Trägheit, Lasten in den Verbindungen, Staging, Trennung, Bewegung. |
| `zorbel-game` | Darstellung mit Bevy: Planet, Schiffe, Bahnlinien, Kamera, Anzeige, Steuerung, Screenshot-Modus. Rendert in f32 relativ zum Kamerafokus. |
| `flight-demo` | Headless-Testflug durch die Atmosphäre in die Umlaufbahn; die Unterstufe fällt ungesteuert zurück. |
| `assets/parts/` | Die fünf Basisteile: Kapsel, Tank, Triebwerk, Steuerfläche, Entkoppler. |

## Grundentscheidungen

- **Starrer Körper statt Federn.** Ein Schiff ist ein einziger starrer
  Körper. Verbindungen sind nur ein Graph; der Bewegungs-Solver kennt sie
  nicht. Wackeln ist dadurch ausgeschlossen.
- **Lasten aus der Bewegung.** Pro Verbindung werden Kraft und Moment
  berechnet, die der Teilbaum darunter braucht, um die Starrkörperbewegung
  mitzumachen (O(n), ein Durchlauf). Zerlegt in Zug/Druck, Scherung,
  Biegung und Torsion, verglichen mit der Festigkeit des schwächeren
  Partners. Über 100 % bricht die Verbindung und das Schiff wird in zwei
  starre Körper getrennt.
- **Gravitation erzeugt keine inneren Lasten** (homogenes Feld), deshalb
  wird mit der Eigenbeschleunigung gerechnet.
- **Eine Physik, ein Integrator.** `dynamics::step` gilt überall: am Boden,
  in der Luft, im Vakuum, im Zeitraffer. Verfahren Kick-Drift-Kick: die
  Gravitation wird exakt über die Keplerbahn gelöst, alle anderen Kräfte
  (Schub, Luft, Reaktionsrad) wirken als Stöße. Ohne diese Kräfte ist der
  Schritt exakt die Keplerbahn und darf beliebig lang sein.
- **Atmosphäre immer im Modell.** Jeder Himmelskörper hat eine; ohne
  Lufthülle ist es `Atmosphere::VACUUM` (Bodendruck 0). Der Druck geht am
  Oberrand glatt (C²) auf null. Luftkräfte und Triebwerksschub hängen
  stetig von Dichte und Druck ab. Ein Test zeigt: der Flug durch immer
  dünnere Luft konvergiert linear gegen den Vakuumflug.
- **Nur die Schrittweite passt sich an** (`suggest_dt`): kurz bei starken
  Kräften oder schneller Drehung, lang im kräftefreien Flug, und ein
  langer Schritt endet spätestens am Atmosphärenrand, damit er nicht durch
  die Lufthülle springt. Das ist Numerik, keine andere Physik.
- **Aerodynamik pro Teil.** Jedes Teil bekommt seine Luftkraft am eigenen
  Ort: Stirnflächen mit Verdeckung im Stapel, Kegelmantel, Queranströmung,
  Platten mit Anstellwinkel, Transschall-Anstieg. Stabilität, Dämpfung und
  Biegelasten folgen daraus ohne Sonderregeln. Steuerflächen bestimmen ihre
  Wirkrichtung selbst aus der aktuellen Anströmung.
- **Triebwerke** haben einen festen Massenstrom; der Schub folgt dem Isp,
  der linear vom Umgebungsdruck abhängt (am Boden 188 kN, im Vakuum 215 kN).
- **Boden als Zwangsbedingung, nicht als Feder.** Pro Schritt wird zuerst
  geprüft, ob der Boden das Schiff exakt zur Ruhe bringen kann (Normalkraft
  ≥ 0, Reibung im Reibkegel, Druckmittelpunkt in der Standfläche). Wenn ja,
  steht es ohne jede Drift; sonst (Rutschen, Kippen, Abheben) löst ein
  iteratives Verfahren die Kontaktstöße. Die Stützkraft geht wie jede äußere
  Kraft in die Lastberechnung ein; Stöße werden über die Nachgiebigkeit der
  Struktur (50 ms) verteilt, damit die Last nicht von der Schrittweite
  abhängt. Jedes Teil hat eine Aufprallgrenze (`crash_tolerance`).
- **Planetenrotation.** Der Planet dreht sich um seine Achse; Startplatz,
  Luft und Boden ruhen im mitdrehenden System `BodyFixed`. Die Luftkräfte
  hängen von der Anströmung relativ zur mitdrehenden Luft ab, auch an jedem
  Teil einzeln. Ein Start nach Osten am Äquator bringt 175 m/s mit.
  Drehrate 0 ist kein Sonderfall.
- **Typisierte Bezugssysteme.** Vektoren verschiedener Systeme zu
  addieren ist ein Compile-Fehler (siehe Doctest in `sim-core/src/frame.rs`).
- **Determinismus.** `BTreeMap` statt `HashMap`, Schrittweiten nur aus dem Zustand berechnet,
  feste Teilereihenfolge (Elternindex < Kindindex).

## Bekannte Vereinfachungen

- Nur ein Himmelskörper (SOI-Wechsel folgt).
- Keine Erwärmung beim Wiedereintritt (Thermik folgt).
- Steuerflächen schlagen ohne Stellgeschwindigkeit aus.
- Aerodynamik einfach gehalten: kein Auftrieb von Rumpfteilen, keine
  Verdeckung durch Radialteile, ein Luftzustand für das ganze Schiff.
- Der Boden ist die glatte Kugel mit mittlerem Radius (kein Gelände).
- Kontaktpunkte sind die Ränder der Stirnflächen (Achteck) und die Ecken
  der Platten, keine Mantelflächen.
- Ein Aufprall über der Grenze zerstört das ganze Schiff, nicht nur das Teil.
- Keine Kollisionen zwischen Schiffen.
- Trägheit von Kegelstümpfen genähert (volumengleicher Zylinder).
- Der Schubvektor nimmt an, dass das Triebwerk entlang der Schiffsachse sitzt.

## Nächste Schritte

1. Zweiter Himmelskörper und SOI-Übergänge.
2. Wiedereintrittserwärmung.
3. Gelände, Landebeine, Zerstörung einzelner Teile beim Aufprall.
