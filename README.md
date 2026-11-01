# Zorbel Space Program

Moderne Umsetzung der Kernmechaniken von Kerbal Space Program in Rust.
Stand: Simulationskern mit Atmosphäre, Teilen, Schiffsgraph, Aerodynamik
pro Teil, Lastberechnung und Staging; ein Universum aus zwei Sternsystemen
mit Wechsel der Einflusssphäre; Darstellung mit Bevy.

## Das Universum

Definiert in `assets/universe.ron` (Galaxie → Systeme → Stern → Planeten → Monde):

| System | Stern | Planeten / Monde | Lage |
|---|---|---|---|
| Zol | Zun | Zorb (Heimat, Atmosphäre) mit Mond Muz (luftlos, gebunden rotierend) | Ursprung |
| Zirkon (Platzhalter) | Zirka | Zirk (dünne Atmosphäre) | 3,7 Lichtjahre entfernt, eigene Bewegung 17 km/s |

Zwischen den Systemen wird später gesprungen, nur aus einem Sprungring um
den Stern (`jump_zone`). Galaxie-Positionen sind ganzzahlig in Mikrometern
gespeichert (i128), damit auch Lichtjahre exakt bleiben.

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
| C | Schwerpunkt (gelb) und Druckpunkt (blau) anzeigen, mit Stabilitätsangabe |
| F | Neustart mit/ohne Finnen |
| R | Neustart |

## Aufbau

| Crate | Inhalt |
|---|---|
| `sim-core` | f64, typisierte Bezugssysteme (`Vec3<Inertial>`, `Vec3<VesselFrame>`, …), Keplerbahnen über universelle Variablen, Universum (Körperbaum, Bahnen, Einflusssphären, Galaxie-Positionen). Keine Engine-Abhängigkeit. |
| `vessel` | Teile aus RON-Dateien, Schiff als Baum, Masse/Schwerpunkt/Trägheit, Lasten in den Verbindungen, Staging, Trennung, Bewegung, Flug im Universum (`flight`: Wechsel der Einflusssphäre). |
| `zorbel-game` | Darstellung mit Bevy: Himmelskörper des aktuellen Systems (Stern als Lichtquelle), fremde Sterne am Himmel, Schiffe, Bahnlinien von Schiffen und Körpern, Kamera, Anzeige, Steuerung, Screenshot-Modus. Rendert in f32 relativ zum Kamerafokus. |
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
- **Drehung exakt, auch im Zeitraffer** (`vessel::rotation`). Zusammengebaute
  Schiffe sind unsymmetrisch; der Trägheitstensor wird pro Schritt auf
  Hauptachsen gebracht, die kräftefreie Drehung als symmetrisches Splitting
  aus drei exakten Drehungen um die Hauptachsen gelöst (2. Ordnung). Der
  Drehimpuls bleibt exakt erhalten, die Energie ohne Drift (Fehler 5·10⁻⁷
  über 10⁶ s). Lange Schritte werden in Teilschritte zerlegt, deren Weite aus
  einer Erhaltungsgröße folgt (|L|/I_min), sonst würde die Energie driften.
  Ein trudelndes Schiff trudelt im Zeitraffer genau wie in Echtzeit, samt
  Dzhanibekov-Effekt, und bremst den Zeitraffer nicht.
- **Lagehalter im Zeitraffer.** Bei langen Schritten bremsen die
  Reaktionsräder in jedem Teilschritt mit ihrem echten Moment, also gleich
  schnell wie in Echtzeit. Steht das Schiff, hält es die Lage.
- **Nur die Schrittweite passt sich an** (`suggest_dt`): kurz bei starken
  Kräften oder schneller Drehung, lang im kräftefreien Flug, und ein
  langer Schritt endet spätestens am Atmosphärenrand, damit er nicht durch
  die Lufthülle springt. Das ist Numerik, keine andere Physik.
- **Physikmodell getrennt von der Darstellung.** Jedes Teil hat in seiner
  Teiledatei eine einfache Hülle (`hull`: Zylinder, Kegelstumpf, Platte), die
  in Polygone zerlegt wird. Eigenschaften, die im Teil selbst stecken (Masse,
  Festigkeit, Aufprallgrenze, Triebwerk), sind Parameter. Alles, was von der
  Anordnung abhängt (Luftkräfte, Verdeckung, Bodenkontakt, später Wärme),
  folgt aus dem Polygonmodell des ganzen Schiffs. Optionale Faktoren
  `drag_factor`/`lift_factor` erlauben Abstimmung, ohne die Formel zu ändern.
- **Einheitliches Luftmodell** (`vessel::aero`, `vessel::geometry`): Verdeckung
  über das ganze Schiff per Tiefenpuffer aus 18 Richtungen (senkrecht
  angeströmt sind nur Spitze und Finnenkanten frei). Anströmung zerlegt in
  Längs- und Queranteil (Allen/Perkins): Druck nach modifizierter
  Newton-Theorie mit Machabhängigkeit, Sog auf Rückflächen, Reibung, und
  Auftrieb nach Schlankkörpertheorie dort, wo das Schiff breiter wird
  (Rumpf-Flügel-Fläche exakt, Seitenverhältnis nach Helmbold). Ruder drehen
  nur die Anströmung des Teils. Druck- und Schwerpunkt ergeben sich daraus;
  eine Rakete ohne Finnen ist instabil und schlägt um.
  Geprüft: Rumpf CNα 2,07 (Munk/Barrowman 2,0), Nasen-Druckpunkt nach
  Barrowman, Zylinder quer Cd 1,20, längs 0,93.
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
- **Patched Conics.** Jedes Schiff gehört zu genau einem Körper und rechnet
  relativ zu ihm. Einflusssphäre r = a·(μ/μ_Eltern)^0,4. Verlässt das Schiff
  sie, oder betritt es die eines Mondes, wird nur der Bezugspunkt umgerechnet
  (Lage und Geschwindigkeit relativ zum Stern bleiben exakt). Die
  Schrittweite endet höchstens 1 s hinter der Grenze; geprüft wird nur, wenn
  die Bahn die Grenze überhaupt erreichen kann, damit Zeitraffer in einer
  stabilen Bahn nicht gebremst wird. Die Körper selbst folgen festen
  Keplerbahnen (schienengeführt wie in KSP).
- **Typisierte Bezugssysteme.** Vektoren verschiedener Systeme zu
  addieren ist ein Compile-Fehler (siehe Doctest in `sim-core/src/frame.rs`).
- **Determinismus.** `BTreeMap` statt `HashMap`, Schrittweiten nur aus dem Zustand berechnet,
  feste Teilereihenfolge (Elternindex < Kindindex).

## Bekannte Vereinfachungen

- Patched Conics: immer nur ein Körper zieht am Schiff.
- Keine Erwärmung beim Wiedereintritt (Thermik folgt).
- Steuerflächen schlagen ohne Stellgeschwindigkeit aus.
- Finnenauftrieb nach Schlankkörpertheorie liegt etwa 1,45-fach über
  Barrowman (Rumpf-Mitnahme ohne Heck dahinter wird überschätzt); Raketen
  sind dadurch etwas stabiler als in der Wirklichkeit.
- Auftrieb machunabhängig; kein Wellenwiderstand nach Flächenregel; ein
  Luftzustand für das ganze Schiff; Hüllen aus je einem Grundkörper.
- Der Boden ist die glatte Kugel mit mittlerem Radius (kein Gelände).
- Kontaktpunkte sind die Ränder der Stirnflächen (Achteck) und die Ecken
  der Platten, keine Mantelflächen.
- Ein Aufprall über der Grenze zerstört das ganze Schiff, nicht nur das Teil.
- Keine Kollisionen zwischen Schiffen.
- Trägheit von Kegelstümpfen genähert (volumengleicher Zylinder).
- Der Schubvektor nimmt an, dass das Triebwerk entlang der Schiffsachse sitzt.

## Nächste Schritte

Vereinbarte Reihenfolge:

1. ~~Universum mit zwei Systemen, Körperbaum, SOI-Übergänge~~ (erledigt).
2. ~~Kräftefreie Drehung für beliebige Trägheit, Lagehalter im Zeitraffer~~
   (erledigt).
3. Navball mit sieben Lagehalter-Modi.
4. Mehrere Ressourcen, Sprungantrieb und Sprungtreibstoff, Sprung zwischen
   den Systemen mit Vorschau der Ankunftsbahn.
5. Schiff als allgemeiner Graph (stabile Teile-IDs, Kantenliste), später
   Andocken und Stationen.

Danach: Wiedereintrittserwärmung, Gelände, Landebeine.
