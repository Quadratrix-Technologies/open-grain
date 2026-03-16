# Open Grain — Research Notes

Physics decisions, literature values, and community sources compiled during development.
This file is meant to be read alongside the code and updated as the model improves.

---

## People This Project Honours

### Ron Mowrey ("Photo Engineer") — 1937–2020
32 years at Eastman Kodak. Invented Ektacolor 30/37 papers, first Kodak blix (1969),
Kodacolor Gold 400, named inventor on Kodachrome K-14 patent. Known universally as
"Photo Engineer" on APUG and Photrio (29,410 posts). Taught hands-on emulsion workshops
at The Photographers' Formulary, Montana. Published: *Photographic Emulsion Making,
Coating and Testing* (book + DVD). Democratised industrial emulsion science for the
DIY community. Died February 15, 2020.

### Grant Milford Haist — 1922–2015
PhD Physical Chemistry, Michigan State 1949. 32 years at Eastman Kodak Research Labs,
28 US patents. Author of *Modern Photographic Processing* Vols. 1 & 2 (Wiley, 1979) —
the definitive reference on photographic processing chemistry. Fellow of the Royal
Photographic Society, PSA, and IS&T.

### L.F.A. Mason (Leslie Frederick Alfred Mason)
Ilford Ltd, UK. Author of *Photographic Processing Chemistry* (Focal Press, 1966; rev.
1975) — described as the first contemporary work devoted exclusively to the chemistry of
silver-halide photographic materials. Covers nucleation, crystal growth, precipitation
chemistry. Cited in US patents.

---

## Key References

### Patents
- **US5422825A** — Supersaturation & halide ion concentration control during
  precipitation. Defines S = [Ag⁺][X⁻] / Ksp. Cascaded control with Vx (halide) and
  Vs (supersaturation) feedback loops. Uses γ = 140 erg/cm² = 0.14 J/m² for AgBr
  surface energy in nucleation/ripening calculations.
- **US5248577A** — Flow rate computation via smoothed halide concentration estimation.
  Least-squares parameter estimation, adaptive control law for optimal halide feed rate.
  Mass balance: C_x(t) = [initial + feed − consumed] / volume.

### Nucleation & Crystal Growth
- **Margolis & Gutoff (1974)** — AIChE J. 20(3). Foundational AgBr precipitation
  simulation paper. CNT-based nucleation model. Key finding: simulation is relatively
  insensitive to pre-exponential A, but sensitive to growth constants.
- **Jagannathan & Wey (1985)** — J. Crystal Growth 73, 226–232. Experimental AgBr
  nucleation study. Nucleus diameter 8–15 nm confirmed by EM. Diffusion-controlled
  growth confirmed.
- **Nielsen & Söhnel (1971)** — J. Crystal Growth 11, 233–242. Interfacial tensions for
  electrolyte crystals from nucleation data. γ vs. log(Ksp) linear correlation.
- **Sugimoto & Shiba (1999)** — J. Phys. Chem. B. Absolute interfacial energy of AgCl,
  AgBr, AgI from Gibbs-Thomson / potentiometric method. "Clean interface" γ ≈ 0.68 J/m²
  for AgBr — likely an upper bound; not the operative value in gelatin precipitation.
- **Sugimoto (gelatin system, ~2000)** — Colloids Surf. A. Spontaneous nucleation of
  monodisperse AgBr from gelatin solution. Measured: Sm = 3.61 (critical supersaturation),
  r* = 1.97 nm (critical nucleus radius), ΔG* = 1.79×10⁻¹⁸ J. These back-calculate to
  γ ≈ 0.14 J/m², consistent with the Kodak patent value.
- **Muhr et al. (1995)** — Chem. Eng. Sci. 50, 345–355. Updated AgBr precipitation
  model. A fitted as free parameter; γ controls nucleation burst shape.

### Ostwald Ripening
- **Lifshitz-Slyozov-Wagner (LSW) theory** — Classical coarsening theory. Predicts
  mean radius grows as r̄ ∝ t^(1/3) during diffusion-limited Ostwald ripening.
  Applicable to the ammonia digest phase.

### Community Sources
- **Photrio — "A Real Formula"** (2006, Ron Mowrey):
  https://www.photrio.com/forum/threads/a-real-formula.21896/
- **UnblinkingEye — Film Emulsion (Kodak AJ-12)**:
  https://unblinkingeye.com/Articles/Emulsion/emulsion.html
- **Photrio — SRAD emulsions thread**:
  https://www.photrio.com/forum/threads/srad-emulsions.42530/
- **Photrio — Emulsion design software**:
  https://www.photrio.com/forum/threads/emulsion-design-software.28475/

---

## Physical Constants for AgBr

| Quantity | Value | Source |
|---|---|---|
| Molar mass | 187.77 g/mol | Fundamental |
| Crystal density | 6.473 g/cm³ | Crystran / PubChem |
| Molar volume Vm | 2.901×10⁻⁵ m³/mol (29.01 cm³/mol) | Derived from above |
| Ksp (25°C) | 5.0×10⁻¹³ mol²/L² | Standard reference |
| Ksp (40°C) | ~6.3×10⁻¹³ mol²/L² | Approx. temp. correction |
| Surface energy γ (gelatin/aqueous) | **0.14 J/m²** | Kodak US5448160; Sugimoto gelatin |
| Surface energy γ (clean interface) | 0.68 J/m² | Sugimoto & Shiba 1999 |
| Critical supersaturation Sm | 3.61 | Sugimoto gelatin nucleation expt. |
| Critical nucleus radius r* | 1.97 nm | Sugimoto gelatin nucleation expt. |
| Nucleation barrier ΔG* at Sm | 1.79×10⁻¹⁸ J | Sugimoto gelatin nucleation expt. |

---

## Nucleation Model

### Current implementation: Classical Nucleation Theory (CNT)
```
J = A × exp(−B_eff / (ln S)²)   [nuclei / L / s]
```

**B_hom** is computed from first principles:
```
B_hom = (16π/3) × (γ/kT)³ × v²
```
where v = Vm / Nₐ = 4.818×10⁻²⁹ m³/molecule.

At T = 318K (45°C), γ = 0.14 J/m²: B_hom ≈ 1260.
B_eff = B_hom × CNT_HET_FACTOR = 1260 × 0.01 = 12.6

**A** — pre-exponential (kinetic/calibration term):
- Theoretical homogeneous: 10²⁵–10³⁰ m⁻³/s
- **Operative value: 1e20 m⁻³/s** — calibrated so that nucleation rate × plume volume
  gives ~10¹³–10¹⁴ nuclei/L over a 10-min SRAD run
- Margolis & Gutoff note insensitivity to A vs B; this fits AJ-12 grain counts

**S_threshold (Sm)**: 3.61 (Sugimoto gelatin experiments — replaces arbitrary 2.0)

**Plume volume correction**: nuclei per step = J × dt × V_plume where
V_plume = (molsBr_this_step × mixingFactor) / halidConc [L].
Without this, J (nuclei/L/s) × dt gives nuclei/L which is dimensionally wrong.

**Dissolution back-reaction**: when S_bulk < 1, a tiny amount of AgBr dissolves
to restore equilibrium. Amount ≈ Ksp/brC mol/L (< 1e-9 mol/L typically).
Formula: ΔC = 2(Ksp − agC·brC) / (agC+brC+√((agC−brC)²+4Ksp))
This ensures correct final pAg (~11.5 for AJ-12) rather than crashing to 99.

### Mixing / Nucleation Cap
Nucleation per timestep is physically capped by the jet plume — the high-S zone
containing only freshly added Br⁻ that hasn't yet mixed into the bulk. Cap:
```
maxNuclei = (molsBrThisStep × mixingFactor) / molesPerNucleus
mixingFactor = mixingTime / dt
```
`mixingTime` (τ, seconds) is the plume dispersal time — how long freshly added halide
persists as a concentrated zone before bulk mixing. Physically motivated by US5422825A's
controlled-addition approach. Typical values:
- ~1–3s: vigorous overhead stirring
- ~5–10s: magnetic stirbar, moderate
- ~30s+: manual / slow mixing

---

## Known Emulsion Presets

### 1. PE SRAD / Kodak AJ-12 (ammonia)
Source: Ron Mowrey (Photo Engineer), Photrio 2006. Identified as Kodak Publication AJ-12.
Also reproduced at UnblinkingEye.com and in Reed & Jones, *Silver Gelatin*.

Original recipe (B into A = silver into halide, opposite of sim convention):
- Solution A: KBr 132g (1.11M), KI 4.5g, Gelatin 30g, Water 1000mL
- Solution B: AgNO₃ 130g (1.53M), Water 500mL + NH₄OH to clear
- B poured into A over 10 minutes at 45°C
- Digest: 30 min at 45°C
- Expected ISO: 3–40 depending on gelatin activity
- Expected grain size: 0.2–1.0 µm (rounded AgBrI cuboctahedra, ~3.3% iodide)

### 2. AJ-12 no-ammonia
Source: Kodak AJ-12 via UnblinkingEye.com
- AgNO₃: 40g in 400mL (0.59M)
- KBr: 32g + KI: 0.8g
- Addition: 20mL per 30s over 10min at 55°C
- Ripening: 10 + 15 min at 55°C
- Finer grain, lower speed than ammonia version

### 3. Ludwik 2024 (small batch)
Source: Photrio 2024, user "ludwik"
- AgNO₃: 10g in 100mL (0.59M)
- KBr: 12g in 100mL (1.01M) + trace KI
- Temperature: 55°C
- Long thermal digest: up to 3 hours
- Reported: ~ISO 30

### 4. LeoniD core-shell (double-jet)
Source: Photrio 2023, user "LeoniD"
- AgNO₃: 20g in 60mL (1.96M) + NH₄OH to clear
- KBr: 20.3g in 31mL (5.51M)
- Sequential KI shell additions (2×0.23g in 5mL each)
- Growth flows: B 6mL/min, C 3.1mL/min
- 2000 rpm nucleation step
- 45°C, 30 min digest

### 5. Fulvio baseline
Source: Photrio 2006, user "Fulvio" + Photo Engineer modifications
- AgNO₃: 20g in 125mL (0.94M)
- KBr: 16g in 125mL (1.07M)
- 50°C, ~10 min ripening
- ~ISO 3 (no ammonia)

---

## Planned Outputs / Validation Targets

- pAg trajectory: should span ~0 → ~11 over precipitation for AJ-12
- Nucleation: single burst in first ~30s, then stops (LaMer model)
- Grain size at end of precipitation: ~50–150 nm diameter (pre-digest)
- Grain size after 30 min ammonia digest: **200–800 nm** diameter
- LSW check: plot r̄³ vs time during digest — should be linear
- Grain count: ~10¹² to 10¹⁴ per litre
- Crystal habit: Mixed/Intermediate to Octahedral {111} for AJ-12 (low final pAg
  with bromide excess + low iodide)

---

## ISO / Speed Notes

- Grain size is primary determinant: speed ∝ grain area (diameter²)
- Tabular (T-grain) more efficient: large cross-section per unit silver volume
- Precipitation determines speed *potential* only
- Chemical sensitization (S+Au): +1–2 stops
- Spectral sensitization (cyanine dyes): determines colour response
- Simulator predicts *relative* speed differences between precipitation profiles

---

## TODO / Known Issues

- [x] Replace LaMer nucleation with CNT formula
- [x] Replace K_GROWTH-based growth with equilibrium precipitation approach
- [x] Add dissolution step to restore S=1 when brC builds up in excess
- [ ] Calibrate CNT_A more precisely against published grain count data
- [ ] Add iodide incorporation model (for double-jet phase 2)
- [ ] Verify LSW t^(1/3) law during digest phase
- [ ] Add temperature dependence to growth and ripening constants
- [ ] Consider adding Br⁻ excess / pAg setpoint control mode (per US5422825A)
