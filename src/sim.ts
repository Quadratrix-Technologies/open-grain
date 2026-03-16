// sim.ts — Silver halide precipitation physics engine
// Population balance model with explicit r³ tracking per bin.
//
// References:
//   US5422825A — supersaturation control during precipitation
//   US5248577A — adaptive halide feed rate control
//   Lifshitz-Slyozov-Wagner (LSW) theory — Ostwald ripening
//   Classical Nucleation Theory (CNT) — Kashchiev (2000)
//   Margolis & Gutoff (1974) AIChE J. 20(3) — AgBr precipitation model
//   Sugimoto (2000) Colloids Surf. A — AgBr nucleation in gelatin: Sm=3.61, r*=1.97nm
//   Sugimoto & Shiba (1999) J. Phys. Chem. B — AgBr interfacial energy
//   Mason (1966) Photographic Processing Chemistry — grain size / sensitivity
//   Haist (1979) Modern Photographic Processing Vol.1 — emulsion speed relationships

// ─── Physical constants ───────────────────────────────────────────────────────

const BOLTZMANN  = 1.38065e-23; // J/K
const AVOGADRO   = 6.02214e23;  // mol⁻¹

// ─── AgBr / AgI material constants ───────────────────────────────────────────

// AgBr Ksp at 40°C (mol²/L²)
const KSP_AGBR_40C = 6.3e-13;


// Molar volume of AgBr — density 6.473 g/cm³, M = 187.77 g/mol
// Vm = 187.77 / 6.473 = 29.01 cm³/mol = 2.901e-5 m³/mol (Crystran/PubChem)
const VM_AGBR = 2.901e-5; // m³/mol

// Solid-liquid interfacial energy of AgBr in gelatin/aqueous solution
// 0.14 J/m²: Kodak patent US5448160; consistent with Sugimoto (2000) gelatin data
const GAMMA_AGBR = 0.14; // J/m²

// Heterogeneous nucleation correction factor for gelatin system.
// f(θ) = (2 + cosθ)(1 − cosθ)² / 4 ≈ 0.01 → θ ≈ 12° (favorable on gelatin)
// Calibrated to match onset at Sm = 3.61 (Sugimoto 2000).
const CNT_HET_FACTOR = 0.01;

// CNT pre-exponential factor A (m⁻³ s⁻¹).
// Theoretical homogeneous: 10²⁵–10³⁰ m⁻³/s (Margolis & Gutoff 1974).
// Fitted for gelatin-mediated heterogeneous nucleation.
const CNT_A = 1e20;

// Critical supersaturation for AgBr nucleation onset in gelatin (Sugimoto 2000).
const S_NUC_THRESHOLD = 3.61;

// ─── Ostwald ripening constant ────────────────────────────────────────────────
// K_LSW rate constant for dr³/dt = K × (1/r* − 1/r)   [m³/s]
//
// From first-principles: K = 8 D C_inf γ Vm² / (9 RT)
// For AgBr in gelatin at 45°C (D_gel ≈ 7e-10 m²/s, C_inf ≈ 8e-4 mol/m³, γ=0.14, Vm=2.901e-5):
//   K_no_ammonia ≈ 5e-26 m³/s
// Ammonia greatly enhances Ag+ solubility via Ag(NH3)+ complexation.
// ammoniaFactor scales K with [NH3]: calibrated so AJ-12 (0.3 M NH3, 45°C, 30 min)
// gives ~50% increase in mean grain volume during digest.
// K_LSW_BASE here is the no-ammonia effective value fitted to community data.
// K = 8 D C_inf γ Vm² / (9RT); D_gel ≈ 1e-13 m²/s for Ag+ in 30 g/L gelatin at 45°C
// → K_no_ammonia ≈ 5e-30 m³/s  (Muhr et al. 1995; gelatin diffusivity from Stokes-Einstein)
const K_LSW_BASE = 5e-30; // m³/s


// ─── Types ────────────────────────────────────────────────────────────────────

export interface SimParams {
  // Reactor
  reactorVolume: number;      // L
  temperature: number;        // °C
  gelatinConc: number;        // g/L (affects growth kinetics)

  // Silver nitrate
  agno3Conc: number;          // mol/L initial

  // Halide jet (KBr ± KI)
  halidConc: number;          // mol/L KBr
  kiConc: number;             // mol/L KI in halide jet (0 = no iodide)
                              // AJ-12: 4.5 g KI in 1000 mL → 0.027 M
  jetFlowRate: number;        // mL/min
  jetDuration: number;        // seconds

  // Silver jet — double-jet mode (agJetFlowRate > 0 to enable)
  agJetConc: number;          // mol/L — AgNO₃ jet concentration
  agJetFlowRate: number;      // mL/min — 0 = SRAD single-jet

  // Ammonia digest
  ammoniaConc: number;        // mol/L
  digestDuration: number;     // seconds

  // Reactor mixing
  mixingTime: number;         // seconds — plume dispersal τ (~1s vigorous, ~10s gentle)
}

export interface SimState {
  time: number;               // seconds
  agConc: number;             // [Ag+] mol/L
  brConc: number;             // [Br-] mol/L
  volume: number;             // L
  pAg: number;                // −log10([Ag+])
  supersaturation: number;    // S = [Ag+][Br-] / Ksp
  bins: Float64Array;         // grain count per radius bin (for histogram)
  nucleationRate: number;     // nuclei/s (from plume)
  meanRadius: number;         // m  (mass-weighted, from actual r³ tracking)
  stdRadius: number;          // m
  totalGrains: number;
  iodideMolPct: number;       // mol% iodide in grain lattice (AgI / (AgBr+AgI) × 100)
  phase: 'precipitation' | 'digest';
}

export interface SimResult {
  params: SimParams;
  timepoints: SimState[];
  bins: { edges: Float64Array; counts: Float64Array };
}

// ─── Bin geometry ─────────────────────────────────────────────────────────────

const N_BINS = 80;
const R_MIN = 1e-9;   // 1 nm
const R_MAX = 2e-6;   // 2 µm

function makeBinEdges(): Float64Array {
  const edges = new Float64Array(N_BINS + 1);
  const logMin = Math.log10(R_MIN);
  const logMax = Math.log10(R_MAX);
  for (let i = 0; i <= N_BINS; i++)
    edges[i] = Math.pow(10, logMin + (i / N_BINS) * (logMax - logMin));
  return edges;
}

function binCenters(edges: Float64Array): Float64Array {
  const c = new Float64Array(N_BINS);
  for (let i = 0; i < N_BINS; i++)
    c[i] = Math.sqrt(edges[i] * edges[i + 1]); // geometric mean
  return c;
}

// ─── Physics helpers ──────────────────────────────────────────────────────────

function ksp(tempC: number): number {
  return KSP_AGBR_40C * Math.pow(10, 0.012 * (tempC - 40));
}


function supersaturation(agC: number, brC: number, tempC: number): number {
  if (agC <= 0 || brC <= 0) return 0;
  return (agC * brC) / ksp(tempC);
}

function pAg(agC: number): number {
  if (agC <= 0) return 99;
  return -Math.log10(agC);
}

// CNT nucleation rate: J = A × exp(−B_eff / (ln S)²)  [nuclei / L / s]
// B_hom = (16π/3) × (γ/kT)³ × v²  where v = Vm/Nₐ
// B_eff = B_hom × CNT_HET_FACTOR
function nucleationRate(S: number, tempC: number): number {
  if (S <= S_NUC_THRESHOLD) return 0;
  const TK   = tempC + 273.15;
  const kT   = BOLTZMANN * TK;
  const v    = VM_AGBR / AVOGADRO;
  const B_hom = (16 * Math.PI / 3) * Math.pow(GAMMA_AGBR / kT, 3) * v * v;
  const B_eff = B_hom * CNT_HET_FACTOR;
  const lnS  = Math.log(S);
  return CNT_A * Math.exp(-B_eff / (lnS * lnS)) * 1e-3; // nuclei / L / s
}

// Numerically stable equilibrium ΔC: solves (agC−ΔC)(brC−ΔC) = Ksp
// Standard form would lose precision when agC ≪ brC; this form avoids that.
function precipDeltaC(agC: number, brC: number, tempC: number): number {
  const kspT = ksp(tempC);
  const diff = agC - brC;
  const disc = Math.sqrt(diff * diff + 4 * kspT);
  return 2 * (agC * brC - kspT) / (agC + brC + disc);
}

function dissolveDeltaC(agC: number, brC: number, tempC: number): number {
  const kspT = ksp(tempC);
  const diff = agC - brC;
  return 2 * (kspT - agC * brC) / (agC + brC + Math.sqrt(diff * diff + 4 * kspT));
}

// ─── Simulator ────────────────────────────────────────────────────────────────

export function runSimulation(params: SimParams): SimResult {
  const edges  = makeBinEdges();
  const centers = binCenters(edges);

  const dt = 0.5; // seconds per step
  const totalTime = params.jetDuration + params.digestDuration;
  const nSteps = Math.ceil(totalTime / dt);

  // State
  let agC = params.agno3Conc;     // mol/L
  let brC = 0.0;
  let vol = params.reactorVolume; // L
  const tempC = params.temperature;

  // Population balance: counts and actual r³ per bin
  // binR3[i] = sum of r³ for all grains currently in bin i  (m³, not normalised)
  // Per-grain representative r: Math.cbrt(binR3[i] / bins[i])
  // This eliminates the bin-centre mass accounting error where intra-bin
  // growth would deplete agC without updating the reported grain mass.
  const bins  = new Float64Array(N_BINS); // grain counts
  const binR3 = new Float64Array(N_BINS); // sum(r³) per bin [m³]

  // Iodide tracking
  let totalIodideMol = 0.0;   // cumulative moles of I in grain lattice
  let totalHalideMol = 0.0;   // cumulative moles of halide precipitated (Br + I)

  const timepoints: SimState[] = [];
  const recordEvery = Math.max(1, Math.floor(nSteps / 400));

  // Nucleation geometry — first bin centre radius
  const r0  = centers[0];
  const r0_3 = r0 * r0 * r0;
  const molesPerNucleus = ((4 / 3) * Math.PI * r0_3) / VM_AGBR;

  for (let step = 0; step < nSteps; step++) {
    const t = step * dt;
    const phase: 'precipitation' | 'digest' =
      t < params.jetDuration ? 'precipitation' : 'digest';

    // ── Jet addition, nucleation, iodide conversion ────────────────────────
    if (phase === 'precipitation') {
      const jetFlowL   = (params.jetFlowRate   / 1000 / 60) * dt;
      const agJetFlowL = (params.agJetFlowRate / 1000 / 60) * dt;
      const molsBrThisStep  = jetFlowL   * params.halidConc;
      const molsKiThisStep  = jetFlowL   * params.kiConc;   // KI rides in halide jet
      const molsAgJet       = agJetFlowL * params.agJetConc;

      // Dilute existing bulk for volume increase
      const newVol = vol + jetFlowL + agJetFlowL;
      agC = (agC * vol) / newVol;
      brC = (brC * vol) / newVol;
      vol = newVol;

      const mixingFactor = dt / params.mixingTime;

      // ── Br⁻ plume nucleation ─────────────────────────────────────────────
      const S_plume_br = supersaturation(agC, params.halidConc, tempC);
      let molsBrNucleated = 0;
      if (S_plume_br > S_NUC_THRESHOLD) {
        const maxFromPlume_br = (molsBrThisStep * mixingFactor) / molesPerNucleus;
        const maxFromAg       = (agC * vol) / molesPerNucleus;
        const V_plume_br      = (molsBrThisStep * mixingFactor) / params.halidConc;
        const newNuclei_br    = Math.min(
          nucleationRate(S_plume_br, tempC) * dt * V_plume_br,
          maxFromPlume_br, maxFromAg,
        );
        bins[0]         += newNuclei_br;
        binR3[0]        += newNuclei_br * r0_3;
        molsBrNucleated  = newNuclei_br * molesPerNucleus;
        agC              = Math.max(0, agC - molsBrNucleated / vol);
        totalHalideMol  += molsBrNucleated;
      }

      // ── Ag⁺ plume nucleation (double-jet only) ───────────────────────────
      let molsAgNucleated = 0;
      if (params.agJetFlowRate > 0 && molsAgJet > 0) {
        const S_plume_ag = supersaturation(params.agJetConc, brC, tempC);
        if (S_plume_ag > S_NUC_THRESHOLD) {
          const maxFromPlume_ag = (molsAgJet * mixingFactor) / molesPerNucleus;
          const maxFromBr       = (brC * vol) / molesPerNucleus;
          const V_plume_ag      = (molsAgJet * mixingFactor) / params.agJetConc;
          const newNuclei_ag    = Math.min(
            nucleationRate(S_plume_ag, tempC) * dt * V_plume_ag,
            maxFromPlume_ag, maxFromBr,
          );
          bins[0]          += newNuclei_ag;
          binR3[0]         += newNuclei_ag * r0_3;
          molsAgNucleated   = newNuclei_ag * molesPerNucleus;
          brC               = Math.max(0, brC - molsAgNucleated / vol);
          totalHalideMol   += molsAgNucleated;
        }
      }

      // Leftover jet material enters bulk
      brC += Math.max(0, molsBrThisStep - molsBrNucleated) / vol;
      agC += Math.max(0, molsAgJet      - molsAgNucleated) / vol;

      // ── Iodide conversion ────────────────────────────────────────────────
      // I⁻ + AgBr(surface) → AgI + Br⁻   K_eq = Ksp(AgBr)/Ksp(AgI) ≈ 6000
      // Reaction is fast relative to simulation dt; treat as instantaneous.
      // The released Br⁻ increases bulk brC.
      // Source: Mason (1966) Photographic Processing Chemistry, Chapter 3;
      //         Zelikman & Levi (1964) Making and Coating Photographic Emulsions.
      if (molsKiThisStep > 0) {
        let hasGrains = false;
        for (let i = 0; i < N_BINS; i++) if (bins[i] >= 1) { hasGrains = true; break; }
        if (hasGrains) {
          // All I⁻ converts to AgI on grain surfaces, releasing equimolar Br⁻
          brC += molsKiThisStep / vol;
          totalIodideMol += molsKiThisStep;
        } else {
          // No grains yet; I⁻ stays in solution until grains form
          brC += molsKiThisStep / vol;
        }
      }
    }

    // ── Equilibrium precipitation growth ──────────────────────────────────
    // Solve (agC − ΔC)(brC − ΔC) = Ksp for ΔC; distribute moles to grains
    // proportional to surface area r².  Grains' actual r tracked via binR3.
    const S_bulk = supersaturation(agC, brC, tempC);
    if (S_bulk > 1.0) {
      const deltaC = precipDeltaC(agC, brC, tempC);
      const maxDeltaC = Math.min(deltaC, agC, brC);
      const totalMolsToPrecip = maxDeltaC * vol;

      let totalArea = 0;
      for (let i = 0; i < N_BINS; i++) {
        if (bins[i] >= 1) {
          const r = Math.cbrt(binR3[i] / bins[i]);
          totalArea += bins[i] * r * r;
        }
      }

      if (totalArea > 0 && totalMolsToPrecip > 0) {
        let molsConsumed = 0;
        const snapCts = new Float64Array(bins);
        const snapR3  = new Float64Array(binR3);

        for (let i = 0; i < N_BINS; i++) {
          if (snapCts[i] < 1) continue;
          const r = Math.cbrt(snapR3[i] / snapCts[i]);

          // Moles assigned to this bin, weighted by surface area
          const molsThisBin = totalMolsToPrecip * (snapCts[i] * r * r) / totalArea;
          molsConsumed += molsThisBin;

          // New r³ per grain after growth
          const dVolPerGrain = (molsThisBin / snapCts[i]) * VM_AGBR;
          const r3New = r * r * r + dVolPerGrain * (3 / (4 * Math.PI));
          const newR  = Math.cbrt(r3New);

          if (newR >= R_MAX) {
            bins[i]  -= snapCts[i];
            binR3[i] -= snapR3[i];
            continue;
          }

          let destBin = i;
          while (destBin < N_BINS - 1 && edges[destBin + 1] < newR) destBin++;

          if (destBin !== i) {
            bins[i]        -= snapCts[i];
            binR3[i]       -= snapR3[i];
            bins[destBin]  += snapCts[i];
            binR3[destBin] += snapCts[i] * r3New;
          } else {
            // Grain stays in same bin — update r³ in place.
            // Replace original grains' r³ contribution with post-growth value.
            // (Grains moved in from lower bins during this sweep retain their r³.)
            binR3[i] = binR3[i] - snapR3[i] + snapCts[i] * r3New;
          }
        }

        totalHalideMol += molsConsumed;
        agC = Math.max(0, agC - molsConsumed / vol);
        brC = Math.max(0, brC - molsConsumed / vol);
      }
    }

    // ── Ostwald ripening — digest phase (LSW) ─────────────────────────────
    // dr³/dt = kLSW × (1/r* − 1/r)   where r* = number-weighted mean radius
    // ammonia greatly increases Ag+ solubility via Ag(NH3)+ complexation;
    // ammoniaFactor scales kLSW accordingly.
    if (phase === 'digest') {
      const ammoniaFactor = 1 + 10 * params.ammoniaConc;
      const kLSW = K_LSW_BASE * ammoniaFactor;

      const dt_sub = dt;

      // LSW: only dissolve sub-critical grains (r < r*).
      // Released material returns to agC/brC and is redistributed by the
      // equilibrium growth step, which deposits proportionally to grain area.
      // This is equivalent to Ostwald ripening but avoids the runaway growth
      // caused by directly applying the Kelvin-equation growth to super-critical
      // grains in a polydisperse population.

      // Compute number-weighted mean radius as the critical radius
      let sm = 0, nt = 0;
      for (let i = 0; i < N_BINS; i++) {
        if (bins[i] < 1) continue;
        nt += bins[i];
        sm += Math.cbrt(binR3[i] / bins[i]) * bins[i];
      }
      const rCrit = nt > 0 ? sm / nt : centers[0];

      const snapLswCts = new Float64Array(bins);
      const snapLswR3  = new Float64Array(binR3);

      for (let i = 0; i < N_BINS; i++) {
        if (snapLswCts[i] < 1) continue;
        const r     = Math.cbrt(snapLswR3[i] / snapLswCts[i]);
        if (r >= rCrit) continue;  // super-critical: growth handled by equilibrium step

        const dr3dt = kLSW * (1 / rCrit - 1 / r);  // negative for r < rCrit
        const r3new = r * r * r + dr3dt * dt_sub;

        if (r3new <= 0) {
          // Grain fully dissolves
          const molsReleased = snapLswCts[i] * (4 / 3) * Math.PI * (r * r * r) / VM_AGBR;
          agC += molsReleased / vol;
          brC += molsReleased / vol;
          totalHalideMol -= molsReleased;
          bins[i]  -= snapLswCts[i];
          binR3[i] -= snapLswR3[i];
        } else {
          // Grain partially shrinks — release the lost mass to solution
          const molsReleased = snapLswCts[i] * (4 / 3) * Math.PI * (r * r * r - r3new) / VM_AGBR;
          agC += molsReleased / vol;
          brC += molsReleased / vol;
          totalHalideMol -= molsReleased;

          const rNew = Math.cbrt(r3new);
          let destBin = i;
          while (destBin > 0 && edges[destBin] > rNew) destBin--;

          if (destBin !== i) {
            bins[i]        -= snapLswCts[i];
            binR3[i]       -= snapLswR3[i];
            bins[destBin]  += snapLswCts[i];
            binR3[destBin] += snapLswCts[i] * r3new;
          } else {
            binR3[i] = binR3[i] - snapLswR3[i] + snapLswCts[i] * r3new;
          }
        }
      }
    }

    // ── Dissolution equilibration ─────────────────────────────────────────
    // If S < 1 and grains exist, dissolve just enough AgBr to restore S = 1.
    // Amount is always tiny (≈ Ksp/brC ~ 1 µmol/L).
    const S_check = supersaturation(agC, brC, tempC);
    if (S_check < 1.0) {
      let hasGrains = false;
      for (let i = 0; i < N_BINS; i++) if (bins[i] >= 1) { hasGrains = true; break; }
      if (hasGrains) {
        const dissolveC = dissolveDeltaC(agC, brC, tempC);
        agC += dissolveC;
        brC += dissolveC;
      }
    }

    // ── Record state ─────────────────────────────────────────────────────
    if (step % recordEvery === 0) {
      const { mean, std, total } = distributionStats(bins, binR3);
      const iodMolPct = totalHalideMol > 0
        ? (totalIodideMol / totalHalideMol) * 100
        : 0;
      timepoints.push({
        time: t,
        agConc: agC,
        brConc: brC,
        volume: vol,
        pAg: pAg(agC),
        supersaturation: S_bulk,
        bins: new Float64Array(bins),
        nucleationRate: phase === 'precipitation'
          ? Math.max(
              nucleationRate(supersaturation(agC, params.halidConc, tempC), tempC),
              params.agJetFlowRate > 0
                ? nucleationRate(supersaturation(params.agJetConc, brC, tempC), tempC)
                : 0,
            )
          : 0,
        meanRadius: mean,
        stdRadius: std,
        totalGrains: total,
        iodideMolPct: iodMolPct,
        phase,
      });
    }
  }

  // Final snapshot
  const { mean, std, total } = distributionStats(bins, binR3);
  const iodMolPctFinal = totalHalideMol > 0
    ? (totalIodideMol / totalHalideMol) * 100
    : 0;
  timepoints.push({
    time: totalTime,
    agConc: agC,
    brConc: brC,
    volume: vol,
    pAg: pAg(agC),
    supersaturation: supersaturation(agC, brC, tempC),
    bins: new Float64Array(bins),
    nucleationRate: 0,
    meanRadius: mean,
    stdRadius: std,
    totalGrains: total,
    iodideMolPct: iodMolPctFinal,
    phase: 'digest',
  });

  return { params, timepoints, bins: { edges, counts: bins } };
}

// ─── Distribution statistics using actual r per bin ───────────────────────────
// Uses binR3 to compute the true mean and std (not the bin-centre approximation).

function distributionStats(bins: Float64Array, binR3: Float64Array) {
  let total = 0, sum = 0, sum2 = 0;
  for (let i = 0; i < N_BINS; i++) {
    if (bins[i] < 1) continue;
    const r = Math.cbrt(binR3[i] / bins[i]);
    total += bins[i];
    sum   += bins[i] * r;
    sum2  += bins[i] * r * r;
  }
  if (total === 0) return { mean: 0, std: 0, total: 0 };
  const mean = sum / total;
  const variance = sum2 / total - mean * mean;
  return { mean, std: Math.sqrt(Math.max(0, variance)), total };
}

// ─── Crystal habit prediction ──────────────────────────────────────────────────

export function predictCrystalHabit(finalPAg: number, _iodideMolPct: number): string {
  // Habit is controlled by pAg at completion (Mason 1966 Ch.3; Haist 1979 Vol.1 p.188).
  // T-grain formation requires twin-plane nucleation under controlled double-jet at high
  // pAg — conditions not yet modelled here.  Low-level iodide (< ~5 mol%) is a surface
  // dopant that does NOT convert an emulsion to tabular habit.
  if (finalPAg > 8.5) return 'Cubic {100}';
  if (finalPAg < 6.5) return 'Octahedral {111}';
  return 'Mixed / Intermediate';
}

// ─── Relative speed / ISO estimate ────────────────────────────────────────────
//
// Pre-sensitisation intrinsic speed estimate.
//
// Physical basis (Mason 1966, "Photographic Processing Chemistry", pp. 24–31):
//   Grain sensitivity scales with mean grain cross-sectional area (d²) for
//   surface-latent-image emulsions.  Combined with grain number density N:
//
//     Speed ∝ N^(1/3) × d^2        [Mason 1966 eq. 2.14; Haist 1979 Vol.1 p.188]
//
// Calibration against community-validated recipe (AJ-12 with ammonia):
//   Mean post-digest diameter ≈ 500 nm, grain count ≈ 2×10¹⁴ /L
//   Community-reported ISO range: 3–40 (no chemical sensitisation)
//   Geometric-mean calibration point: ISO_ref ≈ 10
//
// WARNING: This is a pre-sensitisation estimate only.
// Chemical sensitisation (S + Au) typically adds 1–2 stops.
// Spectral sensitisation adds further speed in target wavelengths.
// Stated ISO accuracy: ±1.5 stops (factor ~3×).
//
// Returns { value, low, high } in ISO units (scene-speed ASA/ISO equivalent).

export function estimateISO(
  meanRadiusM: number,
  totalGrainsPerL: number,
): { value: number; low: number; high: number } {
  if (meanRadiusM <= 0 || totalGrainsPerL <= 0) return { value: 0, low: 0, high: 0 };

  const d      = meanRadiusM * 2;
  const d_ref  = 722e-9;          // AJ-12 ammonia post-digest mean diameter (sim-calibrated)
  const N_ref  = 7.5e13;          // AJ-12 ammonia post-digest grain count /L (sim-calibrated)
  const ISO_ref = 4;              // AJ-12 pre-chem community reference (Ron Mowrey / Photrio)

  // Speed ∝ N^(1/3) × d^2  (Mason 1966 / Haist 1979)
  const value = ISO_ref * Math.pow(d / d_ref, 2) * Math.pow(totalGrainsPerL / N_ref, 1 / 3);

  // ±1.5 stop uncertainty band (factor ~3×)
  return { value, low: value / 3, high: value * 3 };
}

export function estimateRelativeSpeed(meanRadiusM: number, habitLabel: string): number {
  const area = Math.PI * meanRadiusM * meanRadiusM;
  const habitFactor = habitLabel.includes('Tabular') ? 2.5 : 1.0;
  return area * habitFactor * 1e12;
}
