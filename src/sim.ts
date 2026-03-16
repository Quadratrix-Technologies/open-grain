// sim.ts — Silver halide precipitation physics engine
// SRAD (Single Run Ammonia Digest) population balance model
//
// References:
//   US5422825A — supersaturation control during precipitation
//   US5248577A — adaptive halide feed rate control
//   Lifshitz-Slyozov-Wagner (LSW) theory — Ostwald ripening
//   Classical Nucleation Theory (CNT) — Kashchiev (2000)
//   Margolis & Gutoff (1974) AIChE J. 20(3) — AgBr precipitation model
//   Sugimoto (2000) Colloids Surf. A — AgBr nucleation in gelatin: Sm=3.61, r*=1.97nm
//   Sugimoto & Shiba (1999) J. Phys. Chem. B — AgBr interfacial energy

// ─── Physical constants ───────────────────────────────────────────────────────

const BOLTZMANN  = 1.38065e-23; // J/K
const AVOGADRO   = 6.02214e23;  // mol⁻¹

// ─── AgBr material constants (literature values) ──────────────────────────────

// Ksp for AgBr at 40°C (mol²/L²)
const KSP_AGBR_40C = 6.3e-13;

// Molar volume of AgBr — from density 6.473 g/cm³, M = 187.77 g/mol
// Vm = 187.77 / 6.473 = 29.01 cm³/mol = 2.901e-5 m³/mol
// Source: Crystran / PubChem
const VM_AGBR = 2.901e-5; // m³/mol

// Solid-liquid interfacial energy of AgBr in gelatin/aqueous solution
// 0.14 J/m²: Kodak patent US5448160; consistent with Sugimoto (2000) gelatin data
// (clean-interface value from Sugimoto & Shiba 1999 is 0.68 J/m² — not used here)
const GAMMA_AGBR = 0.14; // J/m²

// Heterogeneous nucleation correction factor for gelatin system.
// Gelatin acts as a nucleation template, reducing the CNT barrier by f(θ):
//   f(θ) = (2 + cosθ)(1 − cosθ)² / 4
// Calibrated so nucleation onset matches Sm = 3.61 (Sugimoto 2000 gelatin expts).
// f ≈ 0.01 corresponds to contact angle θ ≈ 12° — favorable nucleation on gelatin.
const CNT_HET_FACTOR = 0.01;

// CNT pre-exponential factor A (m⁻³ s⁻¹).
// Theoretical homogeneous nucleation: 10²⁵–10³⁰ m⁻³/s (Margolis & Gutoff 1974).
// For gelatin-mediated heterogeneous nucleation the effective A is much lower.
// Calibrated so that nucleation rate × plume volume gives ~10¹³–10¹⁴ nuclei/L
// across a 10-min SRAD run (target grain counts from community recipes).
// CNT_HET_FACTOR already reduces the thermodynamic barrier; A is the kinetic prefactor.
const CNT_A = 1e20; // m⁻³/s  (fitted; Margolis & Gutoff note insensitivity to A)

// Critical supersaturation for AgBr nucleation onset in gelatin (Sugimoto 2000).
// Below this threshold nucleation rate is set to zero.
const S_NUC_THRESHOLD = 3.61;

// ─── Growth / ripening constants ──────────────────────────────────────────────

// Ostwald ripening (LSW) rate constant — scales with ammonia concentration
const K_LSW_BASE = 5e-28; // m³/s


// ─── Types ────────────────────────────────────────────────────────────────────

export interface SimParams {
  // Reactor
  reactorVolume: number;      // L
  temperature: number;        // °C
  gelatinConc: number;        // g/L (affects growth kinetics)

  // Silver nitrate
  agno3Conc: number;          // mol/L initial

  // Halide jet (KBr)
  halidConc: number;          // mol/L
  jetFlowRate: number;        // mL/min (constant for now)
  jetDuration: number;        // seconds

  // Silver jet — double-jet mode (set agJetFlowRate > 0 to enable)
  // In double-jet: agno3Conc should be 0 (reactor starts with gelatin only)
  agJetConc: number;          // mol/L — AgNO₃ jet concentration
  agJetFlowRate: number;      // mL/min — 0 = SRAD single-jet, >0 = double-jet

  // Ammonia digest
  ammoniaConc: number;        // mol/L (when added)
  digestDuration: number;     // seconds

  // Reactor mixing
  mixingTime: number;         // seconds — jet plume dispersal time (stirring efficiency)
                              // ~1s = vigorous, ~10s = gentle, ~30s = poor mixing
}

export interface SimState {
  time: number;               // seconds
  agConc: number;             // [Ag+] mol/L
  brConc: number;             // [Br-] mol/L
  volume: number;             // L (grows as jet adds solution)
  pAg: number;                // -log10([Ag+])
  supersaturation: number;    // S = [Ag+][Br-] / Ksp
  bins: Float64Array;         // grain count per radius bin
  nucleationRate: number;     // nuclei/s
  meanRadius: number;         // m
  stdRadius: number;          // m
  totalGrains: number;
  phase: 'precipitation' | 'digest';
}

export interface SimResult {
  params: SimParams;
  timepoints: SimState[];
  bins: { edges: Float64Array; counts: Float64Array }; // final distribution
}

// ─── Bin geometry ─────────────────────────────────────────────────────────────

const N_BINS = 80;
const R_MIN = 1e-9;   // 1 nm
const R_MAX = 2e-6;   // 2 µm — covers typical AgBr grain range

// Log-spaced bin edges
function makeBinEdges(): Float64Array {
  const edges = new Float64Array(N_BINS + 1);
  const logMin = Math.log10(R_MIN);
  const logMax = Math.log10(R_MAX);
  for (let i = 0; i <= N_BINS; i++) {
    edges[i] = Math.pow(10, logMin + (i / N_BINS) * (logMax - logMin));
  }
  return edges;
}

function binCenters(edges: Float64Array): Float64Array {
  const c = new Float64Array(N_BINS);
  for (let i = 0; i < N_BINS; i++) {
    c[i] = Math.sqrt(edges[i] * edges[i + 1]); // geometric mean
  }
  return c;
}

// ─── Physics helpers ──────────────────────────────────────────────────────────

function ksp(tempC: number): number {
  // Approximate temperature dependence of Ksp(AgBr)
  // Ksp increases with temperature
  const dT = tempC - 40;
  return KSP_AGBR_40C * Math.pow(10, 0.012 * dT);
}

function supersaturation(agC: number, brC: number, tempC: number): number {
  if (agC <= 0 || brC <= 0) return 0;
  return (agC * brC) / ksp(tempC);
}

function pAg(agC: number): number {
  if (agC <= 0) return 99;
  return -Math.log10(agC);
}

// Nucleation rate — Classical Nucleation Theory (CNT)
// J = A × exp(−B_eff / (ln S)²)   [nuclei / L / s]
//
// B_hom = (16π/3) × (γ/kT)³ × v²  — thermodynamic barrier, computed at runtime
// B_eff = B_hom × CNT_HET_FACTOR   — corrected for gelatin heterogeneous nucleation
//
// At very high S (jet plume): exp term → 1, rate saturates at CNT_A — no blow-up.
// At S < Sm = 3.61: rate forced to zero (Sugimoto 2000 threshold).
function nucleationRate(S: number, tempC: number): number {
  if (S <= S_NUC_THRESHOLD) return 0;
  const TK   = tempC + 273.15;
  const kT   = BOLTZMANN * TK;
  const v    = VM_AGBR / AVOGADRO;           // m³ per formula unit
  const B_hom = (16 * Math.PI / 3) * Math.pow(GAMMA_AGBR / kT, 3) * v * v;
  const B_eff = B_hom * CNT_HET_FACTOR;
  const lnS  = Math.log(S);
  const rateM3 = CNT_A * Math.exp(-B_eff / (lnS * lnS)); // nuclei / m³ / s
  return rateM3 * 1e-3; // → nuclei / L / s
}

// ─── Simulator ────────────────────────────────────────────────────────────────

export function runSimulation(params: SimParams): SimResult {
  const edges = makeBinEdges();
  const centers = binCenters(edges);

  const dt = 0.5; // seconds per step
  const totalTime = params.jetDuration + params.digestDuration;
  const nSteps = Math.ceil(totalTime / dt);

  // State
  let agC = params.agno3Conc;     // mol/L
  let brC = 0.0;
  let vol = params.reactorVolume; // L
  const bins = new Float64Array(N_BINS); // grain count per bin
  const tempC = params.temperature;

  const timepoints: SimState[] = [];
  const recordEvery = Math.max(1, Math.floor(nSteps / 400)); // ~400 frames max

  for (let step = 0; step < nSteps; step++) {
    const t = step * dt;
    const phase: 'precipitation' | 'digest' = t < params.jetDuration ? 'precipitation' : 'digest';

    // ── Halide jet: deliver this step's Br⁻ ──────────────────────────────────
    // Nucleation happens in the jet plume BEFORE Br⁻ mixes into the bulk.
    // So we:
    //   1. Dilute existing bulk species for the volume increase
    //   2. Nucleate using the plume S (jet conc × bulk Ag⁺) — consuming Ag⁺ only
    //   3. Add the leftover (un-nucleated) Br⁻ to the bulk
    // This prevents nucleation from draining bulk brC and starving growth.

    if (phase === 'precipitation') {
      // ── Volume additions ─────────────────────────────────────────────────
      const jetFlowL   = (params.jetFlowRate   / 1000 / 60) * dt; // Br⁻ jet
      const agJetFlowL = (params.agJetFlowRate / 1000 / 60) * dt; // Ag⁺ jet (0 in SRAD)
      const molsBrThisStep = jetFlowL   * params.halidConc;
      const molsAgJet      = agJetFlowL * params.agJetConc;

      // Dilute existing bulk species for the combined volume increase
      const newVol = vol + jetFlowL + agJetFlowL;
      agC = (agC * vol) / newVol;
      brC = (brC * vol) / newVol;
      vol = newVol;

      // Shared bin geometry for both plume nucleation paths
      const r0              = centers[0];
      const molesPerNucleus = ((4 / 3) * Math.PI * r0 * r0 * r0) / VM_AGBR;
      const mixingFactor    = dt / params.mixingTime;

      // ── Br⁻ plume nucleation (both SRAD and double-jet) ──────────────────
      // Jet Br⁻ meets bulk Ag⁺. Nuclei consume Ag⁺ from bulk; Br⁻ came from jet.
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
        molsBrNucleated  = newNuclei_br * molesPerNucleus;
        agC              = Math.max(0, agC - molsBrNucleated / vol);
      }

      // ── Ag⁺ plume nucleation (double-jet only) ───────────────────────────
      // Jet Ag⁺ meets bulk Br⁻. Nuclei consume Br⁻ from bulk; Ag⁺ came from jet.
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
          molsAgNucleated   = newNuclei_ag * molesPerNucleus;
          brC               = Math.max(0, brC - molsAgNucleated / vol);
        }
      }

      // Leftover from both jets mixes into the bulk
      brC += Math.max(0, molsBrThisStep - molsBrNucleated) / vol;
      agC += Math.max(0, molsAgJet      - molsAgNucleated) / vol;
    }

    const S_bulk = supersaturation(agC, brC, tempC);

    // ── Equilibrium precipitation growth (bulk) ───────────────────────────
    // Solve (agC − ΔC)(brC − ΔC) = Ksp for ΔC, then distribute the
    // resulting AgBr across all grains proportional to their surface area.
    // This is physically exact and avoids the blow-up of dr-based models.
    if (S_bulk > 1.0) {
      const kspT = ksp(tempC);
      // Numerically stable quadratic: ΔC = 2(agC·brC−Ksp) / [(agC+brC)+√((agC−brC)²+4Ksp)]
      // The standard form [(agC+brC)−√disc]/2 suffers catastrophic cancellation when
      // agC << brC (or vice versa) — this alternative form avoids that.
      const diff = agC - brC;
      const disc = Math.sqrt(diff * diff + 4 * kspT);
      const deltaC = 2 * (agC * brC - kspT) / (agC + brC + disc);
      const maxDeltaC = Math.min(deltaC, agC, brC);
      const totalMolsToPrecip = maxDeltaC * vol;

      // Weight grains by surface area (r²) — diffusion-limited growth
      let totalArea = 0;
      for (let i = 0; i < N_BINS; i++) {
        if (bins[i] >= 1) totalArea += bins[i] * centers[i] * centers[i];
      }

      if (totalArea > 0 && totalMolsToPrecip > 0) {
        let molsConsumed = 0;
        // Snapshot prevents double-counting: grains moving to a higher bin
        // during a forward sweep would otherwise be processed again.
        const snap = new Float64Array(bins);

        for (let i = 0; i < N_BINS; i++) {
          if (snap[i] < 1) continue;
          const r = centers[i];

          // Moles deposited on this bin's grains (based on original distribution)
          const molsThisBin = totalMolsToPrecip * (snap[i] * r * r) / totalArea;

          // Volume added per grain → new radius
          const dVolPerGrain = (molsThisBin / snap[i]) * VM_AGBR;
          const r3New = r * r * r + dVolPerGrain * (3 / (4 * Math.PI));
          const newR = Math.cbrt(r3New);

          // Consume moles regardless of whether grain moves bins
          molsConsumed += molsThisBin;

          if (newR >= R_MAX) {
            bins[i] -= snap[i]; // remove these grains (consumed into oversized territory)
            continue;
          }

          // Find destination bin by edge comparison
          let destBin = i;
          while (destBin < N_BINS - 1 && edges[destBin + 1] < newR) destBin++;

          if (destBin !== i) {
            // Subtract only the snapshot count — other grains moved into bins[i]
            // during this sweep must not be erased.
            bins[i] -= snap[i];
            bins[destBin] += snap[i];
          }
        }

        agC = Math.max(0, agC - molsConsumed / vol);
        brC = Math.max(0, brC - molsConsumed / vol);
      }
    }

    // ── Ostwald ripening (digest phase) ───────────────────────────────────
    if (phase === 'digest') {
      const ammoniaFactor = 1 + 10 * params.ammoniaConc; // ammonia accelerates ripening
      const kLSW = K_LSW_BASE * ammoniaFactor;

      // LSW: critical radius r* = 2γVm / (RT ln S), grains below dissolve, above grow
      // Simplified: compute mean radius as r*, dissolve below, grow above
      const { mean } = distributionStats(bins, centers);
      const rCrit = mean > 0 ? mean : centers[0];

      // Snapshot prevents double-counting when grains move between bins
      const snapLSW = new Float64Array(bins);

      for (let i = 0; i < N_BINS; i++) {
        if (snapLSW[i] < 1) continue;
        const r = centers[i];
        // LSW growth rate: dr³/dt = kLSW * (1/rCrit - 1/r)
        const dr3dt = kLSW * (1 / rCrit - 1 / r);
        const r3new = r * r * r + dr3dt * dt;

        if (r3new <= 0) {
          // Grain dissolved — release Ag+ and Br- back to solution
          const volGrain = (4 / 3) * Math.PI * r * r * r;
          const molsReleased = snapLSW[i] * volGrain / VM_AGBR;
          agC += molsReleased / vol;
          brC += molsReleased / vol;
          bins[i] -= snapLSW[i];
        } else {
          const rNew = Math.cbrt(r3new);
          if (rNew >= R_MAX) continue;
          let destBin = i;
          while (destBin < N_BINS - 1 && edges[destBin + 1] < rNew) destBin++;
          while (destBin > 0 && edges[destBin] > rNew) destBin--;
          if (destBin !== i) {
            bins[i] -= snapLSW[i];
            bins[destBin] += snapLSW[i];
          }
        }
      }
    }

    // ── Dissolution equilibration ─────────────────────────────────────────
    // If S < 1 and grains are present, AgBr dissolves back until S = 1.
    // The amount is always tiny (< Ksp^0.5 ~ 1e-6 mol/L) so we update
    // concentrations only; the effect on grain volumes is negligible.
    const S_check = supersaturation(agC, brC, tempC);
    if (S_check < 1.0) {
      const kspT = ksp(tempC);
      let hasGrains = false;
      for (let i = 0; i < N_BINS; i++) if (bins[i] >= 1) { hasGrains = true; break; }
      if (hasGrains) {
        const diff_d = agC - brC;
        const dissolveC = 2 * (kspT - agC * brC) /
          (agC + brC + Math.sqrt(diff_d * diff_d + 4 * kspT));
        agC += dissolveC;
        brC += dissolveC;
      }
    }

    // ── Record state ───────────────────────────────────────────────────────
    if (step % recordEvery === 0) {
      const { mean, std, total } = distributionStats(bins, centers);
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
        phase,
      });
    }
  }

  const { mean, std, total } = distributionStats(bins, centers);
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
    phase: 'digest',
  });

  return {
    params,
    timepoints,
    bins: { edges, counts: bins },
  };
}

function distributionStats(bins: Float64Array, centers: Float64Array) {
  let total = 0;
  let sum = 0;
  let sum2 = 0;
  for (let i = 0; i < N_BINS; i++) {
    total += bins[i];
    sum += bins[i] * centers[i];
    sum2 += bins[i] * centers[i] * centers[i];
  }
  if (total === 0) return { mean: 0, std: 0, total: 0 };
  const mean = sum / total;
  const variance = sum2 / total - mean * mean;
  return { mean, std: Math.sqrt(Math.max(0, variance)), total };
}

// ─── Crystal habit prediction ─────────────────────────────────────────────────

export function predictCrystalHabit(finalPAg: number, iodideMolPct: number): string {
  // Based on pAg regime and iodide content:
  // High pAg (Ag-deficient) + low iodide → cubic (100 faces)
  // Low pAg (Ag-excess) + low iodide → octahedral (111 faces)
  // Low pAg + >2 mol% iodide → tabular (T-grain, 111 faces dominant)
  if (iodideMolPct >= 2) return 'Tabular (T-grain)';
  if (finalPAg > 8.5) return 'Cubic {100}';
  if (finalPAg < 6.5) return 'Octahedral {111}';
  return 'Mixed / Intermediate';
}

// Relative speed estimate: proportional to mean grain cross-section area
export function estimateRelativeSpeed(meanRadiusM: number, habitLabel: string): number {
  const area = Math.PI * meanRadiusM * meanRadiusM;
  const habitFactor = habitLabel.includes('Tabular') ? 2.5 : 1.0; // T-grains more efficient
  return area * habitFactor * 1e12; // normalised to readable units
}
