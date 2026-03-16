// main.ts — Open Grain UI
// Silver halide precipitation simulator
// Honoring Ron Mowrey (PE), Grant Haist, and L.F.A. Mason

import { runSimulation, predictCrystalHabit, estimateRelativeSpeed } from './sim';
import type { SimParams, SimState } from './sim';
import {
  drawChemistryPlot,
  drawGrainHistogram,
  drawGrowthCurve,
  drawNucleationPlot,
  THEME,
} from './plots';

// ─── Layout ───────────────────────────────────────────────────────────────────

const CSS = `
* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  background: ${THEME.bg};
  color: ${THEME.yellow};
  font-family: ${THEME.font};
  font-size: 12px;
  height: 100vh;
  overflow: hidden;
  display: flex;
  flex-direction: column;
}

#titlebar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 6px 12px;
  border-bottom: 1px solid ${THEME.border};
  background: ${THEME.bgPanel};
  flex-shrink: 0;
}

#titlebar h1 {
  font-size: 13px;
  font-weight: 600;
  letter-spacing: 0.08em;
  color: ${THEME.yellow};
}

#titlebar .subtitle {
  color: ${THEME.yellowDim};
  font-size: 10px;
}

#titlebar .controls {
  display: flex;
  gap: 8px;
  align-items: center;
}

button {
  background: transparent;
  border: 1px solid ${THEME.border};
  color: ${THEME.yellow};
  font-family: ${THEME.font};
  font-size: 11px;
  padding: 3px 10px;
  cursor: pointer;
  letter-spacing: 0.05em;
}

button:hover {
  border-color: ${THEME.yellow};
  background: #1a1800;
}

button.primary {
  border-color: ${THEME.yellow};
  color: ${THEME.bg};
  background: ${THEME.yellow};
}

button.primary:hover {
  background: #ffe030;
}

#workspace {
  display: flex;
  flex: 1;
  overflow: hidden;
}

#sidebar {
  width: 240px;
  flex-shrink: 0;
  border-right: 1px solid ${THEME.border};
  background: ${THEME.bgPanel};
  overflow-y: auto;
  padding: 0;
}

.param-group {
  border-bottom: 1px solid ${THEME.border};
  padding: 8px 10px;
}

.param-group-title {
  color: ${THEME.yellowDim};
  font-size: 10px;
  letter-spacing: 0.1em;
  text-transform: uppercase;
  margin-bottom: 8px;
}

.param-row {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-bottom: 5px;
  gap: 6px;
}

.param-row label {
  color: ${THEME.yellowDim};
  font-size: 11px;
  white-space: nowrap;
  flex-shrink: 0;
  min-width: 90px;
}

.param-row input[type="number"] {
  background: ${THEME.bg};
  border: 1px solid ${THEME.border};
  color: ${THEME.yellow};
  font-family: ${THEME.font};
  font-size: 11px;
  padding: 2px 5px;
  width: 80px;
  text-align: right;
}

.param-row input[type="number"]:focus {
  outline: none;
  border-color: ${THEME.yellow};
}

.param-unit {
  color: ${THEME.dimText};
  font-size: 10px;
  width: 28px;
  flex-shrink: 0;
}

#stats-panel {
  padding: 8px 10px;
}

.stats-title {
  color: ${THEME.yellowDim};
  font-size: 10px;
  letter-spacing: 0.1em;
  text-transform: uppercase;
  margin-bottom: 8px;
}

.stat-row {
  display: flex;
  justify-content: space-between;
  margin-bottom: 4px;
  font-size: 11px;
}

.stat-key { color: ${THEME.yellowDim}; }
.stat-val { color: ${THEME.yellow}; }

#plots {
  flex: 1;
  display: grid;
  grid-template-columns: 1fr 1fr;
  grid-template-rows: 1fr 1fr;
  gap: 1px;
  background: ${THEME.border};
  overflow: hidden;
}

.plot-pane {
  background: ${THEME.bgPanel};
  display: flex;
  flex-direction: column;
  overflow: hidden;
  position: relative;
}

.plot-title {
  font-size: 10px;
  color: ${THEME.yellowDim};
  letter-spacing: 0.08em;
  padding: 5px 8px 0 8px;
  text-transform: uppercase;
  flex-shrink: 0;
}

.plot-pane svg {
  flex: 1;
  display: block;
}

#status-bar {
  border-top: 1px solid ${THEME.border};
  background: ${THEME.bgPanel};
  padding: 3px 12px;
  font-size: 10px;
  color: ${THEME.yellowDim};
  flex-shrink: 0;
  display: flex;
  gap: 24px;
}

#status-bar span { color: ${THEME.yellow}; }

::-webkit-scrollbar { width: 6px; }
::-webkit-scrollbar-track { background: ${THEME.bg}; }
::-webkit-scrollbar-thumb { background: ${THEME.border}; }

#preset-group {
  border-bottom: 1px solid ${THEME.border};
  padding: 8px 10px;
}

#preset-group .param-group-title { margin-bottom: 6px; }

#preset-select {
  width: 100%;
  background: ${THEME.bg};
  border: 1px solid ${THEME.border};
  color: ${THEME.yellow};
  font-family: ${THEME.font};
  font-size: 11px;
  padding: 3px 5px;
  margin-bottom: 5px;
  cursor: pointer;
}

#preset-select:focus { outline: none; border-color: ${THEME.yellow}; }

#preset-note {
  color: ${THEME.dimText};
  font-size: 10px;
  line-height: 1.4;
  min-height: 14px;
}

#preset-source {
  color: ${THEME.dimText};
  font-size: 9px;
  margin-top: 3px;
  font-style: italic;
}
`;

// ─── Known emulsion presets ───────────────────────────────────────────────────
// Parameters derived from published community recipes.
// Note: several original recipes add silver into halide (reverse of this sim's
// halide-into-silver convention). Parameters are transposed to match our setup.

interface Preset {
  label: string;
  source: string;
  note?: string;
  params: SimParams;
}

const PRESETS: Preset[] = [
  {
    label: 'PE SRAD / Kodak AJ-12 (ammonia)',
    source: 'Ron Mowrey (Photo Engineer), Photrio 2006 — Kodak AJ-12',
    note: 'Original adds silver to halide; parameters transposed for this sim.',
    params: {
      reactorVolume: 0.5,
      temperature: 45,
      gelatinConc: 30,
      agno3Conc: 1.53,
      halidConc: 1.11,
      jetFlowRate: 100,
      jetDuration: 600,
      agJetConc: 0,
      agJetFlowRate: 0,
      ammoniaConc: 0.3,
      digestDuration: 1800,
      mixingTime: 5,
    },
  },
  {
    label: 'AJ-12 no-ammonia',
    source: 'Kodak AJ-12 via UnblinkingEye.com',
    note: 'No ammonia digest; shorter ripening. Fine-grain result.',
    params: {
      reactorVolume: 0.36,
      temperature: 55,
      gelatinConc: 28,
      agno3Conc: 0.59,
      halidConc: 0.75,
      jetFlowRate: 40,
      jetDuration: 600,
      agJetConc: 0,
      agJetFlowRate: 0,
      ammoniaConc: 0,
      digestDuration: 1500,
      mixingTime: 5,
    },
  },
  {
    label: 'Ludwik 2024 (small batch)',
    source: 'Photrio 2024 — ludwik',
    note: 'Small-scale modern recipe; long thermal digest, ~ISO 30.',
    params: {
      reactorVolume: 0.1,
      temperature: 55,
      gelatinConc: 25,
      agno3Conc: 0.59,
      halidConc: 1.01,
      jetFlowRate: 5,
      jetDuration: 1200,
      agJetConc: 0,
      agJetFlowRate: 0,
      ammoniaConc: 0,
      digestDuration: 5400,
      mixingTime: 5,
    },
  },
  {
    label: 'Fulvio baseline (no-ammonia)',
    source: 'Photrio 2006 — Fulvio / Photo Engineer modifications',
    note: 'Minimal starter formula. ~ISO 3.',
    params: {
      reactorVolume: 0.125,
      temperature: 50,
      gelatinConc: 160,
      agno3Conc: 0.94,
      halidConc: 1.07,
      jetFlowRate: 12.5,
      jetDuration: 600,
      agJetConc: 0,
      agJetFlowRate: 0,
      ammoniaConc: 0,
      digestDuration: 600,
      mixingTime: 5,
    },
  },
  {
    label: 'LeoniD core-shell (double-jet)',
    source: 'Photrio 2023 — LeoniD',
    note: 'Double-jet: simultaneous Ag\u207a and Br\u207b jets into gelatin. KI shell added separately. 45\u00b0C, 30 min digest.',
    params: {
      reactorVolume: 0.05,
      temperature: 45,
      gelatinConc: 40,
      agno3Conc: 0,         // double-jet: no pre-loaded silver
      halidConc: 5.51,      // KBr jet: 20.3 g in 31 mL
      jetFlowRate: 6.0,     // Br\u207b jet flow (B solution)
      jetDuration: 310,     // ~31 mL at 6 mL/min
      agJetConc: 1.96,      // AgNO\u2083 jet: 20 g in 60 mL + NH\u2084OH
      agJetFlowRate: 3.1,   // Ag\u207a jet flow (C solution)
      ammoniaConc: 0.3,
      digestDuration: 1800,
      mixingTime: 2,        // 2000 rpm overhead stirring
    },
  },
];

// ─── Default parameters ───────────────────────────────────────────────────────

const DEFAULT_PARAMS: SimParams = {
  reactorVolume: 1.0,
  temperature: 40,
  gelatinConc: 40,
  agno3Conc: 0.5,
  halidConc: 1.0,
  jetFlowRate: 5.0,
  jetDuration: 300,
  agJetConc: 0,
  agJetFlowRate: 0,
  ammoniaConc: 0.1,
  digestDuration: 600,
  mixingTime: 5,
};

// ─── App state ────────────────────────────────────────────────────────────────

let currentTimepoints: SimState[] = [];
let currentBins: Float64Array = new Float64Array(80);
let currentBinEdges: Float64Array = new Float64Array(81);

// ─── Build UI ─────────────────────────────────────────────────────────────────

function buildUI() {
  const style = document.createElement('style');
  style.textContent = CSS;
  document.head.appendChild(style);

  const app = document.getElementById('app')!;
  app.style.display = 'contents';

  // Title bar
  const titlebar = el('div', 'titlebar');
  titlebar.innerHTML = `
    <div>
      <h1>OPEN GRAIN</h1>
      <div class="subtitle">silver halide precipitation simulator</div>
    </div>
    <div class="controls">
      <button id="btn-reset">RESET</button>
      <button id="btn-run" class="primary">RUN</button>
    </div>
  `;
  app.appendChild(titlebar);

  // Workspace
  const workspace = el('div', 'workspace');

  // Sidebar
  const sidebar = el('div', 'sidebar');
  sidebar.appendChild(buildPresetGroup());
  sidebar.appendChild(buildParamGroup('REACTOR', [
    { key: 'reactorVolume', label: 'volume', unit: 'L', step: 0.1, min: 0.1, max: 10 },
    { key: 'temperature', label: 'temp', unit: '°C', step: 1, min: 20, max: 70 },
    { key: 'gelatinConc', label: 'gelatin', unit: 'g/L', step: 5, min: 5, max: 200 },
    { key: 'mixingTime', label: 'mix time \u03c4', unit: 's', step: 1, min: 1, max: 60 },
  ]));
  sidebar.appendChild(buildParamGroup('SILVER NITRATE', [
    { key: 'agno3Conc', label: '[AgNO\u2083]', unit: 'M', step: 0.05, min: 0, max: 2 },
  ]));
  sidebar.appendChild(buildParamGroup('SILVER JET \u2014 double-jet (0\u2009=\u2009SRAD)', [
    { key: 'agJetConc', label: '[AgNO\u2083] jet', unit: 'M', step: 0.05, min: 0, max: 4 },
    { key: 'agJetFlowRate', label: 'Ag flow rate', unit: 'mL/m', step: 0.1, min: 0, max: 20 },
  ]));
  sidebar.appendChild(buildParamGroup('HALIDE JET', [
    { key: 'halidConc', label: '[KBr]', unit: 'M', step: 0.1, min: 0.1, max: 4 },
    { key: 'jetFlowRate', label: 'flow rate', unit: 'mL/m', step: 0.5, min: 0.5, max: 50 },
    { key: 'jetDuration', label: 'duration', unit: 's', step: 10, min: 30, max: 3600 },
  ]));
  sidebar.appendChild(buildParamGroup('AMMONIA DIGEST', [
    { key: 'ammoniaConc', label: '[NH\u2083]', unit: 'M', step: 0.01, min: 0, max: 1 },
    { key: 'digestDuration', label: 'hold time', unit: 's', step: 30, min: 0, max: 3600 },
  ]));

  const statsPanel = el('div', 'stats-panel');
  statsPanel.innerHTML = `<div class="stats-title">RESULTS</div><div id="stats-content"></div>`;
  sidebar.appendChild(statsPanel);

  workspace.appendChild(sidebar);

  // Plots grid
  const plots = el('div', 'plots');
  plots.appendChild(plotPane('chemistry', 'pAg / Supersaturation'));
  plots.appendChild(plotPane('growth', 'Mean Grain Diameter'));
  plots.appendChild(plotPane('histogram', 'Grain Size Distribution'));
  plots.appendChild(plotPane('nucleation', 'Nucleation Rate'));
  workspace.appendChild(plots);

  app.appendChild(workspace);

  // Status bar
  const statusBar = el('div', 'status-bar');
  statusBar.innerHTML = `
    <div>status: <span id="s-status">idle</span></div>
    <div>grains: <span id="s-grains">—</span></div>
    <div>mean \u00d8: <span id="s-diam">—</span></div>
    <div>final pAg: <span id="s-pag">—</span></div>
  `;
  app.appendChild(statusBar);

  document.getElementById('btn-run')!.addEventListener('click', runSim);
  document.getElementById('btn-reset')!.addEventListener('click', resetSim);
}

function el(tag: string, id: string): HTMLElement {
  const e = document.createElement(tag);
  e.id = id;
  return e;
}

type ParamDef = {
  key: keyof SimParams;
  label: string;
  unit: string;
  step: number;
  min: number;
  max: number;
};

function buildPresetGroup(): HTMLElement {
  const group = document.createElement('div');
  group.id = 'preset-group';
  group.innerHTML = `<div class="param-group-title">KNOWN EMULSIONS</div>`;

  const select = document.createElement('select');
  select.id = 'preset-select';

  const blank = document.createElement('option');
  blank.value = '';
  blank.textContent = '— select preset —';
  select.appendChild(blank);

  PRESETS.forEach((p, i) => {
    const opt = document.createElement('option');
    opt.value = String(i);
    opt.textContent = p.label;
    select.appendChild(opt);
  });

  const note = document.createElement('div');
  note.id = 'preset-note';

  const source = document.createElement('div');
  source.id = 'preset-source';

  select.addEventListener('change', () => {
    const idx = parseInt(select.value);
    if (isNaN(idx)) {
      note.textContent = '';
      source.textContent = '';
      return;
    }
    const preset = PRESETS[idx];
    applyPreset(preset.params);
    note.textContent = preset.note ?? '';
    source.textContent = preset.source;
  });

  group.appendChild(select);
  group.appendChild(note);
  group.appendChild(source);
  return group;
}

function applyPreset(params: SimParams) {
  (Object.keys(params) as (keyof SimParams)[]).forEach(key => {
    const inp = document.getElementById(`p-${key}`) as HTMLInputElement | null;
    if (inp) inp.value = String(params[key]);
  });
}

function buildParamGroup(title: string, params: ParamDef[]): HTMLElement {
  const group = document.createElement('div');
  group.className = 'param-group';
  group.innerHTML = `<div class="param-group-title">${title}</div>`;

  for (const p of params) {
    const row = document.createElement('div');
    row.className = 'param-row';
    row.innerHTML = `
      <label for="p-${p.key}">${p.label}</label>
      <input type="number" id="p-${p.key}"
        value="${DEFAULT_PARAMS[p.key]}"
        step="${p.step}" min="${p.min}" max="${p.max}" />
      <span class="param-unit">${p.unit}</span>
    `;
    group.appendChild(row);
  }
  return group;
}

function plotPane(id: string, title: string): HTMLElement {
  const pane = document.createElement('div');
  pane.className = 'plot-pane';
  pane.innerHTML = `<div class="plot-title">${title}</div>`;
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.id = `svg-${id}`;
  pane.appendChild(svg);
  return pane;
}

// ─── Read params ──────────────────────────────────────────────────────────────

function readParams(): SimParams {
  function v(key: keyof SimParams): number {
    const inp = document.getElementById(`p-${key}`) as HTMLInputElement | null;
    return inp ? parseFloat(inp.value) : (DEFAULT_PARAMS[key] as number);
  }
  return {
    reactorVolume: v('reactorVolume'),
    temperature: v('temperature'),
    gelatinConc: v('gelatinConc'),
    agno3Conc: v('agno3Conc'),
    halidConc: v('halidConc'),
    jetFlowRate: v('jetFlowRate'),
    jetDuration: v('jetDuration'),
    agJetConc: v('agJetConc'),
    agJetFlowRate: v('agJetFlowRate'),
    ammoniaConc: v('ammoniaConc'),
    digestDuration: v('digestDuration'),
    mixingTime: v('mixingTime'),
  };
}

// ─── Run / reset ──────────────────────────────────────────────────────────────

function runSim() {
  setText('s-status', 'running...');
  setTimeout(() => {
    const result = runSimulation(readParams());
    currentTimepoints = result.timepoints;
    currentBins = result.bins.counts;
    currentBinEdges = result.bins.edges;
    renderPlots();
    updateStats(result.timepoints);
    setText('s-status', 'done');
  }, 20);
}

function resetSim() {
  currentTimepoints = [];
  currentBins = new Float64Array(80);
  renderPlots();
  updateStats([]);
  setText('s-status', 'idle');
  setText('s-grains', '—');
  setText('s-diam', '—');
  setText('s-pag', '—');
}

// ─── Render ───────────────────────────────────────────────────────────────────

function renderPlots() {
  draw('svg-chemistry', (e, w, h) => drawChemistryPlot(e, currentTimepoints, w, h));
  draw('svg-growth',    (e, w, h) => drawGrowthCurve(e, currentTimepoints, w, h));
  draw('svg-histogram', (e, w, h) => drawGrainHistogram(e, currentBins, currentBinEdges, w, h));
  draw('svg-nucleation',(e, w, h) => drawNucleationPlot(e, currentTimepoints, w, h));
}

function draw(id: string, fn: (el: SVGSVGElement, w: number, h: number) => void) {
  const svg = document.getElementById(id) as SVGSVGElement | null;
  if (!svg) return;
  const parent = svg.parentElement!;
  const titleEl = parent.querySelector('.plot-title') as HTMLElement | null;
  const titleH = titleEl?.offsetHeight ?? 20;
  const w = Math.floor(parent.clientWidth);
  const h = Math.floor(parent.clientHeight - titleH);
  if (w < 10 || h < 10) return;
  fn(svg, w, h);
}

function updateStats(timepoints: SimState[]) {
  const el = document.getElementById('stats-content')!;
  if (timepoints.length === 0) {
    el.innerHTML = `<div class="stat-row"><span class="stat-key">no data</span></div>`;
    return;
  }
  const f = timepoints[timepoints.length - 1];
  const dNm = (f.meanRadius * 2e9).toFixed(1);
  const sNm = (f.stdRadius * 2e9).toFixed(1);
  const cv = f.meanRadius > 0 ? ((f.stdRadius / f.meanRadius) * 100).toFixed(1) : '—';
  const habit = predictCrystalHabit(f.pAg, 0);
  const speed = estimateRelativeSpeed(f.meanRadius, habit).toFixed(2);

  el.innerHTML = `
    <div class="stat-row"><span class="stat-key">mean \u00d8</span><span class="stat-val">${dNm} nm</span></div>
    <div class="stat-row"><span class="stat-key">std dev</span><span class="stat-val">${sNm} nm</span></div>
    <div class="stat-row"><span class="stat-key">CV</span><span class="stat-val">${cv}%</span></div>
    <div class="stat-row"><span class="stat-key">total grains</span><span class="stat-val">${fmt(f.totalGrains)}</span></div>
    <div class="stat-row"><span class="stat-key">final pAg</span><span class="stat-val">${f.pAg.toFixed(2)}</span></div>
    <div class="stat-row"><span class="stat-key">habit</span><span class="stat-val">${habit}</span></div>
    <div class="stat-row"><span class="stat-key">rel. speed</span><span class="stat-val">${speed}</span></div>
  `;

  setText('s-grains', fmt(f.totalGrains));
  setText('s-diam', `${dNm} nm`);
  setText('s-pag', f.pAg.toFixed(2));
}

function setText(id: string, val: string) {
  const e = document.getElementById(id);
  if (e) e.textContent = val;
}

function fmt(n: number): string {
  if (n >= 1e9) return `${(n / 1e9).toFixed(2)}B`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(2)}M`;
  if (n >= 1e3) return `${(n / 1e3).toFixed(1)}K`;
  return n.toFixed(0);
}

// ─── Resize ───────────────────────────────────────────────────────────────────

let resizeTimer: ReturnType<typeof setTimeout>;
window.addEventListener('resize', () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => { if (currentTimepoints.length > 0) renderPlots(); }, 100);
});

// ─── Init ─────────────────────────────────────────────────────────────────────

buildUI();
