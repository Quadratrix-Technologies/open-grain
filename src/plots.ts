// plots.ts — D3 chart rendering for Open Grain
import * as d3 from 'd3';
import type { SimState } from './sim';

// ─── Theme ────────────────────────────────────────────────────────────────────

export const THEME = {
  bg: '#0d0d0b',
  bgPanel: '#111109',
  border: '#3a3000',
  yellow: '#F7CE00',      // Kodak Yellow
  yellowDim: '#7a6700',
  red: '#E03C31',         // Kodak Red
  redDim: '#6b1c17',
  white: '#e8e0c8',
  dimText: '#5a5030',
  font: "'IBM Plex Mono', 'Courier New', monospace",
};

const MARGIN = { top: 20, right: 20, bottom: 40, left: 55 };

// ─── Helpers ──────────────────────────────────────────────────────────────────

function clearSvg(el: SVGSVGElement) {
  while (el.firstChild) el.removeChild(el.firstChild);
}

function baseSvg(el: SVGSVGElement, width: number, height: number) {
  const svg = d3.select(el)
    .attr('width', width)
    .attr('height', height);

  svg.append('rect')
    .attr('width', width)
    .attr('height', height)
    .attr('fill', THEME.bgPanel);

  return svg;
}

function chartArea(svg: d3.Selection<SVGSVGElement, unknown, null, undefined>, width: number, height: number) {
  const w = width - MARGIN.left - MARGIN.right;
  const h = height - MARGIN.top - MARGIN.bottom;
  const g = svg.append('g').attr('transform', `translate(${MARGIN.left},${MARGIN.top})`);
  return { g, w, h };
}

function styleAxis(sel: d3.Selection<SVGGElement, unknown, null, undefined>) {
  sel.selectAll('line, path').attr('stroke', THEME.border);
  sel.selectAll('text')
    .attr('fill', THEME.yellowDim)
    .attr('font-family', THEME.font)
    .attr('font-size', '10px');
}

function axisLabel(
  g: d3.Selection<SVGGElement, unknown, null, undefined>,
  text: string,
  x: number,
  y: number,
  rotate = 0
) {
  g.append('text')
    .attr('transform', rotate ? `translate(${x},${y}) rotate(${rotate})` : `translate(${x},${y})`)
    .attr('text-anchor', 'middle')
    .attr('fill', THEME.yellowDim)
    .attr('font-family', THEME.font)
    .attr('font-size', '10px')
    .text(text);
}

// ─── pAg / Supersaturation plot ───────────────────────────────────────────────

export function drawChemistryPlot(
  el: SVGSVGElement,
  timepoints: SimState[],
  width: number,
  height: number
) {
  clearSvg(el);
  if (timepoints.length < 2) return;

  const svg = baseSvg(el, width, height);
  const { g, w, h } = chartArea(svg, width, height);

  const times = timepoints.map(d => d.time);
  const xScale = d3.scaleLinear().domain([0, d3.max(times)!]).range([0, w]);

  // pAg axis (left)
  const pAgVals = timepoints.map(d => d.pAg).filter(v => v < 20);
  const pAgScale = d3.scaleLinear()
    .domain([d3.max(pAgVals)! + 0.5, Math.max(0, d3.min(pAgVals)! - 0.5)])
    .range([0, h]);

  // Supersaturation axis (right) — clip at 95th percentile so initial spike
  // doesn't flatten the rest of the curve
  const sVals = timepoints.map(d => d.supersaturation);
  const sortedS = [...sVals].sort((a, b) => a - b);
  const sP95 = sortedS[Math.floor(sortedS.length * 0.95)] ?? 1;
  const sScale = d3.scaleLinear()
    .domain([0, sP95 * 1.1])
    .range([h, 0])
    .clamp(true);

  // Phase divider
  const digestStart = timepoints.find(d => d.phase === 'digest')?.time;
  if (digestStart != null) {
    g.append('line')
      .attr('x1', xScale(digestStart)).attr('x2', xScale(digestStart))
      .attr('y1', 0).attr('y2', h)
      .attr('stroke', THEME.border)
      .attr('stroke-dasharray', '4,4');
    g.append('text')
      .attr('x', xScale(digestStart) + 4).attr('y', 12)
      .attr('fill', THEME.dimText)
      .attr('font-family', THEME.font)
      .attr('font-size', '9px')
      .text('DIGEST');
  }

  // Axes
  const xAxis = d3.axisBottom(xScale).ticks(6).tickFormat(d => `${d}s`);
  const yAxisLeft = d3.axisLeft(pAgScale).ticks(5);
  const yAxisRight = d3.axisRight(sScale).ticks(4);

  const xG = g.append('g').attr('transform', `translate(0,${h})`).call(xAxis);
  const yGL = g.append('g').call(yAxisLeft);
  const yGR = g.append('g').attr('transform', `translate(${w},0)`).call(yAxisRight);

  styleAxis(xG); styleAxis(yGL); styleAxis(yGR);

  axisLabel(g, 'time', w / 2, h + 34);
  axisLabel(g, 'pAg', -h / 2, -42, -90);
  axisLabel(g, 'S', w + 42, h / 2, -90);

  // pAg line
  const pAgLine = d3.line<SimState>()
    .x(d => xScale(d.time))
    .y(d => pAgScale(Math.min(d.pAg, 20)))
    .defined(d => d.pAg < 20);

  g.append('path')
    .datum(timepoints)
    .attr('fill', 'none')
    .attr('stroke', THEME.yellow)
    .attr('stroke-width', 1.5)
    .attr('d', pAgLine);

  const sLine = d3.line<SimState>()
    .x(d => xScale(d.time))
    .y(d => sScale(d.supersaturation));

  g.append('path')
    .datum(timepoints)
    .attr('fill', 'none')
    .attr('stroke', THEME.red)
    .attr('stroke-width', 1.5)
    .attr('d', sLine);

  // Legend
  const leg = g.append('g').attr('transform', `translate(${w - 100}, 4)`);
  leg.append('line').attr('x1', 0).attr('x2', 18).attr('y1', 6).attr('y2', 6)
    .attr('stroke', THEME.yellow).attr('stroke-width', 1.5);
  leg.append('text').attr('x', 22).attr('y', 10).attr('fill', THEME.yellow)
    .attr('font-family', THEME.font).attr('font-size', '10px').text('pAg');

  leg.append('line').attr('x1', 0).attr('x2', 18).attr('y1', 22).attr('y2', 22)
    .attr('stroke', THEME.red).attr('stroke-width', 1.5);
  leg.append('text').attr('x', 22).attr('y', 26).attr('fill', THEME.red)
    .attr('font-family', THEME.font).attr('font-size', '10px').text('S');
}

// ─── Grain size distribution histogram ───────────────────────────────────────

export function drawGrainHistogram(
  el: SVGSVGElement,
  bins: Float64Array,
  edges: Float64Array,
  width: number,
  height: number
) {
  clearSvg(el);

  const svg = baseSvg(el, width, height);
  const { g, w, h } = chartArea(svg, width, height);

  const N = bins.length;
  const centers = new Float64Array(N);
  for (let i = 0; i < N; i++) centers[i] = Math.sqrt(edges[i] * edges[i + 1]);

  const xScale = d3.scaleLog()
    .domain([edges[0], edges[N]])
    .range([0, w]);

  const yMax = d3.max(bins) ?? 1;
  const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([h, 0]);

  const xAxis = d3.axisBottom(xScale)
    .tickValues([1e-9, 1e-8, 1e-7, 1e-6])
    .tickFormat(d => {
      const v = +d;
      if (v >= 1e-6) return `${v * 1e6}µm`;
      if (v >= 1e-9) return `${v * 1e9}nm`;
      return '';
    });
  const yAxis = d3.axisLeft(yScale).ticks(4);

  const xG = g.append('g').attr('transform', `translate(0,${h})`).call(xAxis);
  const yG = g.append('g').call(yAxis);
  styleAxis(xG); styleAxis(yG);

  axisLabel(g, 'grain radius', w / 2, h + 34);
  axisLabel(g, 'count', -h / 2, -42, -90);

  // Bars
  for (let i = 0; i < N; i++) {
    if (bins[i] < 1) continue;
    const x0 = xScale(edges[i]);
    const x1 = xScale(edges[i + 1]);
    const barW = Math.max(0, x1 - x0 - 1);

    g.append('rect')
      .attr('x', x0)
      .attr('y', yScale(bins[i]))
      .attr('width', barW)
      .attr('height', h - yScale(bins[i]))
      .attr('fill', THEME.yellow)
      .attr('opacity', 0.8);
  }
}

// ─── Mean grain diameter over time ───────────────────────────────────────────

export function drawGrowthCurve(
  el: SVGSVGElement,
  timepoints: SimState[],
  width: number,
  height: number
) {
  clearSvg(el);
  if (timepoints.length < 2) return;

  const svg = baseSvg(el, width, height);
  const { g, w, h } = chartArea(svg, width, height);

  const times = timepoints.map(d => d.time);
  const diameters = timepoints.map(d => d.meanRadius * 2 * 1e9); // nm

  const xScale = d3.scaleLinear().domain([0, d3.max(times)!]).range([0, w]);
  const yScale = d3.scaleLinear()
    .domain([0, (d3.max(diameters) ?? 1) * 1.1])
    .range([h, 0]);

  const xAxis = d3.axisBottom(xScale).ticks(6).tickFormat(d => `${d}s`);
  const yAxis = d3.axisLeft(yScale).ticks(5).tickFormat(d => `${d}nm`);

  const xG = g.append('g').attr('transform', `translate(0,${h})`).call(xAxis);
  const yG = g.append('g').call(yAxis);
  styleAxis(xG); styleAxis(yG);

  axisLabel(g, 'time', w / 2, h + 34);
  axisLabel(g, 'mean diameter', -h / 2, -48, -90);

  // Phase divider
  const digestStart = timepoints.find(d => d.phase === 'digest')?.time;
  if (digestStart != null) {
    g.append('line')
      .attr('x1', xScale(digestStart)).attr('x2', xScale(digestStart))
      .attr('y1', 0).attr('y2', h)
      .attr('stroke', THEME.border)
      .attr('stroke-dasharray', '4,4');
  }

  // Std deviation band
  const stdData = timepoints.filter(d => d.meanRadius > 0);
  const area = d3.area<SimState>()
    .x(d => xScale(d.time))
    .y0(d => yScale(Math.max(0, (d.meanRadius - d.stdRadius) * 2 * 1e9)))
    .y1(d => yScale((d.meanRadius + d.stdRadius) * 2 * 1e9));

  g.append('path')
    .datum(stdData)
    .attr('fill', THEME.yellow)
    .attr('opacity', 0.15)
    .attr('d', area);

  // Mean line
  const line = d3.line<SimState>()
    .x(d => xScale(d.time))
    .y(d => yScale(d.meanRadius * 2 * 1e9))
    .defined(d => d.meanRadius > 0);

  g.append('path')
    .datum(timepoints)
    .attr('fill', 'none')
    .attr('stroke', THEME.yellow)
    .attr('stroke-width', 1.5)
    .attr('d', line);
}

// ─── Nucleation rate over time ────────────────────────────────────────────────

export function drawNucleationPlot(
  el: SVGSVGElement,
  timepoints: SimState[],
  width: number,
  height: number
) {
  clearSvg(el);
  if (timepoints.length < 2) return;

  const svg = baseSvg(el, width, height);
  const { g, w, h } = chartArea(svg, width, height);

  const precipPoints = timepoints.filter(d => d.phase === 'precipitation');
  if (precipPoints.length < 2) return;

  const times = precipPoints.map(d => d.time);
  const rates = precipPoints.map(d => d.nucleationRate);

  const xScale = d3.scaleLinear().domain([0, d3.max(times)!]).range([0, w]);
  // Clip y at 95th percentile so an initial spike doesn't flatten the rest of the curve
  const sortedRates = [...rates].sort((a, b) => a - b);
  const p95 = sortedRates[Math.floor(sortedRates.length * 0.95)] ?? 1;
  const yScale = d3.scaleLinear().domain([0, p95 * 1.1]).range([h, 0]).clamp(true);

  const xAxis = d3.axisBottom(xScale).ticks(6).tickFormat(d => `${d}s`);
  const yAxis = d3.axisLeft(yScale).ticks(4);

  const xG = g.append('g').attr('transform', `translate(0,${h})`).call(xAxis);
  const yG = g.append('g').call(yAxis);
  styleAxis(xG); styleAxis(yG);

  axisLabel(g, 'time', w / 2, h + 34);
  axisLabel(g, 'nucleation rate', -h / 2, -48, -90);

  const area = d3.area<SimState>()
    .x(d => xScale(d.time))
    .y0(h)
    .y1(d => yScale(d.nucleationRate));

  g.append('path')
    .datum(precipPoints)
    .attr('fill', THEME.red)
    .attr('opacity', 0.3)
    .attr('d', area);

  const line = d3.line<SimState>()
    .x(d => xScale(d.time))
    .y(d => yScale(d.nucleationRate));

  g.append('path')
    .datum(precipPoints)
    .attr('fill', 'none')
    .attr('stroke', THEME.red)
    .attr('stroke-width', 1.5)
    .attr('d', line);
}
