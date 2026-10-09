import * as d3 from "d3";
import "./style.css";
import { benjaminiHochberg, buildKegg, mannWhitney, pathwayAbundance, pathwayEvidence, sharedKoSheets, sharedLinks } from "./core";
import type { KeggRef, Scheme } from "./core";

// ------------------------------------------------------------------------------------------- interface
// The app draws itself into #pacfm-app, so the same build runs standalone (index.html) and inside the lab site.
// data-assets on the root says where kegg_reference.json and the demo files are (default: next to the page).
const root = document.getElementById("pacfm-app")!;
const ASSETS = (root.dataset.assets ?? "./").replace(/\/?$/, "/");
// data-api (lab site): pathway results are computed on the server, which keeps the KEGG reference; the browser only
// receives pathway-level rows and links. Without it (standalone), the app computes locally from a KEGG reference file.
const API = root.dataset.api ?? "";
const CSRF = root.dataset.csrf ?? "";
const themeParam = new URLSearchParams(location.search).get("theme");       // ?theme=light|dark overrides the system
if (themeParam === "light" || themeParam === "dark") root.dataset.theme = themeParam;
// small line icons for the ribbon (24px viewBox, stroked with currentColor)
const ICONS: Record<string, string> = {
  table: '<rect x="3" y="4" width="18" height="16" rx="2"/><path d="M3 10h18M3 15h18M9 4v16"/>',
  groups: '<circle cx="8" cy="9" r="3"/><circle cx="16" cy="9" r="3"/><path d="M3 19c.8-3 2.8-4.5 5-4.5s4.2 1.5 5 4.5M11 19c.8-3 2.8-4.5 5-4.5s4.2 1.5 5 4.5"/>',
  book: '<path d="M4 5a2 2 0 0 1 2-2h13v16H6a2 2 0 0 0-2 2z"/><path d="M4 21V5"/>',
  play: '<circle cx="12" cy="12" r="9"/><path d="M10 8.5v7l6-3.5z"/>',
  minus: '<circle cx="11" cy="11" r="7"/><path d="M8 11h6M20 20l-4-4"/>',
  plus: '<circle cx="11" cy="11" r="7"/><path d="M8 11h6M11 8v6M20 20l-4-4"/>',
  fit: '<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/>',
  image: '<rect x="3" y="4" width="18" height="16" rx="2"/><circle cx="9" cy="10" r="2"/><path d="M21 16l-5-5-8 9"/>',
  download: '<path d="M12 4v11M7 10l5 5 5-5M5 20h14"/>',
  split: '<rect x="3" y="4" width="18" height="16" rx="2"/><path d="M14 4v16M16.5 9h2M16.5 12h2M16.5 15h2"/><circle cx="8.5" cy="12" r="3"/>',
  chart: '<rect x="3" y="4" width="18" height="16" rx="2"/><circle cx="12" cy="12" r="4.5"/>',
  rows: '<rect x="3" y="4" width="18" height="16" rx="2"/><path d="M3 9h18M3 14h18M3 19h18"/>',
  expand: '<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5M9 9l-4-4M15 9l4-4M9 15l-4 4M15 15l4 4"/>',
  collapse: '<path d="M9 4v5H4M15 4v5h5M9 20v-5H4M15 20v-5h5"/>',
};
const ico = (name: string) => `<svg class="pf-ico" viewBox="0 0 24 24" aria-hidden="true">${ICONS[name]}</svg>`;

// An Office-style ribbon: tabs, each showing groups of controls with the group's name underneath.
root.innerHTML = `
  ${root.hasAttribute("data-standalone") ? `<header class="pf-head"><h1>PACFM 2 <span>pathway analysis with shared-enzyme awareness for functional metagenomics</span></h1></header>` : ""}
  <div class="pf-ribbon">
    <div class="pf-tabs" role="tablist" aria-label="PACFM 2 toolbar">
      <button type="button" role="tab" id="pf-tab-data" aria-controls="pf-panel-data" aria-selected="true">Data</button>
      <button type="button" role="tab" id="pf-tab-analysis" aria-controls="pf-panel-analysis" aria-selected="false" tabindex="-1">Analysis</button>
      <button type="button" role="tab" id="pf-tab-view" aria-controls="pf-panel-view" aria-selected="false" tabindex="-1">View</button>
      <button type="button" role="tab" id="pf-tab-export" aria-controls="pf-panel-export" aria-selected="false" tabindex="-1">Export</button>
    </div>
    <div class="pf-panel" role="tabpanel" id="pf-panel-data" aria-labelledby="pf-tab-data">
      <div class="pf-group">
        <div class="pf-group-body">
          <label class="pf-big" title="KO abundance table: TSV/CSV, HUMAnN, eggNOG-mapper or KofamScan output">${ico("table")}<span>KO table</span>
            <input type="file" id="pf-koFile" accept=".tsv,.txt,.csv,.tab,.annotations" /></label>
          <label class="pf-big" title="Sample groups (optional): two columns, sample and group">${ico("groups")}<span>Sample groups</span>
            <input type="file" id="pf-groupFile" accept=".tsv,.txt,.csv" /></label>
          <label class="pf-big" title="KEGG reference (optional; the site provides one)">${ico("book")}<span>KEGG reference</span>
            <input type="file" id="pf-keggFile" accept=".json" /></label>
        </div>
        <div class="pf-group-name">Open</div>
      </div>
      <div class="pf-group">
        <div class="pf-group-body"><button id="pf-demo" type="button" class="pf-big">${ico("play")}<span>Load demo data</span></button></div>
        <div class="pf-group-name">Example</div>
      </div>
    </div>
    <div class="pf-panel" role="tabpanel" id="pf-panel-analysis" aria-labelledby="pf-tab-analysis" hidden>
      <div class="pf-group">
        <div class="pf-group-body pf-stack">
          <label>Scheme
            <select id="pf-scheme">
              <option value="fraccov" selected>coverage-weighted fractional (recommended)</option>
              <option value="emcov">coverage-weighted EM</option>
              <option value="em">evidence-proportional (EM)</option>
              <option value="fractional">fractional (1/n)</option>
              <option value="cutoff">association cutoff (PACFM 2016)</option>
              <option value="naive">naive (sum of all KOs)</option>
            </select>
          </label>
        </div>
        <div class="pf-group-name">Counting</div>
      </div>
      <div class="pf-group" id="pf-cutoffRow">
        <div class="pf-group-body pf-stack"><label>KOs in at most <output id="pf-cutoffVal">2</output> pathway(s)
          <input type="range" id="pf-cutoff" min="1" max="20" value="2" /></label></div>
        <div class="pf-group-name">Cutoff</div>
      </div>
      <div class="pf-group">
        <div class="pf-group-body pf-stack"><label class="pf-check"><input type="checkbox" id="pf-gate" checked /> evidence gate</label>
          <span class="muted">Group results need a pathway-unique KO<br>or ≥ 50 % KO coverage.</span></div>
        <div class="pf-group-name">Group test</div>
      </div>
    </div>
    <div class="pf-panel" role="tabpanel" id="pf-panel-view" aria-labelledby="pf-tab-view" hidden>
      <div class="pf-group">
        <div class="pf-group-body pf-stack"><label>Profile <select id="pf-view"></select></label></div>
        <div class="pf-group-name">Profile</div>
      </div>
      <div class="pf-group">
        <div class="pf-group-body pf-stack">
          <label title="How shared KOs are drawn between pathways">Style
            <select id="pf-linkStyle">
              <option value="sheets" selected>sheets (pathway pairs)</option>
              <option value="lines">bundled lines (KOs)</option>
              <option value="both">sheets and lines</option>
              <option value="none">no links</option>
            </select></label>
          <label id="pf-sheetRow" title="Draw a sheet between two pathways when they share at least this many detected KOs">sheets: ≥ <output id="pf-sheetVal">10</output> shared KOs
            <input type="range" id="pf-sheetMin" min="1" max="40" value="10" /></label>
          <label id="pf-lineRow" title="Draw lines for KOs shared by at most this many pathways">lines: KOs in ≤ <output id="pf-linkVal">5</output> pathways
            <input type="range" id="pf-linkMax" min="2" max="30" value="5" /></label>
        </div>
        <div class="pf-group-name">Links</div>
      </div>
      <div class="pf-group">
        <div class="pf-group-body pf-stack"><label class="pf-check" title="Hide pathways without signal in this profile"><input type="checkbox" id="pf-hideEmpty" checked /> hide empty pathways</label></div>
        <div class="pf-group-name">Pathways</div>
      </div>
      <div class="pf-group">
        <div class="pf-group-body" role="group" aria-label="Layout">
          <button type="button" class="pf-big" data-layout="both" aria-pressed="true" title="Figure and table side by side">${ico("split")}<span>Figure + table</span></button>
          <button type="button" class="pf-big" data-layout="figure" aria-pressed="false" title="Collapse the table: figure in full width">${ico("chart")}<span>Figure only</span></button>
          <button type="button" class="pf-big" data-layout="table" aria-pressed="false" title="Collapse the figure: table in full width">${ico("rows")}<span>Table only</span></button>
          <button type="button" class="pf-big" data-fullscreen title="Show the figure in full screen (Esc to leave)">${ico("expand")}<span>Full screen</span></button>
        </div>
        <div class="pf-group-name">Layout</div>
      </div>
      <div class="pf-group">
        <div class="pf-group-body">
          <button type="button" class="pf-big" data-zoom="out" title="Zoom out (Ctrl + wheel)">${ico("minus")}<span>Zoom out</span></button>
          <output class="pf-zoomval" aria-live="polite">100 %</output>
          <button type="button" class="pf-big" data-zoom="in" title="Zoom in (Ctrl + wheel)">${ico("plus")}<span>Zoom in</span></button>
          <button type="button" class="pf-big" data-zoom="fit" title="Fit the figure to the window">${ico("fit")}<span>Fit</span></button>
        </div>
        <div class="pf-group-name">Zoom</div>
      </div>
    </div>
    <div class="pf-panel" role="tabpanel" id="pf-panel-export" aria-labelledby="pf-tab-export" hidden>
      <div class="pf-group">
        <div class="pf-group-body">
          <button id="pf-svgOut" type="button" class="pf-big">${ico("image")}<span>Figure (SVG)</span></button>
          <button id="pf-tsvOut" type="button" class="pf-big">${ico("download")}<span>Table (TSV)</span></button>
        </div>
        <div class="pf-group-name">Export</div>
      </div>
    </div>
  </div>
  <p id="pf-status" class="pf-status">${API ? "Tables are sent to the server for analysis and not stored." : "Files are read in your browser and never uploaded."}</p>
  <div class="pf-main">
    <section id="pf-figure" class="pf-figure">
      <svg id="pf-radial" class="pf-radial" role="img" aria-label="Radial pathway plot"></svg>
      <div class="pf-zoomfloat" role="group" aria-label="Zoom">
        <button type="button" data-zoom="out" aria-label="Zoom out">−</button>
        <output class="pf-zoomval">100 %</output>
        <button type="button" data-zoom="in" aria-label="Zoom in">+</button>
        <button type="button" data-zoom="fit" aria-label="Fit to window" title="Fit">${ico("fit")}</button>
        <span class="pf-sep" aria-hidden="true"></span>
        <button type="button" data-fullscreen aria-label="Full screen" title="Full screen (Esc to leave)">${ico("expand")}</button>
        <button type="button" data-layout-toggle="figure" aria-label="Figure in full width" title="Hide the table: figure in full width">${ico("chart")}</button>
      </div>
      <div id="pf-tooltip" class="pf-tooltip" hidden></div>
    </section>
    <section class="pf-tablepane"><div class="pf-panehead"><h2>Pathways</h2>
      <span><button type="button" class="pf-mini" data-layout-toggle="table" title="Hide the figure: table in full width">${ico("rows")}<span>full width</span></button>
      <button type="button" class="pf-mini" data-layout="figure" title="Hide the table">${ico("collapse")}<span>hide</span></button></span></div>
      <div id="pf-table"></div></section>
  </div>`;

// ------------------------------------------------------------------------------------------- state
let kegg: KeggRef | null = null;
let keggVersion = "";
let samples: string[] = [];
let data = new Map<string, Map<string, number>>();      // sample -> KO -> abundance
let groups = new Map<string, string>();                 // sample -> group
let lastRows: Row[] = [];

interface Row {
  pathway: string; name: string; levelA: string; levelB: string; value: number; naive: number;
  uniqueShare: number; coverage: number; gate?: boolean; weak?: boolean; pNaive?: number; qNaive?: number; p?: number; q?: number; naiveOnly?: boolean;
}

const $ = <T extends HTMLElement = HTMLElement>(id: string) => root.querySelector(`#pf-${id}`) as T;
const status = (s: string) => ($("status").textContent = s);
const val = (id: string) => ($(id) as HTMLInputElement | HTMLSelectElement).value;

// ------------------------------------------------------------------------------------------- input
function loadKeggJson(text: string) {
  const j = JSON.parse(text);
  kegg = buildKegg(j.kp, j.h);
  keggVersion = j.version ?? "";
}

const splitLine = (line: string) => (line.includes("\t") ? line.split("\t") : line.split(","));

/** KO x sample tables (plain or HUMAnN), eggNOG-mapper annotations (one sample), KofamScan detail output (one sample) */
function parseKoTable(text: string, fileName: string) {
  const lines = text.split(/\r?\n/).filter(l => l.length && !l.startsWith("##"));
  const out = new Map<string, Map<string, number>>();
  const base = fileName.replace(/\.[^.]+$/, "");
  const bump = (s: string, ko: string, v = 1) => {
    const m = out.get(s) ?? out.set(s, new Map()).get(s)!;
    m.set(ko, (m.get(ko) ?? 0) + v);
  };
  const header = splitLine(lines[0].replace(/^#\s*/, ""));
  const keggCol = header.findIndex(h => h === "KEGG_ko");
  if (keggCol >= 0) {
    for (const l of lines.slice(1)) {
      if (l.startsWith("#")) continue;
      for (const ko of (l.split("\t")[keggCol] ?? "").split(",")) if (/^ko:K\d{5}$/.test(ko)) bump(base, ko.slice(3));
    }
    return out;
  }
  if (lines.some(l => /^\*\s+\S+\s+K\d{5}\s/.test(l))) {
    const seen = new Set<string>();
    for (const l of lines) {
      const m = l.match(/^\*\s+(\S+)\s+(K\d{5})\s/);
      if (m && !seen.has(m[1] + m[2])) { seen.add(m[1] + m[2]); bump(base, m[2]); }
    }
    return out;
  }
  const cols = header.slice(1).map(c => c.replace(/_Abundance.*$|_RPKs?$/, ""));
  for (const l of lines.slice(1)) {
    const f = splitLine(l);
    if (f[0].includes("|")) continue;                       // HUMAnN stratified rows
    const ko = f[0].match(/K\d{5}/)?.[0];
    if (!ko) continue;
    cols.forEach((c, i) => { const v = Number(f[i + 1]); if (Number.isFinite(v) && v > 0) bump(c, ko, v); });
  }
  return out;
}

function parseGroups(text: string) {
  const g = new Map<string, string>();
  for (const l of text.split(/\r?\n/).slice(1)) {
    const f = splitLine(l);
    if (f.length >= 2 && f[0]) g.set(f[0].trim(), f[1].trim());
  }
  return g;
}

function setData(d: Map<string, Map<string, number>>) {
  data = d; samples = [...d.keys()];
  const v = $("view") as HTMLSelectElement;
  const gnames = [...new Set(samples.map(s => groups.get(s)).filter(Boolean))] as string[];
  v.innerHTML = [`<option value="__all">mean of all samples (${samples.length})</option>`,
    ...gnames.map(g => `<option value="__group:${g}">mean of group ${g}</option>`),
    ...samples.map(s => `<option value="${s}">${s}</option>`)].join("");
  const kos = new Set(samples.flatMap(s => [...d.get(s)!.keys()]));
  const mapped = kegg ? [...kos].filter(k => kegg!.koPathways.has(k)).length : 0;
  if (kegg || !API) status(`${samples.length} sample(s), ${kos.size.toLocaleString()} KOs (${mapped.toLocaleString()} in KEGG pathways)` +
    (gnames.length ? `; groups ${gnames.join(" vs ")}` : "") + `. ${keggVersion}`);
  render();
}

// ------------------------------------------------------------------------------------------- computation
function relative(m: Map<string, number>) {
  const t = [...m.values()].reduce((a, b) => a + b, 0) || 1;
  return new Map([...m].map(([k, v]) => [k, v / t]));
}

function viewProfile(view: string): Map<string, number> {
  const pick = view === "__all" ? samples : view.startsWith("__group:") ? samples.filter(s => groups.get(s) === view.slice(8)) : [view];
  const acc = new Map<string, number>();
  for (const s of pick) for (const [k, v] of relative(data.get(s)!)) acc.set(k, (acc.get(k) ?? 0) + v / pick.length);
  return acc;
}

const share = (m: Map<string, number>) => {
  const t = [...m.values()].reduce((a, b) => a + b, 0) || 1;
  return new Map([...m].map(([k, v]) => [k, v / t]));
};

type Link = [string, string, number];
type Sheet = [string, string, number, number, number];       // a, b, shared KOs, share of a's KOs, share of b's KOs
interface Result { rows: Row[]; links: Link[]; sheets: Sheet[]; mapped?: number; kos?: number; version?: string }

function computeLocal(): Result {
  if (!kegg || !samples.length) return { rows: [], links: [], sheets: [] };
  const k = kegg;
  const scheme = val("scheme") as Scheme;
  const cutoff = Number(val("cutoff"));
  const prof = viewProfile(val("view"));
  const cur = share(pathwayAbundance(k, prof, scheme, cutoff));
  const nai = share(pathwayAbundance(k, prof, "naive"));
  const ev = pathwayEvidence(k, prof);
  const rows: Row[] = [...nai.keys()].map(p => {
    const i = k.info.get(p)!, e = ev.get(p);
    return { pathway: p, name: i.name, levelA: i.levelA, levelB: i.levelB, value: cur.get(p) ?? 0, naive: nai.get(p)!,
             uniqueShare: e?.uniqueShare ?? NaN, coverage: e?.coverage ?? 0 };
  });
  const gnames = [...new Set(samples.map(s => groups.get(s)).filter(Boolean))] as string[];
  if (gnames.length === 2) {
    const per = (sch: Scheme) => samples.map(s => share(pathwayAbundance(k, relative(data.get(s)!), sch, cutoff)));
    const pn = per("naive"), pc = per(scheme);
    const g = samples.map(s => groups.get(s));
    const test = (P: Map<string, number>[], p: string) =>
      mannWhitney(P.filter((_, i) => g[i] === gnames[0]).map(m => m.get(p) ?? 0), P.filter((_, i) => g[i] === gnames[1]).map(m => m.get(p) ?? 0));
    rows.forEach(r => { r.pNaive = test(pn, r.pathway); r.p = test(pc, r.pathway); });
    const qn = benjaminiHochberg(rows.map(r => r.pNaive!)), qc = benjaminiHochberg(rows.map(r => r.p!));
    // evidence gate: a pathway-unique KO detected in at least half of the samples, or median KO coverage >= 0.5
    const evs = samples.map(s => pathwayEvidence(k, relative(data.get(s)!)));
    const median = (v: number[]) => { const a = v.slice().sort((x, y) => x - y), m = a.length >> 1; return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2; };
    const useGate = ($("gate") as HTMLInputElement).checked;
    rows.forEach((r, i) => {
      const e = evs.map(m => m.get(r.pathway));
      r.gate = e.filter(x => (x?.uniqueCoverage ?? 0) > 0).length >= samples.length / 2 || median(e.map(x => x?.coverage ?? 0)) >= 0.5;
      r.qNaive = qn[i]; r.q = qc[i];
      const sig = qc[i] < 0.05 && (!useGate || r.gate);
      r.weak = useGate && qc[i] < 0.05 && !r.gate;
      r.naiveOnly = qn[i] < 0.05 && !sig;
    });
  }
  const links = [...sharedLinks(k, prof, Number(val("linkMax")))].map(([key, v]) => [...key.split("|"), v] as Link);
  return { rows, links, sheets: sharedKoSheets(k, prof, Number(val("sheetMin"))) };
}

/** lab-site mode: send the KO table (gzip JSON) and the settings to the server; only pathway-level results come back */
async function computeRemote(): Promise<Result> {
  if (!samples.length) return { rows: [], links: [], sheets: [] };
  const kos = [...new Set(samples.flatMap(s => [...data.get(s)!.keys()]))];
  const body = JSON.stringify({ samples, kos, values: samples.map(s => kos.map(k => data.get(s)!.get(k) ?? 0)),
    groups: Object.fromEntries(groups), scheme: val("scheme"), cutoff: Number(val("cutoff")), view: val("view"),
    link_max: Number(val("linkMax")), sheet_min: Number(val("sheetMin")), gate: ($("gate") as HTMLInputElement).checked });
  const gz = await new Response(new Blob([body]).stream().pipeThrough(new CompressionStream("gzip"))).arrayBuffer();
  const r = await fetch(API, { method: "POST", body: gz, credentials: "same-origin",
    headers: { "Content-Type": "application/json", "Content-Encoding": "gzip", "X-CSRFToken": CSRF } });
  const j = await r.json();
  if (!r.ok) throw new Error(j.error ?? r.statusText);
  for (const row of j.rows as Row[]) if (row.uniqueShare === null) row.uniqueShare = NaN;
  return j as Result;
}

let renderSeq = 0;
async function compute(): Promise<Result | null> {
  const seq = ++renderSeq;
  if (kegg || !API) return computeLocal();
  status("Computing pathways on the server…");
  try {
    const res = await computeRemote();
    if (seq !== renderSeq) return null;               // a newer request has been sent meanwhile
    if (res.kos !== undefined) status(`${samples.length} sample(s), ${res.kos.toLocaleString()} KOs (${(res.mapped ?? 0).toLocaleString()} in KEGG pathways). ${res.version ?? ""}`);
    return res;
  } catch (e) {
    if (seq === renderSeq) status(`Server error: ${(e as Error).message}`);
    return null;
  }
}

// ------------------------------------------------------------------------------------------- radial figure
const ORDER_A = ["Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes",
                 "Organismal Systems", "Human Diseases", "Drug Development"];
const COL_A = d3.scaleOrdinal<string>().domain(ORDER_A).range(["#2a78d6", "#1baf7a", "#8e6bd6", "#d6a12a", "#898781", "#eb6834", "#52514e"]);

interface Node { name: string; children?: Node[]; row?: Row }

async function render() {
  const res = await compute();
  if (!res) return;
  const rows = res.rows; lastRows = rows;
  drawTable(rows);
  const svg = layer;                                  // the zoomable layer inside the svg (see zoom below)
  svg.selectAll("*").remove();
  if (!rows.length) return;
  const hide = ($("hideEmpty") as HTMLInputElement).checked;
  const shown = rows.filter(r => !hide || r.naive > 0).sort((a, b) =>
    (ORDER_A.indexOf(a.levelA) - ORDER_A.indexOf(b.levelA)) || a.levelB.localeCompare(b.levelB) || a.name.localeCompare(b.name));
  const W = 900, R = W / 2, r0 = 190, rBar = 255, maxBar = 120;
  // the legend gets its own column left of the circle (LEGEND_W), so it never covers bars or labels
  figure.attr("viewBox", `${-R - LEGEND_W} ${-R} ${W + LEGEND_W} ${W}`).attr("font-family", "system-ui, sans-serif");
  const tree: Node = { name: "root", children: [] };
  for (const r of shown) {
    let a = tree.children!.find(c => c.name === r.levelA); if (!a) tree.children!.push(a = { name: r.levelA, children: [] });
    let b = a.children!.find(c => c.name === r.levelB); if (!b) a.children!.push(b = { name: r.levelB, children: [] });
    b.children!.push({ name: r.pathway, row: r });
  }
  const root = d3.hierarchy<Node>(tree);
  const clustered = d3.cluster<Node>().size([2 * Math.PI, r0]).separation(() => 1)(root);
  const leaves = clustered.leaves();
  const byId = new Map(leaves.map(l => [l.data.name, l]));
  const step = (2 * Math.PI) / leaves.length;
  const arc = d3.arc<{ innerRadius: number; outerRadius: number; startAngle: number; endAngle: number }>();
  for (const [depth, inner, outer] of [[1, r0 + 34, r0 + 50], [2, r0 + 14, r0 + 30]] as const) {
    svg.append("g").selectAll("path").data(clustered.descendants().filter(d => d.depth === depth)).join("path")
      .attr("d", d => { const ls = d.leaves(); return arc({ innerRadius: inner, outerRadius: outer,
        startAngle: ls[0].x - step / 2, endAngle: ls[ls.length - 1].x + step / 2 }); })
      .attr("fill", d => COL_A((depth === 1 ? d : d.parent!).data.name)).attr("opacity", depth === 1 ? 0.9 : 0.55)
      .attr("stroke", "var(--bg, #fff)").attr("stroke-width", 1)
      .on("mousemove", (ev, d) => tip(ev, `<b>${d.data.name}</b>`)).on("mouseleave", hideTip);
  }
  const vmax = d3.max(shown, r => Math.max(r.value, r.naive)) || 1;
  const h = d3.scaleSqrt().domain([0, vmax]).range([0, maxBar]);
  const shade = d3.scaleSequential(d3.interpolateBlues).domain([-0.15, 1]);
  const bar = (r: Row, v: number) => { const l = byId.get(r.pathway)!; return arc({ innerRadius: rBar, outerRadius: rBar + h(v),
    startAngle: l.x - step * 0.42, endAngle: l.x + step * 0.42 }); };
  const g = svg.append("g");
  g.selectAll("path.naive").data(shown).join("path").attr("class", "naive").attr("d", r => bar(r, r.naive))
    .attr("fill", "none").attr("stroke", r => (r.naiveOnly ? "#eb6834" : "#898781")).attr("stroke-width", r => (r.naiveOnly ? 1.6 : 0.6));
  g.selectAll("path.cur").data(shown).join("path").attr("class", "cur").attr("d", r => bar(r, r.value))
    .attr("fill", r => (Number.isFinite(r.uniqueShare) ? shade(r.uniqueShare) : "#cccccc"))
    .on("mousemove", (ev, r) => tip(ev, `<b>${r.name}</b> (map${r.pathway})<br>${r.levelB}<br>current ${(100 * r.value).toFixed(2)} % · naive ${(100 * r.naive).toFixed(2)} %` +
      `<br>from pathway-unique KOs: ${Number.isFinite(r.uniqueShare) ? (100 * r.uniqueShare).toFixed(0) + " %" : "– (no unique KO)"}` +
      `<br>KO coverage ${(100 * r.coverage).toFixed(0)} %` + (r.q !== undefined ? `<br>q naive ${r.qNaive!.toExponential(1)} · q current ${r.q.toExponential(1)}` + (r.gate ? "" : "<br>evidence gate not passed") : "")))
    .on("mouseleave", hideTip);
  const links = res.links
    .map(([a, b, v]) => ({ a: byId.get(a), b: byId.get(b), v }))
    .filter(l => l.a && l.b).sort((x, y) => y.v - x.v).slice(0, 600);
  const style = val("linkStyle");
  if (style === "lines" || style === "both") {
    const lw = d3.scaleSqrt().domain([0, d3.max(links, l => l.v) || 1]).range([0.2, 3]);
    const line = d3.lineRadial<d3.HierarchyPointNode<Node>>().curve(d3.curveBundle.beta(0.85)).radius(d => d.y).angle(d => d.x);
    svg.insert("g", ":first-child").attr("fill", "none").selectAll("path").data(links).join("path")
      .attr("d", l => line(l.a!.path(l.b!))).attr("stroke", "#2a78d6").attr("stroke-opacity", 0.25).attr("stroke-width", l => lw(l.v));
  }
  // sheets: one ribbon per pathway pair sharing >= N detected KOs; width ~ shared KOs (within the pathway's slot),
  // colour = top-level category of the pathway with the larger share of its KOs in the overlap
  const rowById = new Map(shown.map(r => [r.pathway, r]));
  const sheets = style === "lines" || style === "none" ? [] : res.sheets.filter(([a, b]) => byId.has(a) && byId.has(b));
  const sw = d3.scaleSqrt().domain([0, d3.max(sheets, x => x[2]) || 1]).range([0.5, 3.2]);    // in pathway slots
  const ribbon = d3.ribbon<unknown, { startAngle: number; endAngle: number; radius: number }>();
  const end = (p: string, n: number) => { const x = byId.get(p)!.x, w = (step * sw(n)) / 2; return { startAngle: x - w, endAngle: x + w, radius: r0 - 3 }; };
  const sheetG = svg.insert("g", ":first-child").attr("class", "pf-sheets");
  const pct = (f: number) => `${(100 * f).toFixed(0)} %`;
  const sheetPaths = sheetG.selectAll("path").data(sheets).join("path")
    .attr("d", x => ribbon({ source: end(x[0], x[2]), target: end(x[1], x[2]) } as never) as unknown as string)
    .attr("fill", x => COL_A(rowById.get(x[0])!.levelA)).attr("fill-opacity", 0.28)
    .attr("stroke", x => COL_A(rowById.get(x[0])!.levelA)).attr("stroke-opacity", 0.45).attr("stroke-width", 0.4)
    .on("mousemove", (ev, x) => {
      d3.select(ev.currentTarget as SVGPathElement).raise().attr("fill-opacity", 0.75);
      const A = rowById.get(x[0])!, B = rowById.get(x[1])!;
      tip(ev, `<b>${A.name}</b> ↔ <b>${B.name}</b><br>${x[2]} shared KOs (detected)<br>` +
              `${pct(x[3])} of the KOs of ${A.name} · ${pct(x[4])} of ${B.name}<br>${A.levelA} · ${B.levelA}`);
    })
    .on("mouseleave", ev => { d3.select(ev.currentTarget as SVGPathElement).attr("fill-opacity", 0.28); hideTip(); });
  // hovering a pathway bar highlights its sheets
  g.selectAll<SVGPathElement, Row>("path.cur")
    .on("mouseenter.sheets", (_, r) => sheetPaths.attr("fill-opacity", x => (x[0] === r.pathway || x[1] === r.pathway ? 0.8 : 0.06))
                                                 .attr("stroke-opacity", x => (x[0] === r.pathway || x[1] === r.pathway ? 0.9 : 0.1)))
    .on("mouseleave.sheets", () => sheetPaths.attr("fill-opacity", 0.28).attr("stroke-opacity", 0.45));
  const top = new Set(shown.slice().sort((a, b) => b.value - a.value).slice(0, 24).map(r => r.pathway));
  svg.append("g").selectAll("text").data(shown.filter(r => top.has(r.pathway))).join("text")
    .attr("transform", r => { const l = byId.get(r.pathway)!; const a = (l.x * 180) / Math.PI - 90; const rr = rBar + h(Math.max(r.value, r.naive)) + 4;
      return `rotate(${a}) translate(${rr},0)` + (l.x > Math.PI ? " rotate(180)" : ""); })
    .attr("text-anchor", r => (byId.get(r.pathway)!.x > Math.PI ? "end" : "start")).attr("dy", "0.32em")
    .attr("font-size", 9).attr("fill", "var(--ink2)").text(r => (r.name.length > 34 ? r.name.slice(0, 32) + "…" : r.name));
  drawLegend(shown, shade, style, sheets.length);
}

// ------------------------------------------------------------------------------------------- legend
// Drawn in the svg itself (outside the zoom layer), so it is part of screenshots and of the SVG export.
const LEGEND_W = 330;
const SCHEME_LABEL: Record<string, string> = { naive: "naive summation", cutoff: "association cutoff", fractional: "fractional (1/n)",
  em: "evidence-proportional (EM)", fraccov: "coverage-weighted fractional", emcov: "coverage-weighted EM" };
function drawLegend(shown: Row[], shade: d3.ScaleSequential<string>, style: string, nSheets: number) {
  figure.selectAll("g.pf-legend").remove();
  const ink = "var(--ink, #0b0b0b)", ink2 = "var(--ink2, #52514e)";
  const lg = figure.append("g").attr("class", "pf-legend").attr("font-family", "system-ui, sans-serif")
    .attr("transform", `translate(${-450 - LEGEND_W + 8},-442)`).style("pointer-events", "none");
  const box = lg.append("rect").attr("x", -6).attr("y", -6).attr("rx", 4)
    .attr("fill", "var(--bg, #fff)").attr("fill-opacity", 0.88).attr("stroke", "var(--grid, #e4e3df)");
  let y = 0;
  const title = (t: string) => { lg.append("text").attr("x", 0).attr("y", y + 11).attr("font-size", 13).attr("font-weight", 600).attr("fill", ink).text(t); y += 19; };
  const item = (draw: (g: d3.Selection<SVGGElement, unknown, null, undefined>) => void, label: string) => {
    const g = lg.append("g").attr("transform", `translate(0,${y})`); draw(g);
    g.append("text").attr("x", 32).attr("y", 9.5).attr("font-size", 11.5).attr("fill", ink2).text(label); y += 17;
  };
  title("KEGG top-level category");
  const cats = ORDER_A.filter(a => shown.some(r => r.levelA === a));
  for (const c of cats) item(g => g.append("rect").attr("width", 24).attr("height", 10).attr("rx", 2).attr("fill", COL_A(c)).attr("opacity", 0.9), c);
  y += 6; title("Pathway bars");
  item(g => [0, 0.5, 1].forEach((f, i) => g.append("rect").attr("x", i * 8).attr("width", 8).attr("height", 10).attr("fill", shade(f))),
       "fill: share from pathway-specific KOs (0 → 100 %)");
  item(g => g.append("rect").attr("width", 24).attr("height", 10).attr("fill", "#cccccc"), "grey: pathway without specific KOs");
  item(g => g.append("rect").attr("x", 1).attr("y", 1).attr("width", 22).attr("height", 8).attr("fill", "none").attr("stroke", "#898781"),
       "outline: naive summation");
  if (shown.some(r => r.naiveOnly))
    item(g => g.append("rect").attr("x", 1).attr("y", 1).attr("width", 22).attr("height", 8).attr("fill", "none").attr("stroke", "#eb6834").attr("stroke-width", 1.6),
         "orange outline: significant with naive summation only");
  if (style !== "none") { y += 6; title("Centre"); }
  if (style === "sheets" || style === "both")
    item(g => g.append("path").attr("d", "M0,10 C8,10 8,0 24,0 L24,6 C12,6 10,10 0,10Z").attr("fill", COL_A(cats[0] ?? ORDER_A[0])).attr("fill-opacity", 0.45),
         `sheets: pathway pairs sharing ≥ ${val("sheetMin")} detected KOs (${nSheets})`);
  if (style === "lines" || style === "both")
    item(g => g.append("path").attr("d", "M0,9 C10,9 14,1 24,1").attr("fill", "none").attr("stroke", "#2a78d6").attr("stroke-opacity", 0.6),
         `lines: KOs shared by ≤ ${val("linkMax")} pathways`);
  y += 4;
  lg.append("text").attr("x", 0).attr("y", y + 10).attr("font-size", 11).attr("fill", ink2)
    .text(`Scheme: ${SCHEME_LABEL[val("scheme")] ?? val("scheme")}${val("scheme") === "cutoff" ? ` (KOs in ≤ ${val("cutoff")} pathways)` : ""}`);
  y += 15;
  const bb = (lg.node() as SVGGElement).getBBox();
  box.attr("width", bb.width + 12).attr("height", y + 8);
}

function tip(ev: MouseEvent, html: string) {
  const t = $("tooltip"); t.innerHTML = html; t.hidden = false;
  const b = $("figure").getBoundingClientRect();
  t.style.left = `${ev.clientX - b.left + 12}px`; t.style.top = `${ev.clientY - b.top + 12}px`;
}
const hideTip = () => { $("tooltip").hidden = true; };

// ------------------------------------------------------------------------------------------- table
let sortKey: keyof Row = "value", sortDesc = true;
function drawTable(rows: Row[]) {
  const hasG = rows.some(r => r.q !== undefined);
  type Col = [keyof Row, string, (r: Row) => string];
  const cols: Col[] = [
    ["name", "pathway", r => r.name + (r.weak ? '<span class="badge weak">weak evidence</span>' : r.naiveOnly ? '<span class="badge">naive only</span>' : "")],
    ["value", "current %", r => (100 * r.value).toFixed(2)], ["naive", "naive %", r => (100 * r.naive).toFixed(2)],
    ["uniqueShare", "unique KO %", r => (Number.isFinite(r.uniqueShare) ? (100 * r.uniqueShare).toFixed(0) : "–")],
    ["coverage", "KO cov %", r => (100 * r.coverage).toFixed(0)]];
  if (hasG) cols.push(["qNaive", "q naive", r => r.qNaive!.toExponential(1)], ["q", "q current", r => r.q!.toExponential(1)]);
  const sorted = rows.slice().sort((a, b) => {
    const x = a[sortKey] as number | string, y = b[sortKey] as number | string;
    return (x < y ? -1 : x > y ? 1 : 0) * (sortDesc ? -1 : 1);
  });
  $("table").innerHTML = `<table><thead><tr>${cols.map(([key, t]) => `<th data-k="${String(key)}">${t}</th>`).join("")}</tr></thead><tbody>` +
    sorted.slice(0, 400).map(r => `<tr class="${r.weak ? "weak" : r.naiveOnly ? "flag" : ""}">${cols.map(([, , f]) => `<td>${f(r)}</td>`).join("")}</tr>`).join("") +
    "</tbody></table>";
  $("table").querySelectorAll("th").forEach(th => th.addEventListener("click", () => {
    const key = (th as HTMLElement).dataset.k as keyof Row; sortDesc = sortKey === key ? !sortDesc : true; sortKey = key; drawTable(lastRows); }));
}

// ------------------------------------------------------------------------------------------- zoom
// d3.zoom on the svg moves one layer that holds the whole figure: buttons (ribbon and the figure's
// corner), Ctrl/⌘ + wheel, and dragging once zoomed in. A plain wheel still scrolls the page.
const figure = d3.select<SVGSVGElement, unknown>($("radial") as unknown as SVGSVGElement);
const layer = figure.append("g").attr("class", "pf-layer");
const ZOOM_MIN = 0.5, ZOOM_MAX = 8;
const zoom = d3.zoom<SVGSVGElement, unknown>()
  .scaleExtent([ZOOM_MIN, ZOOM_MAX])
  .filter(ev => (ev.type === "wheel" ? ev.ctrlKey || ev.metaKey : !ev.button))
  // ~15 % per wheel step (d3's default multiplies Ctrl+wheel by 10, made for trackpad pinch: far too fast for a mouse)
  .wheelDelta(ev => -ev.deltaY * (ev.deltaMode === 1 ? 0.05 : ev.deltaMode ? 1 : 0.002))
  .on("zoom", ev => {
    layer.attr("transform", ev.transform.toString());
    root.querySelectorAll<HTMLOutputElement>(".pf-zoomval").forEach(o => (o.textContent = `${Math.round(ev.transform.k * 100)} %`));
    root.querySelectorAll<HTMLButtonElement>("[data-zoom=in]").forEach(b => (b.disabled = ev.transform.k >= ZOOM_MAX - 1e-6));
    root.querySelectorAll<HTMLButtonElement>("[data-zoom=out]").forEach(b => (b.disabled = ev.transform.k <= ZOOM_MIN + 1e-6));
    $("figure").classList.toggle("is-zoomed", ev.transform.k > 1.001);
  });
figure.call(zoom).on("dblclick.zoom", null);
root.querySelectorAll<HTMLButtonElement>("[data-zoom]").forEach(b => b.addEventListener("click", () => {
  const t = figure.transition().duration(200);
  if (b.dataset.zoom === "in") zoom.scaleBy(t, 1.25);
  else if (b.dataset.zoom === "out") zoom.scaleBy(t, 0.8);
  else zoom.transform(t, d3.zoomIdentity);
}));

// ------------------------------------------------------------------------------------------- ribbon tabs
const tabs = [...root.querySelectorAll<HTMLButtonElement>(".pf-tabs [role=tab]")];
function selectTab(tab: HTMLButtonElement, focus = false) {
  for (const t of tabs) {
    const on = t === tab;
    t.setAttribute("aria-selected", String(on));
    t.tabIndex = on ? 0 : -1;
    root.querySelector<HTMLElement>(`#${t.getAttribute("aria-controls")}`)!.hidden = !on;
  }
  if (focus) tab.focus();
}
tabs.forEach((t, i) => {
  t.addEventListener("click", () => selectTab(t));
  t.addEventListener("keydown", ev => {
    const step = ev.key === "ArrowRight" ? 1 : ev.key === "ArrowLeft" ? -1 : 0;
    if (step) { ev.preventDefault(); selectTab(tabs[(i + step + tabs.length) % tabs.length], true); }
  });
});

// ------------------------------------------------------------------------------------------- export
function download(name: string, text: string, type: string) {
  const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click();
}
$("svgOut").addEventListener("click", () => {
  const copy = ($("radial") as unknown as SVGSVGElement).cloneNode(true) as SVGSVGElement;   // exported unzoomed
  copy.querySelector(".pf-layer")?.removeAttribute("transform");
  download("pacfm2.svg", new XMLSerializer().serializeToString(copy), "image/svg+xml");
});
$("tsvOut").addEventListener("click", () => {
  const head = ["pathway", "name", "level_a", "level_b", "current", "naive", "unique_share", "coverage", "q_naive", "q_current", "evidence_gate"];
  download("pacfm2_pathways.tsv", [head.join("\t"), ...lastRows.map(r => [r.pathway, r.name, r.levelA, r.levelB, r.value, r.naive,
    r.uniqueShare, r.coverage, r.qNaive ?? "", r.q ?? "", r.gate ?? ""].join("\t"))].join("\n"), "text/tab-separated-values");
});

// ------------------------------------------------------------------------------------------- layout and full screen
// data-layout on .pf-main: both (default) | figure (table collapsed) | table (figure collapsed)
const main = root.querySelector<HTMLElement>(".pf-main")!;
function setLayout(layout: string) {
  main.dataset.layout = layout;
  root.querySelectorAll<HTMLButtonElement>("button[data-layout]").forEach(b => {
    if (b.closest(".pf-group-body")) b.setAttribute("aria-pressed", String(b.dataset.layout === layout)); });
  // the quick toggles switch back to both panes when pressed again
  root.querySelectorAll<HTMLButtonElement>("button[data-layout-toggle]").forEach(b => {
    const on = b.dataset.layoutToggle === layout;
    b.setAttribute("aria-pressed", String(on));
    const label = b.querySelector("span"); if (label) label.textContent = on ? "both panes" : "full width"; });
  try { localStorage.setItem("pacfm.layout", layout); } catch { /* storage unavailable */ }
  if (layout !== "table") zoom.transform(figure, d3.zoomIdentity);     // re-fit the figure to its new width
}
root.querySelectorAll<HTMLButtonElement>("button[data-layout]").forEach(b => b.addEventListener("click", () => setLayout(b.dataset.layout!)));
root.querySelectorAll<HTMLButtonElement>("button[data-layout-toggle]").forEach(b => b.addEventListener("click", () =>
  setLayout(main.dataset.layout === b.dataset.layoutToggle ? "both" : b.dataset.layoutToggle!)));
{ // ?layout=figure|table|both (shareable) wins over the layout remembered in this browser
  let want = new URLSearchParams(location.search).get("layout");
  if (!want) try { want = localStorage.getItem("pacfm.layout"); } catch { /* storage unavailable */ }
  if (want === "figure" || want === "table") setLayout(want);
}

// full screen shows the figure section (plot, legend, zoom buttons, tooltips); Esc or the button leaves it
const figSection = $("figure");
function toggleFullscreen() {
  if (document.fullscreenElement) { void document.exitFullscreen(); return; }
  if (main.dataset.layout === "table") setLayout("both");               // the figure must be visible first
  void figSection.requestFullscreen?.();
}
root.querySelectorAll<HTMLButtonElement>("button[data-fullscreen]").forEach(b => {
  b.addEventListener("click", toggleFullscreen);
  if (!document.fullscreenEnabled) b.disabled = true;
});
document.addEventListener("fullscreenchange", () => {
  const on = document.fullscreenElement === figSection;
  root.querySelectorAll<HTMLButtonElement>("button[data-fullscreen]").forEach(b => {
    b.setAttribute("aria-pressed", String(on)); b.innerHTML = b.classList.contains("pf-big")
      ? `${ico(on ? "collapse" : "expand")}<span>${on ? "Leave full screen" : "Full screen"}</span>` : ico(on ? "collapse" : "expand"); });
  zoom.transform(figure, d3.zoomIdentity);
});

// ------------------------------------------------------------------------------------------- wiring
const onFile = (id: string, fn: (text: string, name: string) => void) =>
  $(id).addEventListener("change", async e => { const f = (e.target as HTMLInputElement).files?.[0]; if (f) fn(await f.text(), f.name); });
onFile("koFile", (t, n) => setData(parseKoTable(t, n)));
onFile("groupFile", t => { groups = parseGroups(t); setData(data); });
onFile("keggFile", t => { loadKeggJson(t); if (samples.length) setData(data); else status(`KEGG reference loaded. ${keggVersion}`); });
async function loadDemo() {
  groups = parseGroups(await (await fetch(ASSETS + "demo_groups.tsv")).text());
  setData(parseKoTable(await (await fetch(ASSETS + "demo_ko_table.tsv")).text(), "demo"));
}
$("demo").addEventListener("click", loadDemo);
for (const id of ["scheme", "view", "hideEmpty", "gate"]) $(id).addEventListener("change", render);
let debounce = 0;
const later = () => { clearTimeout(debounce); debounce = window.setTimeout(render, API && !kegg ? 250 : 0); };
$("cutoff").addEventListener("input", () => { $("cutoffVal").textContent = val("cutoff"); later(); });
$("linkMax").addEventListener("input", () => { $("linkVal").textContent = val("linkMax"); later(); });
$("sheetMin").addEventListener("input", () => { $("sheetVal").textContent = val("sheetMin"); later(); });
const showLinkRows = () => { const st = val("linkStyle");
  $("sheetRow").style.display = st === "lines" || st === "none" ? "none" : ""; $("lineRow").style.display = st === "sheets" || st === "none" ? "none" : ""; };
$("linkStyle").addEventListener("change", () => { showLinkRows(); render(); });
showLinkRows();
$("scheme").addEventListener("change", () => { $("cutoffRow").style.display = val("scheme") === "cutoff" ? "" : "none"; });
$("cutoffRow").style.display = "none";

const wantDemo = new URLSearchParams(location.search).has("demo");
if (API) {
  status("Pathways are computed on the server; tables are sent for analysis and not stored.");
  if (wantDemo) loadDemo();
} else {
  fetch(ASSETS + "kegg_reference.json").then(r => (r.ok ? r.text() : Promise.reject(new Error("missing"))))
    .then(t => { loadKeggJson(t); status(`KEGG reference loaded (${kegg!.pathwayKos.size} pathways). ${keggVersion}`);
                 if (wantDemo) return loadDemo(); })
    .catch(() => status("No KEGG reference found: load one with 'KEGG reference' (build it with build_kegg_reference.py)."));
}
