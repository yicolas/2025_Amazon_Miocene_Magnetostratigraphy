#!/usr/bin/env python3
"""
Build a self-contained research-data dashboard over the geology repos in
~/Documents/GitHub.

Re-run this any time (after a `git pull`, after new measurements land) and it
regenerates research_dashboard.html from whatever is currently on disk.

    python3 build_dashboard.py

Output: research_dashboard.html  (open in any browser, no server needed)
"""

import csv
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime

# ---------------------------------------------------------------- config
ROOT = os.environ.get("GEO_ROOT", os.path.expanduser("~/Documents/GitHub"))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "research_dashboard.html")

REPOS = [
    "2024_Caborca_Tapeats",
    "2024_Laos_Australasian_Impactor",
    "2025_Amazon_Miocene_Magnetostratigraphy",
    "2026_MadreDeDios",
    "2026_pyRockmag",
]

DATA_EXT = {".csv", ".xlsx", ".dat", ".rmg", ".txt"}
SKIP_DIRS = {".git", "node_modules", ".ipynb_checkpoints", "__pycache__"}

AMZ = os.path.join(ROOT, "2025_Amazon_Miocene_Magnetostratigraphy")
ODP_CSV = os.path.join(AMZ, "ODP_202_1237_220_300_mag_data_corrected.csv")
UCHANNEL_DIR = os.path.join(AMZ, "OSU_uchannel_magnetometer_NRM_data")
CROSSPARAM = os.path.join(AMZ, "cross_parameter_comparison.csv")
PYROCKMAG = os.path.join(ROOT, "2026_pyRockmag", "pyrockmag_stats_1sample.csv")

# Demag step used for the polarity interpretation. 25 mT is the most-cleaned
# step present in the ODP file (steps are 0 / 20 / 25 mT).
POLARITY_DEMAG_MT = "25"


# ---------------------------------------------------------------- helpers
def sh(cmd, cwd=None):
    try:
        return subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, timeout=30
        ).stdout.strip()
    except Exception:
        return ""


def walk_files(repo_path):
    for dirpath, dirnames, filenames in os.walk(repo_path):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            yield os.path.join(dirpath, fn)


# ---------------------------------------------------------------- 1. repo inventory
def repo_inventory():
    repos = []
    for name in REPOS:
        path = os.path.join(ROOT, name)
        if not os.path.isdir(path):
            continue

        ext_counts = Counter()
        data_files = 0
        newest_mtime = 0.0
        total_bytes = 0

        for fp in walk_files(path):
            ext = os.path.splitext(fp)[1].lower()
            if ext not in DATA_EXT:
                continue
            ext_counts[ext.lstrip(".")] += 1
            data_files += 1
            try:
                st = os.stat(fp)
                total_bytes += st.st_size
                newest_mtime = max(newest_mtime, st.st_mtime)
            except OSError:
                pass

        last_commit = sh(
            ["git", "log", "-1", "--format=%ad|%s", "--date=short"], cwd=path
        )
        commit_date, _, commit_msg = last_commit.partition("|")
        n_commits = sh(["git", "rev-list", "--count", "HEAD"], cwd=path)

        repos.append(
            {
                "name": name,
                "data_files": data_files,
                "ext_counts": dict(ext_counts.most_common()),
                "size_mb": round(total_bytes / 1e6, 1),
                "last_data_change": (
                    datetime.fromtimestamp(newest_mtime).strftime("%Y-%m-%d")
                    if newest_mtime
                    else "—"
                ),
                "last_commit": commit_date or "—",
                "last_commit_msg": commit_msg[:60] or "—",
                "commits": int(n_commits) if n_commits.isdigit() else 0,
            }
        )
    return repos


# ---------------------------------------------------------------- 2. ODP magnetostratigraphy
def odp_series():
    """Inclination / declination / intensity vs depth at the cleaned demag step."""
    if not os.path.exists(ODP_CSV):
        return None

    by_level = defaultdict(list)
    with open(ODP_CSV, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            try:
                lvl = r["Demag_Level(mT)"].strip()
                by_level[lvl].append(
                    {
                        "depth": float(r["Depth(mbsf)"]),
                        "incl": float(r["Incl(deg)"]),
                        "decl": float(r["Decl(deg)"]),
                        "intensity": float(r["Intensity(A/m)"]),
                        "core": r["Core"].strip(),
                    }
                )
            except (ValueError, KeyError):
                continue

    if POLARITY_DEMAG_MT not in by_level:
        return None

    pts = sorted(by_level[POLARITY_DEMAG_MT], key=lambda p: p["depth"])

    normal = sum(1 for p in pts if p["incl"] > 0)
    reversed_ = sum(1 for p in pts if p["incl"] < 0)

    # Intensity decay across demag steps (median per step) — a QC view.
    decay = []
    for lvl in sorted(by_level, key=lambda x: float(x)):
        vals = sorted(p["intensity"] for p in by_level[lvl])
        if vals:
            decay.append(
                {"level": float(lvl), "median": vals[len(vals) // 2], "n": len(vals)}
            )

    return {
        "points": [
            {
                "d": round(p["depth"], 3),
                "i": round(p["incl"], 2),
                "dec": round(p["decl"], 2),
                "int": p["intensity"],
            }
            for p in pts
        ],
        "demag_levels": sorted((float(k) for k in by_level), key=float),
        "polarity_step": float(POLARITY_DEMAG_MT),
        "n_normal": normal,
        "n_reversed": reversed_,
        "depth_min": round(min(p["depth"] for p in pts), 2),
        "depth_max": round(max(p["depth"] for p in pts), 2),
        "decay": decay,
        "total_measurements": sum(len(v) for v in by_level.values()),
    }


# ---------------------------------------------------------------- 3. u-channel sections
def uchannel_summary():
    """Per-section summary of the OSU u-channel NRM runs."""
    if not os.path.isdir(UCHANNEL_DIR):
        return []

    out = []
    for fn in sorted(os.listdir(UCHANNEL_DIR)):
        if not fn.upper().endswith(".DAT"):
            continue
        fp = os.path.join(UCHANNEL_DIR, fn)
        incls, depths, ints = [], [], []
        try:
            with open(fp, newline="", encoding="utf-8", errors="replace") as fh:
                rdr = csv.DictReader(fh, delimiter="\t")
                for r in rdr:
                    try:
                        incls.append(float(r["Inclination: Unrotated"]))
                        depths.append(float(r["Depth"]))
                        ints.append(float(r["Intensity"]))
                    except (ValueError, KeyError, TypeError):
                        continue
        except OSError:
            continue

        if not incls:
            continue

        out.append(
            {
                "file": fn,
                "section": fn.split("_NRM")[0],
                "n": len(incls),
                "depth_min": round(min(depths), 2),
                "depth_max": round(max(depths), 2),
                "pct_normal": round(100 * sum(1 for i in incls if i > 0) / len(incls), 1),
                "median_intensity": sorted(ints)[len(ints) // 2],
            }
        )
    return out


# ---------------------------------------------------------------- 4. age models
def age_models():
    if not os.path.exists(CROSSPARAM):
        return []
    with open(CROSSPARAM, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- 5. rock mag stats
def rockmag():
    if not os.path.exists(PYROCKMAG):
        return []
    with open(PYROCKMAG, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- assemble
def main():
    payload = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "root": ROOT,
        "repos": repo_inventory(),
        "odp": odp_series(),
        "uchannel": uchannel_summary(),
        "age_models": age_models(),
        "rockmag": rockmag(),
    }

    html = HTML_TEMPLATE.replace("__DATA__", json.dumps(payload))
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(f"wrote {OUT}")
    print(f"  repos:        {len(payload['repos'])}")
    print(f"  data files:   {sum(r['data_files'] for r in payload['repos'])}")
    if payload["odp"]:
        o = payload["odp"]
        print(f"  ODP points:   {len(o['points'])} at {o['polarity_step']} mT "
              f"({o['n_normal']} normal / {o['n_reversed']} reversed)")
    print(f"  u-channels:   {len(payload['uchannel'])}")


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Geology Research Data Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
<style>
  :root {
    --bg: #0f1216; --panel: #171b21; --line: #262c35;
    --text: #e6e9ee; --muted: #8b95a3;
    --normal: #4b9fd5; --reversed: #d9705f; --accent: #c9a227;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; padding: 28px 32px 60px;
    background: var(--bg); color: var(--text);
    font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, sans-serif;
  }
  h1 { font-size: 21px; margin: 0 0 4px; font-weight: 600; letter-spacing: -0.01em; }
  .sub { color: var(--muted); font-size: 12.5px; margin-bottom: 26px; }
  h2 { font-size: 14px; margin: 34px 0 12px; font-weight: 600;
       text-transform: uppercase; letter-spacing: 0.07em; color: var(--muted); }
  .kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px,1fr)); gap: 12px; }
  .kpi { background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 14px 16px; }
  .kpi .v { font-size: 25px; font-weight: 650; letter-spacing: -0.02em; }
  .kpi .l { color: var(--muted); font-size: 11.5px; margin-top: 3px; }
  .panel { background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 16px 18px; }
  .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  @media (max-width: 900px) { .grid2 { grid-template-columns: 1fr; } }
  table { width: 100%; border-collapse: collapse; font-size: 12.5px; }
  th { text-align: left; color: var(--muted); font-weight: 500; padding: 7px 10px;
       border-bottom: 1px solid var(--line); font-size: 11px;
       text-transform: uppercase; letter-spacing: 0.05em; white-space: nowrap; }
  td { padding: 7px 10px; border-bottom: 1px solid rgba(38,44,53,.55); }
  tr:last-child td { border-bottom: none; }
  .num { text-align: right; font-variant-numeric: tabular-nums; }
  .repo { font-weight: 550; }
  .tag { display: inline-block; background: #222831; color: var(--muted);
         border-radius: 4px; padding: 1px 6px; margin-right: 4px; font-size: 11px;
         font-variant-numeric: tabular-nums; }
  .chartbox { position: relative; height: 460px; }
  .chartbox.short { height: 240px; }
  .legend { display: flex; gap: 16px; font-size: 12px; color: var(--muted); margin-bottom: 10px; }
  .dot { display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 5px; }
  .note { color: var(--muted); font-size: 12px; margin-top: 10px; line-height: 1.55; }
  code { background: #222831; padding: 1px 5px; border-radius: 3px; font-size: 12px; }
</style>
</head>
<body>

<h1>Geology Research Data</h1>
<div class="sub" id="sub"></div>

<div class="kpis" id="kpis"></div>

<h2>Repositories</h2>
<div class="panel">
  <table id="repotable">
    <thead><tr>
      <th>Repo</th><th class="num">Data files</th><th>Formats</th>
      <th class="num">Size</th><th class="num">Commits</th>
      <th>Last commit</th><th>Newest data file</th>
    </tr></thead>
    <tbody></tbody>
  </table>
</div>

<h2 id="odphdr">ODP Leg 202, Site 1237B — magnetostratigraphy</h2>
<div class="grid2">
  <div class="panel">
    <div class="legend">
      <span><span class="dot" style="background:var(--normal)"></span>Normal polarity (I &gt; 0)</span>
      <span><span class="dot" style="background:var(--reversed)"></span>Reversed polarity (I &lt; 0)</span>
    </div>
    <div class="chartbox"><canvas id="inclChart"></canvas></div>
    <div class="note" id="inclnote"></div>
  </div>
  <div class="panel">
    <div class="legend"><span>NRM intensity vs depth (log scale)</span></div>
    <div class="chartbox"><canvas id="intChart"></canvas></div>
    <div class="note">Intensity at the same cleaned demag step. Sharp drops often
    track lithology changes rather than field behaviour — worth cross-checking
    against the core description.</div>
  </div>
</div>

<div class="grid2" style="margin-top:16px">
  <div class="panel">
    <div class="legend"><span>Median intensity by AF demag step (QC)</span></div>
    <div class="chartbox short"><canvas id="decayChart"></canvas></div>
    <div class="note">A clean, monotonic drop from NRM through the AF steps is what
    you want to see. A flat or rising tail suggests a stubborn overprint.</div>
  </div>
  <div class="panel">
    <div class="legend"><span>Declination vs depth</span></div>
    <div class="chartbox short"><canvas id="declChart"></canvas></div>
  </div>
</div>

<h2>OSU u-channel NRM runs</h2>
<div class="panel">
  <table id="uctable">
    <thead><tr>
      <th>Section</th><th class="num">Points</th><th class="num">Depth (m)</th>
      <th class="num">% normal polarity</th><th class="num">Median intensity (A/m)</th>
    </tr></thead>
    <tbody></tbody>
  </table>
  <div class="note">% normal is a crude first pass — sign of unrotated inclination only,
  no cleaning applied. Sections near 50% are likely straddling a reversal.</div>
</div>

<h2>Age–depth models (cross-parameter comparison)</h2>
<div class="panel">
  <table id="amtable"><thead><tr></tr></thead><tbody></tbody></table>
</div>

<h2>Rock magnetic parameters (pyRockmag)</h2>
<div class="panel">
  <table id="rmtable"><thead><tr></tr></thead><tbody></tbody></table>
</div>

<script>
const D = __DATA__;

const fmt = n => n === null || n === undefined || n === '' ? '—'
  : (typeof n === 'number' ? n.toLocaleString(undefined,{maximumFractionDigits:2}) : n);
const sci = v => { const n = parseFloat(v); return isNaN(n) ? '—' : n.toExponential(2); };

/* ---------- header + KPIs ---------- */
const totalData = D.repos.reduce((s,r)=>s+r.data_files,0);
const totalSize = D.repos.reduce((s,r)=>s+r.size_mb,0);
document.getElementById('sub').textContent =
  `${D.repos.length} repositories · generated ${D.generated} · re-run build_dashboard.py to refresh`;

const kpis = [
  ['Data files', totalData.toLocaleString()],
  ['ODP measurements', D.odp ? D.odp.total_measurements.toLocaleString() : '—'],
  ['Depth interval', D.odp ? `${D.odp.depth_min}–${D.odp.depth_max}` : '—', 'mbsf'],
  ['Normal / reversed', D.odp ? `${D.odp.n_normal} / ${D.odp.n_reversed}` : '—'],
  ['U-channel sections', D.uchannel.length],
  ['Data on disk', totalSize.toFixed(0), 'MB'],
];
document.getElementById('kpis').innerHTML = kpis.map(([l,v,u])=>
  `<div class="kpi"><div class="v">${v}${u?`<span style="font-size:13px;color:var(--muted);font-weight:400"> ${u}</span>`:''}</div><div class="l">${l}</div></div>`
).join('');

/* ---------- repo table ---------- */
document.querySelector('#repotable tbody').innerHTML = D.repos.map(r => `
  <tr>
    <td class="repo">${r.name}</td>
    <td class="num">${r.data_files.toLocaleString()}</td>
    <td>${Object.entries(r.ext_counts).map(([e,c])=>`<span class="tag">${e} ${c}</span>`).join('') || '—'}</td>
    <td class="num">${r.size_mb ? r.size_mb + ' MB' : '—'}</td>
    <td class="num">${r.commits || '—'}</td>
    <td>${r.last_commit}<div style="color:var(--muted);font-size:11px">${r.last_commit_msg}</div></td>
    <td>${r.last_data_change}</td>
  </tr>`).join('');

/* ---------- ODP charts ---------- */
const gridColor = 'rgba(255,255,255,.06)';
const baseOpts = {
  responsive: true, maintainAspectRatio: false,
  plugins: { legend: { display: false } },
  scales: {
    y: { reverse: true, title: { display: true, text: 'Depth (mbsf)', color: '#8b95a3' },
         grid: { color: gridColor }, ticks: { color: '#8b95a3' } },
    x: { grid: { color: gridColor }, ticks: { color: '#8b95a3' } }
  }
};

if (D.odp) {
  const P = D.odp.points;
  document.getElementById('odphdr').textContent =
    `ODP Leg 202, Site 1237B — magnetostratigraphy (${D.odp.polarity_step} mT demag step)`;
  document.getElementById('inclnote').textContent =
    `${P.length.toLocaleString()} measurements at the ${D.odp.polarity_step} mT step. ` +
    `Available steps: ${D.odp.demag_levels.join(', ')} mT. Polarity assigned by sign of inclination only — ` +
    `no statistical cleaning, so treat this as a first look, not an interpretation.`;

  new Chart(document.getElementById('inclChart'), {
    type: 'scatter',
    data: { datasets: [{
      data: P.map(p => ({ x: p.i, y: p.d })),
      pointRadius: 1.6,
      backgroundColor: P.map(p => p.i > 0 ? 'rgba(75,159,213,.75)' : 'rgba(217,112,95,.75)')
    }]},
    options: { ...baseOpts,
      plugins: { legend:{display:false}, tooltip: { callbacks: {
        label: c => `incl ${c.parsed.x.toFixed(1)}°  @  ${c.parsed.y.toFixed(2)} mbsf` }}},
      scales: { ...baseOpts.scales,
        x: { ...baseOpts.scales.x, min: -90, max: 90,
             title: { display: true, text: 'Inclination (°)', color: '#8b95a3' }}}
    }
  });

  new Chart(document.getElementById('intChart'), {
    type: 'scatter',
    data: { datasets: [{
      data: P.filter(p => p.int > 0).map(p => ({ x: p.int, y: p.d })),
      pointRadius: 1.4, backgroundColor: 'rgba(201,162,39,.6)'
    }]},
    options: { ...baseOpts,
      scales: { ...baseOpts.scales,
        x: { type: 'logarithmic', grid: { color: gridColor }, ticks: { color: '#8b95a3' },
             title: { display: true, text: 'Intensity (A/m)', color: '#8b95a3' }}}
    }
  });

  new Chart(document.getElementById('declChart'), {
    type: 'scatter',
    data: { datasets: [{
      data: P.map(p => ({ x: p.dec, y: p.d })),
      pointRadius: 1.2, backgroundColor: 'rgba(139,149,163,.55)'
    }]},
    options: { ...baseOpts,
      scales: { ...baseOpts.scales,
        x: { ...baseOpts.scales.x, min: 0, max: 360,
             title: { display: true, text: 'Declination (°)', color: '#8b95a3' }}}
    }
  });

  new Chart(document.getElementById('decayChart'), {
    type: 'line',
    data: {
      labels: D.odp.decay.map(d => d.level + ' mT'),
      datasets: [{
        data: D.odp.decay.map(d => d.median),
        borderColor: '#c9a227', backgroundColor: 'rgba(201,162,39,.12)',
        pointRadius: 4, tension: .25, fill: true
      }]
    },
    options: {
      responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: false }, tooltip: { callbacks: {
        label: c => `median ${c.parsed.y.toExponential(2)} A/m  (n=${D.odp.decay[c.dataIndex].n})` }}},
      scales: {
        y: { type: 'logarithmic', grid: { color: gridColor }, ticks: { color: '#8b95a3' },
             title: { display: true, text: 'Median intensity (A/m)', color: '#8b95a3' }},
        x: { grid: { color: gridColor }, ticks: { color: '#8b95a3' }}
      }
    }
  });
}

/* ---------- u-channel table ---------- */
document.querySelector('#uctable tbody').innerHTML = D.uchannel.map(u => `
  <tr>
    <td class="repo">${u.section}</td>
    <td class="num">${u.n.toLocaleString()}</td>
    <td class="num">${u.depth_min} – ${u.depth_max}</td>
    <td class="num">${u.pct_normal}%</td>
    <td class="num">${u.median_intensity.toExponential(2)}</td>
  </tr>`).join('') || '<tr><td colspan="5" style="color:var(--muted)">No u-channel files found.</td></tr>';

/* ---------- generic tables ---------- */
function genericTable(sel, rows, limit) {
  const t = document.querySelector(sel);
  if (!rows.length) {
    t.querySelector('tbody').innerHTML =
      '<tr><td style="color:var(--muted)">No data found.</td></tr>';
    return;
  }
  const cols = Object.keys(rows[0]);
  t.querySelector('thead tr').innerHTML =
    cols.map(c => `<th class="${isNaN(parseFloat(rows[0][c])) ? '' : 'num'}">${c}</th>`).join('');
  t.querySelector('tbody').innerHTML = rows.slice(0, limit).map(r => `
    <tr>${cols.map(c => {
      const v = r[c];
      const n = parseFloat(v);
      const isNum = v !== '' && !isNaN(n);
      const disp = isNum && (Math.abs(n) < 1e-3 && n !== 0) ? n.toExponential(2) : (v === '' ? '—' : v);
      return `<td class="${isNum ? 'num' : ''}">${disp}</td>`;
    }).join('')}</tr>`).join('');
}
genericTable('#amtable', D.age_models, 20);
genericTable('#rmtable', D.rockmag, 20);
</script>
</body>
</html>
"""


if __name__ == "__main__":
    sys.exit(main())
