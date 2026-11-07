#!/usr/bin/env python3
"""
Peru Amazon Miocene magnetostratigraphy — polarity synthesis
============================================================
Standalone, reproducible extraction of the "Corrected polarity" pipeline that
previously lived in Uchannel_notebook.ipynb (bottom sections). Nothing is
re-fit: this consumes the existing despiked + 300-1000 mT PCA ChRM directions
in OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed/dir/pca_corrected/.

Pipeline
--------
1. Load the despiked/300-1000 PCA directions per group.
2. Clean: drop MAD > MAD_MAX, trim EDGE_TRIM_CM off each group end.
3. Foam: keep only positions inside viable_positions (drop foam spacers),
   then compress foam gaps so composite depth = real sediment.
4. Southern-hemisphere polarity sign: at ~4 S a GAD normal field gives
   NEGATIVE inclination, so POLARITY_SIGN = -1  (normal = inc < 0).
5. Running-median de-noise the inclination, merge into robust polarity zones
   >= MIN_ZONE_CM2, stack groups into a composite section per site.
6. Emit tidy CSVs (per-sample directions + polarity zones) and figures.

Outputs (written to ./out):
  samples_corrected.csv        per-sample directions with polarity label
  polarity_zones.csv           merged polarity zones per site
  polarity_<SITE>.png          inclination (raw+smoothed) + MAD + polarity col
  TAM_vs_SC_polarity.png       side-by-side corrected polarity columns
  polarity_vs_GPTS.png         observed polarity vs mid-Miocene GPTS
"""
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- parameters
PCA_DIR = os.environ.get(
    "PCA_DIR",
    "/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy/"
    "OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed/dir/pca_corrected",
)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

MAD_MAX      = 12.0   # drop fits worse than this (deg)
EDGE_TRIM_CM = 5.0    # trim off each end of every group (u-channel edge effect)
POLARITY_SIGN = -1    # southern hemisphere (~4 S): NORMAL = negative inclination
SMOOTH_CM    = 7      # running-median window on inclination (cm)
MIN_ZONE_CM2 = 12     # minimum robust polarity-zone thickness (cm)

# ---------------------------------------------------------------- core geometry
groups = {
    "Group 1": ["TAM 1", "TAM 2", "TAM 3B"],
    "Group 2": ["TAM 3A", "TAM 4", "TAM 5E", "TAM 5D", "TAM 5C"],
    "Group 3": ["TAM 5B", "TAM 5A", "TAM 6", "TAM 7"],
    "Group 4": ["SC 1C", "SC 1B", "SC 1A", "SC 2", "SC 3C"],
    "Group 5": ["SC 3B", "SC 3A", "SC 4D", "SC 4C", "SC 4B", "SC 4A"],
    "Group 6": ["SC 5", "SC 6", "SC 7C", "SC 7B", "SC 7A"],
    "Group 7": ["SC 8"],
}
SC_cores = ["SC 1C","SC 1B","SC 1A","SC 2","SC 3C","SC 3B","SC 3A","SC 4D","SC 4C","SC 4B",
            "SC 4A","SC 5","SC 6","SC 7C","SC 7B","SC 7A","SC 8"]
SC_len   = [18.5,13,8,50,34,27.5,20,24,16,26.5,22.5,23,50,14,13.5,13,50]
SC_above = [1,2.5,4,7.5,1.2,1.5,3.6,4.5,.25,1.5,.5,2.5,3,2,1,1.5,15]
SC_below = [1,5,4,15.5,3.5,3,2.5,7,.5,2,.25,3.5,2,1.5,1.5,2,8]
TAM_cores= ["TAM 1","TAM 2","TAM 3B","TAM 3A","TAM 4","TAM 5E","TAM 5D","TAM 5C","TAM 5B",
            "TAM 5A","TAM 6","TAM 7"]
TAM_len  = [75,50.5,13.6,36,74,12,11,10,8.5,7.8,50,49]
TAM_above= [9,12.3,4.5,0,10,2.5,3,1.5,2.5,2,2.1,1]
TAM_below= [9.2,7.5,0,0,7,2.1,2,1.5,2.3,1.5,2.5,.5]

geom = {}
for c,l,a,b in list(zip(SC_cores,SC_len,SC_above,SC_below)) + \
               list(zip(TAM_cores,TAM_len,TAM_above,TAM_below)):
    geom[c] = dict(length=l, above=a, below=b)

def build_viable():
    """Viable (non-foam) position ranges per group, in with-foam u-channel cm."""
    vp = {}
    for gname, cores in groups.items():
        x, ranges = 0.0, []
        for c in cores:
            g = geom[c]
            vs, ve = x + g["above"], x + g["length"] - g["below"]
            if vs < ve:
                ranges.append([vs, ve])
            x += g["length"]
        vp[gname] = ranges
    return vp
viable_positions = build_viable()

GROUP_NAME = {"TAM_1_2_3B":"Group 1","TAM_3A_4_5E_5D_5C":"Group 2","TAM_5B_5A_6_7":"Group 3",
              "SC_1C_1B_1A_2_3C":"Group 4","SC_3B_3A_4D_4C_4B_4A":"Group 5",
              "SC_5_6_7C_7B_7A":"Group 6","SC_8":"Group 7"}
SITE_GROUPS = {
    "TAM": ["TAM_1_2_3B", "TAM_3A_4_5E_5D_5C", "TAM_5B_5A_6_7"],
    "SC":  ["SC_1C_1B_1A_2_3C", "SC_3B_3A_4D_4C_4B_4A", "SC_5_6_7C_7B_7A", "SC_8"],
}

# ---------------------------------------------------------------- load + clean
def _load_pca(path):
    rows = [r.rstrip("\n").split("\t") for r in open(path, errors="replace")]
    h = rows[0]
    def ci(*names):
        for n in names:
            if n in h: return h.index(n)
    pi, di, ii, mi = ci("Position(cm)","position"), ci("Decl"), ci("Incl"), ci("MAD")
    out = []
    for r in rows[1:]:
        try: out.append((float(r[pi]), float(r[di]), float(r[ii]), float(r[mi])))
        except: pass
    return pd.DataFrame(out, columns=["pos","dec","inc","mad"])

def group_key(fname):
    base = os.path.basename(fname).split(" - ")[1].split(" (")[0]
    return base.split("_NRM")[0]

fits = {}
for f in sorted(glob.glob(os.path.join(PCA_DIR, "PCA*"))):
    df = _load_pca(f)
    if df.empty: continue
    lo, hi = df["pos"].min()+EDGE_TRIM_CM, df["pos"].max()-EDGE_TRIM_CM
    df = df[(df["mad"] < MAD_MAX) & (df["pos"] >= lo) & (df["pos"] <= hi)].reset_index(drop=True)
    fits[group_key(f)] = df

# ---------------------------------------------------------------- foam + smooth
def _runmed(vals, pos, win):
    out = np.empty(len(vals))
    for i in range(len(vals)):
        m = (pos >= pos[i]-win/2) & (pos <= pos[i]+win/2)
        out[i] = np.median(vals[m])
    return out

def is_normal(inc):
    return (POLARITY_SIGN*inc) > 0

def _foam_keep(gkey, df):
    ranges = viable_positions.get(GROUP_NAME.get(gkey))
    if not ranges: return df.reset_index(drop=True)
    keep = df["pos"].apply(lambda p: any(a <= p <= b for a,b in ranges))
    return df[keep].reset_index(drop=True)

fits2 = {}
for g, df in fits.items():
    if df.empty: continue
    d = _foam_keep(g, df).copy()
    if len(d) < 3: continue
    d["inc_s"] = _runmed(d["inc"].values, d["pos"].values, SMOOTH_CM)
    step = np.diff(d["pos"].values, prepend=d["pos"].values[0])
    step[step > 3] = 1.0                 # collapse foam gaps (>3 cm) to nominal spacing
    d["sed"] = np.cumsum(step)
    d["polarity"] = np.where(d["inc_s"].apply(is_normal), "N", "R")
    fits2[g] = d

def zones_corr(df, min_cm=MIN_ZONE_CM2, depthcol="pos", col="inc_s"):
    s = np.array([1 if is_normal(v) else -1 for v in df[col]]); pos = df[depthcol].values
    b = [0] + [k for k in range(1, len(s)) if s[k] != s[k-1]] + [len(s)]
    return [(pos[b[j]], pos[b[j+1]-1], "N" if s[b[j]] > 0 else "R")
            for j in range(len(b)-1) if pos[b[j+1]-1] - pos[b[j]] >= min_cm]

# ---------------------------------------------------------------- composites
def composite2(site):
    parts, off = [], 0.0
    for g in SITE_GROUPS[site]:
        if g not in fits2: continue
        d = fits2[g].copy()
        d["cdepth"] = (d["sed"] - d["sed"].min()) + off
        d["site"] = site
        parts.append(d)
        off += (d["sed"].max() - d["sed"].min()) + 1
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
comp2 = {s: composite2(s) for s in SITE_GROUPS}

# ---------------------------------------------------------------- GPTS + age
GPTS_NORMAL = [
    (12.049, 12.174, "C5An.1n"), (12.272, 12.474, "C5An.2n"),
    (12.735, 12.770, "C5Ar.1n"), (12.829, 12.887, "C5Ar.2n"),
    (13.032, 13.183, "C5AAn"),   (13.363, 13.608, "C5ABn"),
    (13.739, 14.070, "C5ACn"),   (14.163, 14.609, "C5ADn"),
]
AGE_MODEL = {"TAM": (12.935, 778.0), "SC": (12.945, 786.0)}  # (age at top Ma, sed cm/Ma)
def d2age(site, cd): top, sed = AGE_MODEL[site]; return top + (cd - comp2[site]["cdepth"].min())/sed

# ---------------------------------------------------------------- console report
print("="*70)
print(f"PCA source : {PCA_DIR}")
print(f"filters    : MAD<{MAD_MAX:g}, edge-trim {EDGE_TRIM_CM:g}cm, foam-removed, "
      f"smooth {SMOOTH_CM}cm, zones>= {MIN_ZONE_CM2}cm, sign {POLARITY_SIGN:+d}")
print("="*70)
print("\nPer-group corrected polarity (N = normal = inc<0 at ~4 S):")
for g, df in fits2.items():
    z = zones_corr(df)
    seq = "".join(t[2] for t in z) or "(none)"
    print(f"  {g:22} {seq:8} ({len(z)} zones)  n={len(df):3}  "
          f"medMAD={df['mad'].median():4.1f}  medIncl={df['inc'].median():+6.1f}")
print("\nPer-site composite corrected polarity:")
zone_rows = []
for s, df in comp2.items():
    if df.empty: continue
    z = zones_corr(df, depthcol="cdepth")
    print(f"  {s:4} {''.join(t[2] for t in z):10}  "
          f"{len(df)} pts, {df['cdepth'].max():.0f} cm sediment, {len(z)} zones")
    for a, b, p in z:
        zone_rows.append(dict(site=s, top_cm=round(a,1), bottom_cm=round(b,1),
                              thickness_cm=round(b-a,1), polarity=p,
                              top_age_Ma=round(d2age(s,a),3), bottom_age_Ma=round(d2age(s,b),3)))

# ---------------------------------------------------------------- CSV outputs
samp = pd.concat([d.assign(group=g) for g, d in fits2.items()], ignore_index=True)
samp = samp[["group","pos","sed","inc","dec","mad","inc_s","polarity"]]
samp.columns = ["group","uchannel_pos_cm","sed_depth_cm","inc","dec","mad","inc_smooth","polarity"]
samp.to_csv(os.path.join(OUT,"samples_corrected.csv"), index=False)
pd.DataFrame(zone_rows).to_csv(os.path.join(OUT,"polarity_zones.csv"), index=False)
print(f"\nwrote samples_corrected.csv ({len(samp)} rows) and "
      f"polarity_zones.csv ({len(zone_rows)} zones)")

# ---------------------------------------------------------------- figures
def fig_site(site):
    df = comp2[site]
    z = zones_corr(df, depthcol="cdepth")
    fig, (a1, aM, a2) = plt.subplots(
        1, 3, figsize=(7.5, 9), gridspec_kw={"width_ratios":[5,2,1]}, sharey=True)
    a1.plot(df["inc"], df["cdepth"], lw=0.6, color="0.72", label="raw ChRM fit")
    a1.plot(df["inc_s"], df["cdepth"], lw=1.5, color="k", label=f"{SMOOTH_CM} cm median")
    a1.axvline(0, color="crimson", ls="--", lw=0.9)
    a1.set_xlim(-90, 90); a1.invert_yaxis()
    a1.set_xlabel("Inclination (deg)"); a1.set_ylabel("Composite sediment depth (cm)")
    a1.legend(fontsize=8, loc="lower right")
    a1.set_title(f"Site {site}   (normal = negative inclination, ~4 S)")
    aM.plot(df["mad"], df["cdepth"], lw=0.8, color="#2b6cb0")
    aM.axvline(MAD_MAX, color="0.5", ls=":", lw=0.9)
    aM.set_xlabel("MAD (deg)"); aM.set_xlim(0, MAD_MAX+2); aM.set_title("fit quality")
    for a, b, p in z:
        a2.axhspan(a, b, color=("k" if p == "N" else "w"), ec="k")
        a2.text(0.5, (a+b)/2, p, ha="center", va="center",
                color=("w" if p == "N" else "k"), fontsize=9, fontweight="bold")
    a2.set_xlim(0, 1); a2.set_xticks([]); a2.set_title("polarity")
    plt.tight_layout()
    p = os.path.join(OUT, f"polarity_{site}.png")
    fig.savefig(p, dpi=200, bbox_inches="tight"); plt.close(fig)
    return p

def fig_compare():
    sites = [s for s in comp2 if not comp2[s].empty]
    fig, axes = plt.subplots(1, len(sites), figsize=(2.4*len(sites), 9), sharey=False)
    if len(sites) == 1: axes = [axes]
    for ax, s in zip(axes, sites):
        df = comp2[s]
        for a, b, p in zones_corr(df, depthcol="cdepth"):
            ax.axhspan(a, b, color=("k" if p == "N" else "w"), ec="k")
            ax.text(0.5,(a+b)/2,p,ha="center",va="center",
                    color=("w" if p=="N" else "k"),fontweight="bold")
        ax.set_title(f"Site {s}"); ax.set_xlim(0,1); ax.set_xticks([]); ax.invert_yaxis()
        ax.set_ylabel("Composite sediment depth (cm)")
    fig.suptitle("TAM vs SC corrected polarity (black = Normal)", y=1.005)
    plt.tight_layout()
    p = os.path.join(OUT, "TAM_vs_SC_polarity.png")
    fig.savefig(p, dpi=200, bbox_inches="tight"); plt.close(fig)
    return p

def fig_gpts():
    sites = [s for s in comp2 if not comp2[s].empty]
    fig, axes = plt.subplots(1, len(sites)+1, figsize=(2.4*(len(sites)+1), 9), sharey=True)
    # full window: youngest (top of any record) to oldest (bottom of any record)
    young = min(d2age(s, comp2[s]["cdepth"].min()) for s in sites)
    old   = max(d2age(s, comp2[s]["cdepth"].max()) for s in sites)
    amin, amax = young, old
    axg = axes[0]
    axg.set_title("GPTS\n(GTS2020)"); axg.set_ylabel("Age (Ma)")
    axg.set_ylim(amax, amin); axg.set_xlim(0, 1); axg.set_xticks([])
    axg.axhspan(amin, amax, color="w", ec="k")
    for a, b, nm in GPTS_NORMAL:
        axg.axhspan(a, b, color="k")
        if amin <= (a+b)/2 <= amax:
            axg.text(1.08, (a+b)/2, nm, va="center", fontsize=7)
    for ax, s in zip(axes[1:], sites):
        df = comp2[s]
        for a, b, p in zones_corr(df, depthcol="cdepth"):
            aa, bb = d2age(s, a), d2age(s, b)
            ax.axhspan(min(aa,bb), max(aa,bb), color=("k" if p == "N" else "w"), ec="k")
        ax.set_title(f"Site {s}\ntop {AGE_MODEL[s][0]} Ma\n{AGE_MODEL[s][1]:.0f} cm/Ma")
        ax.set_xlim(0,1); ax.set_xticks([])
    fig.suptitle("Observed polarity vs mid-Miocene GPTS (black = Normal) — "
                 "AGE_MODEL is a provisional tie", y=1.01)
    plt.tight_layout()
    p = os.path.join(OUT, "polarity_vs_GPTS.png")
    fig.savefig(p, dpi=200, bbox_inches="tight"); plt.close(fig)
    return p

paths = [fig_site(s) for s in comp2 if not comp2[s].empty]
paths += [fig_compare(), fig_gpts()]
print("wrote figures:")
for p in paths: print("  ", os.path.basename(p))
print("\nDONE")
