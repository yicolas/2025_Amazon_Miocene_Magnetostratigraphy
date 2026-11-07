#!/usr/bin/env python3
"""Digitise the Mora-Rojas et al. (2023) La Venta composite standard section (CS) from the vector
supplement s6 ("Extended lithological resolution of the composite standard section", 1:1000),
and the reinterpreted polarity column from supplement s15.

Outputs (next to this script):
  laventa_cs_beds.csv      base_m, top_m, litho, grain, width_x, codes, + feature flags
  laventa_cs_polarity.csv  base_m, top_m, polarity, chron

The CS metres are the authors' metres (fossil localities in s10, U-Pb samples in Table 4 and the
Monkey Beds at StL 9, ~557-567 m, all use them), so everything in Fig. 2 can sit on one scale.

usage: python build_laventa_cs.py "<mora-rojas appendix folder>"
needs: pymupdf, numpy, pandas, pillow
"""
import sys, os, re
import numpy as np, pandas as pd
import pymupdf as fitz

APP = sys.argv[1] if len(sys.argv) > 1 else "."
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- s6: lithology column ----------------
page = fitz.open(os.path.join(APP, "g2023v45a6_s6.pdf"))[0]
# metre axis: label centres give 2.8347 pt per metre, 0 m at y = 3490.5 pt
Y0, PT_PER_M = 3490.5, 2.8347
def pt2m(y): return (Y0 - y) / PT_PER_M
X_LEFT, X_RIGHT, Z = 376.0, 522.0, 4
clip = fitz.Rect(X_LEFT, 180, X_RIGHT, 3700)
pix = page.get_pixmap(matrix=fitz.Matrix(Z, Z), clip=clip, alpha=False)
A = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, 3).astype(int)

FILLS = {  # vector fill colours of s6 (RGB 0-255) -> lithology
    "sandstone": (255, 255, 0), "mudstone": (204, 201, 201), "claystone": (255, 227, 224),
    "siltstone": (230, 230, 230), "conglomerate": (255, 214, 0), "mudstone_red": (240, 128, 128)}
FC = np.array(list(FILLS.values())); FN = list(FILLS)
# grain-size ticks (right edge of each class, page x) read from the s6 header
GS = [(406, "clay"), (426, "silt"), (436, "very fine sand"), (446, "fine sand"), (456, "medium sand"),
      (466, "coarse sand"), (476, "very coarse sand"), (486, "granule"), (496, "pebble"), (516, "cobble")]

rows = []
for r in range(A.shape[0]):
    line = A[r]
    d = np.abs(line[:, None, :] - FC[None, :, :]).sum(-1)          # distance to each fill colour
    best = d.argmin(1); ok = d.min(1) < 30
    lo = int((381 - X_LEFT) * Z); seg = best[lo:lo + 12][ok[lo:lo + 12]]
    if not len(seg): rows.append((r, None, np.nan)); continue
    k = np.bincount(seg).argmax()
    # bed width = rightmost pixel of the run that starts at the column's left edge: walk right through
    # fill colour and black texture/outline, stop at the first >=3 px of white (outside the bed)
    white = (line >= 245).all(1); j = lo; run = 0; last = lo
    while j < len(line) and run < 3:
        if white[j]: run += 1
        else: run = 0; last = j
        j += 1
    rows.append((r, FN[k], X_LEFT + last / Z))
R = pd.DataFrame(rows, columns=["r", "litho", "xr"])
R["m"] = pt2m(180 + R.r / Z)
R = R[(R.m > -0.5) & (R.m < 1130)].sort_values("m").reset_index(drop=True)
R["litho"] = R.litho.ffill().bfill()
def gclass(x):
    return next((n for xr, n in GS if x <= xr + 1.0), "boulder")
R["grain"] = [gclass(x) if np.isfinite(x) else None for x in R.xr]
R["grain"] = R.grain.ffill().bfill()
# run-length encode rows into beds (litho + grain class), absorbing <0.3 m slivers (bed-boundary lines)
beds = []
key = (R.litho + "|" + R.grain)
start = 0
for i in range(1, len(R) + 1):
    if i == len(R) or key[i] != key[start]:
        beds.append([R.m[start], R.m[i - 1], R.litho[start], R.grain[start], R.xr[start:i].median()]); start = i
B = pd.DataFrame(beds, columns=["base_m", "top_m", "litho", "grain", "width_x"])
B["top_m"] = list(B.base_m[1:]) + [B.top_m.iloc[-1]]
changed = True
while changed:                                  # merge slivers into the thicker neighbour
    changed = False
    th = B.top_m - B.base_m
    for i in th[th < 0.3].index:
        j = i - 1 if i > 0 and (i == len(B) - 1 or th[i - 1] >= th.get(i + 1, 0)) else i + 1
        if j not in B.index: continue
        B.loc[j, "base_m"] = min(B.base_m[i], B.base_m[j]); B.loc[j, "top_m"] = max(B.top_m[i], B.top_m[j])
        B = B.drop(i).sort_values("base_m").reset_index(drop=True); changed = True; break
# merge consecutive beds that became identical
out = []
for _, b in B.iterrows():
    if out and out[-1]["litho"] == b.litho and out[-1]["grain"] == b.grain: out[-1]["top_m"] = b.top_m
    else: out.append(b.to_dict())
B = pd.DataFrame(out); B[["base_m", "top_m"]] = B[["base_m", "top_m"]].round(2)
# the red fill is only used for red mudstones; sand-width red rows are outline/annotation artefacts
B.loc[(B.litho == "mudstone_red") & ~B.grain.isin(["clay", "silt"]), "litho"] = "sandstone"

# lithofacies codes (text column) -> features of the bed they label
W = page.get_text("words")
codes = [((w[1] + w[3]) / 2, w[4]) for w in W if 640 < w[0] < 745 and "_" in w[4]]
B["codes"] = ""
for y, c in codes:
    m = pt2m(y); i = B.index[(B.base_m <= m) & (B.top_m > m)]
    if len(i): B.loc[i[0], "codes"] = (B.loc[i[0], "codes"] + " " + c).strip()
TOK = {"planar_lamination": r"Ppl", "mottling": r"var|mot", "epsilon_xstrat": r"Epsi", "planar_xbed": r"Xpt|Xp\b|_Xp",
       "trough_xbed": r"Xt", "pebbles": r"Cms|Clen", "concretion": r"CaNod", "ripples": r"uCr",
       "fish": r"Fish", "red": r"Red"}
for k, rx in TOK.items(): B[k] = B.codes.str.contains(rx).astype(int)
B["grading"] = np.where(B.codes.str.contains(r"_nG"), "FU", np.where(B.codes.str.contains(r"_iG"), "CU", ""))
B.to_csv(os.path.join(HERE, "laventa_cs_beds.csv"), index=False)
print(f"beds: {len(B)}  0-{B.top_m.max():.0f} m;", B.litho.value_counts().to_dict())

# ---------------- s15: reinterpreted polarity (first magnetostratigraphy column) ----------------
pg = fitz.open(os.path.join(APP, "g2023v45a6_s15.pdf"))[0]
pm = pg.get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
G = np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.h, pm.w, 3).mean(-1)
# calibration from the thickness axis of s15 at 300 dpi (40 m -> y 2257, 560 m -> y 1313)
def py2m(y): return (2329.6 - y) / 1.8154
col = G[:, 225] < 100; tr = []; prev = None
for y in range(int(2329.6 + 3), int(2329.6 - 1130 * 1.8154), -1):
    if col[y] != prev: tr.append((py2m(y), "N" if col[y] else "R")); prev = col[y]
P = [(round(tr[i][0], 1), round(tr[i + 1][0], 1), tr[i][1]) for i in range(len(tr) - 1)]
P = [p for p in P if p[1] - p[0] > 1.0 and not (p[2] == "R" and p[1] < 10)]
P[0] = (0.0,) + P[0][1:]
CHRONS = ["C5ABn", "C5AAr", "C5AAn", "C5Ar.3r", "C5Ar.2n", "C5Ar.2r", "C5Ar.1n", "C5Ar.1r", "C5An.2n", "C5An.1r", "C5An.1n"]
PP = pd.DataFrame(P, columns=["base_m", "top_m", "polarity"]); PP["chron"] = CHRONS[:len(PP)]
PP.to_csv(os.path.join(HERE, "laventa_cs_polarity.csv"), index=False)
print(PP.to_string(index=False))
