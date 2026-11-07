#!/usr/bin/env python3
"""Parse the 27 clean XRF Results.txt into one tidy CSV (run on device).
Usage: python3 parse_xrf.py <XRF_root> <out_csv>
- uses Results.txt only (result.txt are stale/mislabeled copies)
- header-aware per file; keeps a fixed element set where present
- validity==1 kept as 'valid'; mm->cm; group tag; dedupe exact (group,pos) dups
"""
import glob, os, sys, csv

ROOT, OUT = sys.argv[1], sys.argv[2]
ELEMENTS = ["Al","Si","P","S","Cl","K","Ca","Ti","V","Cr","Mn","Fe","Ni","Cu","Zn",
            "Rb","Sr","Y","Zr","Nb","Ba","Pb"]
EXTRA = ["cps","MSE","Mo inc","Mo coh"]
# XRF folder GROUP -> paleomag group key
GRP = {"GROUP1":"TAM_1_2_3B","GROUP2":"TAM_3A_4_5E_5D_5C","GROUP3":"TAM_5B_5A_6_7",
       "GROUP4":"SC_1C_1B_1A_2_3C","GROUP5":"SC_3B_3A_4D_4C_4B_4A",
       "GROUP6":"SC_5_6_7C_7B_7A","GROUP7":"SC_8"}
SITE = lambda g: "TAM" if g.startswith("TAM") else "SC"

files = sorted(glob.glob(os.path.join(ROOT,"GROUP*","**","Results.txt"), recursive=True))
rows_out = []
seen = set()
for f in files:
    parts = f.replace(ROOT,"").strip("/").split("/")
    gfolder = parts[0]
    seg = os.path.dirname(f).split("/")[-1]
    grp = GRP.get(gfolder, gfolder)
    lines = open(f, encoding="latin-1").read().splitlines()
    if len(lines) < 4: continue
    hdr = lines[2].split("\t")
    idx = {name: (hdr.index(name) if name in hdr else None)
           for name in ["position (mm)","validity"] + ELEMENTS + EXTRA}
    for ln in lines[3:]:
        c = ln.split("\t")
        if idx["position (mm)"] is None or len(c) <= idx["position (mm)"]: continue
        try: pos_mm = float(c[idx["position (mm)"]])
        except: continue
        key = (grp, round(pos_mm,1))
        if key in seen:   # dedupe exact duplicate positions (RUN2 re-scans)
            continue
        seen.add(key)
        val = c[idx["validity"]] if idx["validity"] is not None and len(c)>idx["validity"] else ""
        rec = {"site":SITE(grp),"group":grp,"segment":seg,
               "pos_mm":pos_mm,"pos_cm":pos_mm/10.0,
               "valid": 1 if val.strip()=="1" else 0}
        for name in EXTRA+ELEMENTS:
            j = idx[name]
            v = ""
            if j is not None and len(c)>j:
                try: v = float(c[j])
                except: v = ""
            rec[name.replace(" ","_")] = v
        rows_out.append(rec)

cols = ["site","group","segment","pos_mm","pos_cm","valid","cps","MSE","Mo_inc","Mo_coh"] + ELEMENTS
with open(OUT,"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in rows_out: w.writerow(r)
nvalid = sum(1 for r in rows_out if r["valid"]==1)
print(f"parsed {len(files)} files -> {len(rows_out)} rows ({nvalid} valid) -> {OUT}")
by = {}
for r in rows_out: by.setdefault(r["group"],[0,0]); by[r["group"]][0]+=1; by[r["group"]][1]+= (r["valid"]==1)
for g,(n,v) in by.items(): print(f"  {g:22} {n:4} rows, {v:4} valid")
