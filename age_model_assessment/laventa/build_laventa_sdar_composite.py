#!/usr/bin/env python3
"""Stack the Mora-Rojas et al. (2023) SDAR bed logs (GSA suppl. g2023v45a6_s3) for the
La Venta polygon (LV_1..LV_6) into one composite section, and tag each bed with the
accessory features named in its notes / grading columns.

Composite metres: LV_n base = top of LV_{n-1} (LV_1 base = 0; the Saldana tuff at -2..0 is dropped).
Check: LV_3 'crocodile bone - 280 m', 'turtle shell bone - 287 m' and LV_4 'Toxodont teeth - 333.5 m'
fall inside the tagged beds with these offsets, so the stacking reproduces the authors' metres.

usage:  python build_laventa_sdar_composite.py "<path to>/g2023v45a6_s3.zip"
writes: laventa_sdar_composite.csv (next to this script)
"""
import sys, os, re, io, zipfile, warnings
import pandas as pd
warnings.filterwarnings("ignore")

ZIP = sys.argv[1] if len(sys.argv) > 1 else "g2023v45a6_s3.zip"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "laventa_sdar_composite.csv")
SECTIONS = ["1_LV_1", "2_LV_2", "3_LV_3", "4_LV_4", "5_LV_5", "6_LV_6"]

FEATURES = {   # feature -> regex on the bed notes
    "trough_xbed":      r"trough",
    "planar_xbed":      r"cross-bedding tangential|planar cross-bedding|tabular cross|cross-bedding indicated|cross bedding tangential",
    "epsilon_xstrat":   r"epsilon",
    "planar_lamination":r"parallel lamination|planar lamination|parallel stratification|planar-parallel",
    "pebbles":          r"conglomeratic|conglomerate lens",
    "concretion":       r"nodule",
    "mottling":         r"mottl|varicolou?red",
    "bone":             r"\bbone|teeth|tooth",
    "plant":            r"plant remains|leaves",
    "tuff_ash":         r"\bash\b",
}
FOSSIL = r"(crocodile|turtle|toxodont|rodent)[^-–−]*[-–−]\s*([\d.]+)\s*m"

rows = []; base = 0.0
with zipfile.ZipFile(ZIP) as z:
    for sec in SECTIONS:
        name = next(n for n in z.namelist() if n.endswith(f"{sec}_section_SDAR.xlsx"))
        b = pd.read_excel(io.BytesIO(z.read(name)), sheet_name="beds")
        b = b[b.top > 0].copy()                       # drop the pre-Honda tuff below LV_1 datum
        for _, r in b.iterrows():
            note = str(r.notes) if pd.notna(r.notes) else ""
            g = str(r.grading).lower() if pd.notna(r.grading) else ""
            row = dict(section=sec.split("_", 1)[1], base_m=round(base + max(r.base, 0), 2),
                       top_m=round(base + r.top, 2), litho=r.prim_litho, grain_size=r.grain_size,
                       grading={"normal": "FU", "inverse": "CU"}.get(g, ""), notes=note)
            for k, rx in FEATURES.items():
                row[k] = int(bool(re.search(rx, note, re.I)))
            m = re.search(FOSSIL, note, re.I)
            row["fossil"] = f"{m.group(1).lower()}@{m.group(2)}" if m else ""
            if r.prim_litho == "tuff": row["tuff_ash"] = 1
            rows.append(row)
        base += b.top.max()
df = pd.DataFrame(rows); df.to_csv(OUT, index=False)
print(f"wrote {OUT}: {len(df)} beds, 0-{df.top_m.max():.1f} m")
for _, r in df[df.fossil != ""].iterrows():
    print(f"  {r.fossil:16s} in bed {r.base_m:.1f}-{r.top_m:.1f} m ({r.section})")
