#!/usr/bin/env python3
"""
Updated TAM/SC scores under the Purussaurus + palyno bracket
============================================================
Two scoring layers:
 (1) constraint filter  - which cross_parameter models survive the fossil+palyno bracket
 (2) pattern-fit score  - for a site with N-R-N, how well do the observed normal/reversed
                          thicknesses match each candidate reversed chron's duration ratios
                          (a proper magnetostratigraphic test; independent of the old score)
Output (./out4): updated_scores.csv
"""
import os, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"out4"); os.makedirs(OUT,exist_ok=True)
ZONES=pd.read_csv(os.path.join(HERE,"out","polarity_zones.csv"))

# GTS2020 normals
NORM={"C5Ar.1n":(12.735,12.770),"C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),
      "C5ABn":(13.363,13.608),"C5ACn":(13.739,14.070),"C5ADn":(14.163,14.609)}
sn=sorted(NORM.items(),key=lambda kv:kv[1])
# reversed chrons = gaps; each with (top,base) and flanking normal durations
REV=[]
for (na,(ta,ba)),(nb,(tb,bb)) in zip(sn[:-1],sn[1:]):
    REV.append(dict(name=f"{na}r", top=ba, base=tb,
                    Nabove=(na, ba-ta), Nbelow=(nb, bb-tb)))
# rename to conventional-ish
NAME={"C5Ar.1nr":"C5Ar(.1r/.2r)","C5Ar.2nr":"C5Ar.3r/C5AAr(g)","C5AAnr":"C5AAr","C5ABnr":"C5ABr","C5ACnr":"C5ACr"}

# constraints
PALYNO=(12.7,14.2); GENUS_CEIL=13.8; TIGHT_CEIL=13.0

def zone_thick(site):
    z=ZONES[ZONES.site==site].reset_index(drop=True)
    return z

# ---- SC has a clean N-R-N; measure the three thicknesses
zsc=zone_thick("SC")
# N above = first N block(s) above the R; R; N below = normal blocks below R
R=zsc[zsc.polarity=="R"].iloc[0]
Nab=zsc[(zsc.polarity=="N")&(zsc.bottom_cm<=R.top_cm)]
Nbe=zsc[(zsc.polarity=="N")&(zsc.top_cm>=R.bottom_cm)]
sc_Nabove=Nab.thickness_cm.sum(); sc_R=R.thickness_cm; sc_Nbelow=Nbe.thickness_cm.sum()
sc_dtop=R.top_cm; sc_dbase=R.bottom_cm; sc_bottom=zsc.bottom_cm.max()

print("SC observed pattern:  N_above=%.0f cm | R=%.0f cm | N_below=%.0f cm  (section 0-%.0f cm)"
      %(sc_Nabove,sc_R,sc_Nbelow,sc_bottom))
print("TAM: only R(top)-N preserved (top normal lost with invalid Group 1) -> under-constrained\n")

rows=[]
for rv in REV:
    if not (PALYNO[0] <= (rv["top"]+rv["base"])/2 <= PALYNO[1]):
        pass
    dur_R=rv["base"]-rv["top"]
    s=sc_R/dur_R                              # sed rate from R thickness
    Nab_pred=rv["Nabove"][1]*s; Nbe_pred=rv["Nbelow"][1]*s
    misfit=abs(sc_Nabove-Nab_pred)+abs(sc_Nbelow-Nbe_pred)
    top_age=rv["top"]-sc_dtop/s               # section top (youngest)
    bot_age=rv["base"]+(sc_bottom-sc_dbase)/s # section bottom (oldest)
    rows.append(dict(reversed_chron=rv["name"], age_top=round(rv["top"],3),age_base=round(rv["base"],3),
        sed_cm_Ma=round(s), pattern_misfit_cm=round(misfit),
        sec_top_Ma=round(top_age,2), sec_bot_Ma=round(bot_age,2),
        palyno_ok=(top_age>=PALYNO[0] and bot_age<=PALYNO[1]),
        fossil_genus_ok=(bot_age<=GENUS_CEIL), fossil_tight_ok=(bot_age<=TIGHT_CEIL)))
sc=pd.DataFrame(rows).sort_values("pattern_misfit_cm")
sc.insert(0,"site","SC")
print("SC — pattern-fit score per candidate reversed chron (lower misfit = better):")
print(sc.to_string(index=False))

# best
best=sc.iloc[0]
print(f"\n  >> SC best pattern fit: {best.reversed_chron} @ {best.age_top}-{best.age_base} Ma, "
      f"sed {best.sed_cm_Ma} cm/Ma, section {best.sec_top_Ma}-{best.sec_bot_Ma} Ma, misfit {best.pattern_misfit_cm} cm")

# ---- TAM: if it shares SC's best reversed chron, what sed rate does its R need?
ztam=zone_thick("TAM"); Rt=ztam[ztam.polarity=="R"].iloc[0]
tam_R=Rt.thickness_cm
for cand in ["C5ABr","C5AAr"]:
    rv=[r for r in REV if r["name"].startswith(cand[:4]) and r["name"]==cand]
    # find by name
    match=[r for r in REV if r["name"]==cand]
    if match:
        d=match[0]["base"]-match[0]["top"]; s=tam_R/d
        print(f"  TAM R={tam_R:.0f}cm as {cand} (dur {d:.3f} Ma) -> sed {s:.0f} cm/Ma "
              f"(top normal lost, so weak; consistency check only)")

sc.to_csv(os.path.join(OUT,"updated_scores.csv"),index=False)
print("\nwrote updated_scores.csv")
