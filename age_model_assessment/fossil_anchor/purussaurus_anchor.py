#!/usr/bin/env python3
"""
Purussaurus fossil anchor -> age bracket on the TAM column
==========================================================
Combines a vertebrate (Purussaurus) CEILING with the palyno T15/T16 FLOOR to
bracket the section, and tests which mid-Miocene reversed chrons survive the
bracket vs the best-fit inclination age models.

Output (./out3): fig_purussaurus_bracket.png, surviving_chrons.csv
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"out3"); os.makedirs(OUT,exist_ok=True)

# GTS2020 normals (Ma)
NORM={"C5An.1n":(12.049,12.174),"C5An.2n":(12.272,12.474),"C5Ar.1n":(12.735,12.770),
 "C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),"C5ABn":(13.363,13.608),
 "C5ACn":(13.739,14.070),"C5ADn":(14.163,14.609)}
# reversed intervals (gaps between successive normals) + name (reversed chron below each normal)
sn=sorted(NORM.items(),key=lambda kv:kv[1])
REV=[]
for (n0,(t0,b0)),(n1,(t1,b1)) in zip(sn[:-1],sn[1:]):
    REV.append((b0,t1,f"{n0}r"))   # base of younger normal .. top of next-older normal

# constraints
PALYNO_FLOOR = 12.7            # T15/T16 boundary: section is OLDER than this (>=)
PURU_GENUS_CEIL = 13.8         # genus Purussaurus FAD ~13.8 Ma (Laventan base): section YOUNGER than this
PURU_NEIV = (11.8,13.8)        # P. neivensis species range (Laventan)
PEBAS_ASSEMBLAGE_CEIL = 13.0   # tight, context-based (Iquitos Pebas croc bonebeds ~13 Ma)

# best-fit inclination reversal ages (from age_scenarios), coeval family
MODEL_R = {"TAM (778 cm/Ma)":(12.935,13.00),"SC (786 cm/Ma)":(13.08,13.135)}

def survives(lo,hi):  # reversed intervals overlapping [lo,hi]
    return [(a,b,nm) for (a,b,nm) in REV if b>lo and a<hi]

brackets={
 "palyno only (>=12.7)": (PALYNO_FLOOR, 14.7),
 "palyno + Purussaurus genus (12.7-13.8)": (PALYNO_FLOOR, PURU_GENUS_CEIL),
 "palyno + Pebas-assemblage ceiling (12.7-13.0)": (PALYNO_FLOOR, PEBAS_ASSEMBLAGE_CEIL),
}
rows=[]
for name,(lo,hi) in brackets.items():
    s=survives(lo,hi)
    rows.append(dict(bracket=name,lo=lo,hi=hi,n_candidates=len(s),
                     candidates="; ".join(f"{nm}({a}-{b})" for a,b,nm in s)))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"surviving_chrons.csv"),index=False)
print("SURVIVING REVERSED-CHRON CANDIDATES PER BRACKET")
for r in rows:
    print(f"\n{r['bracket']}  -> {r['n_candidates']} candidates")
    print("   ",r['candidates'])
print("\nBest-fit model reversal ages:")
for k,(a,b) in MODEL_R.items(): print(f"   {k}: {a}-{b} Ma")

# ---------------------------------------------------------------- figure
fig,ax=plt.subplots(figsize=(9,6))
lo,hi=11.9,14.7
# GPTS column (normals black)
xg=0.0
for nm,(t,b) in NORM.items():
    ax.axhspan(t,b,xmin=0.0,xmax=0.12,color="k")
    ax.text(-0.01,(t+b)/2,nm,ha="right",va="center",fontsize=7)
ax.text(0.06,hi+0.05,"GPTS\n(GTS2020)",ha="center",va="bottom",fontsize=8)
# reversed candidates labels in middle
for a,b,nm in REV:
    if b>lo and a<hi:
        ax.text(0.16,(a+b)/2,nm,fontsize=6,va="center",color="0.4")

# constraint bars
ax.axhspan(lo,PALYNO_FLOOR,xmin=0.30,xmax=0.42,color="#c0392b",alpha=0.25)
ax.axhline(PALYNO_FLOOR,xmin=0.30,xmax=0.42,color="#c0392b",lw=2)
ax.text(0.36,PALYNO_FLOOR-0.02,"palyno T15/T16\nfloor 12.7 Ma\n(older than)",ha="center",va="top",fontsize=7,color="#c0392b")

ax.axhspan(PURU_GENUS_CEIL,hi,xmin=0.46,xmax=0.58,color="#27ae60",alpha=0.20)
ax.axhline(PURU_GENUS_CEIL,xmin=0.46,xmax=0.58,color="#27ae60",lw=2)
ax.text(0.52,PURU_GENUS_CEIL+0.02,"Purussaurus genus\nFAD ceiling ~13.8 Ma\n(younger than)",ha="center",va="bottom",fontsize=7,color="#1e8449")
ax.axhline(PEBAS_ASSEMBLAGE_CEIL,xmin=0.46,xmax=0.58,color="#27ae60",lw=1,ls=":")
ax.text(0.59,PEBAS_ASSEMBLAGE_CEIL,"Pebas assemblage\n~13.0 (tight)",ha="left",va="center",fontsize=6,color="#1e8449")

# allowed bracket band (12.7-13.8)
ax.axhspan(PALYNO_FLOOR,PURU_GENUS_CEIL,xmin=0.64,xmax=0.72,color="#f1c40f",alpha=0.35)
ax.text(0.68,(PALYNO_FLOOR+PURU_GENUS_CEIL)/2,"ALLOWED\n12.7-13.8 Ma",ha="center",va="center",fontsize=7,rotation=90)

# best-fit model reversals
cols={"TAM (778 cm/Ma)":"#2980b9","SC (786 cm/Ma)":"#8e44ad"}
x0=0.78
for i,(k,(a,b)) in enumerate(MODEL_R.items()):
    xx=x0+i*0.10
    ax.axhspan(a,b,xmin=xx,xmax=xx+0.08,color=cols[k],alpha=0.9)
    ax.text(xx+0.04,a-0.02,k.split()[0]+"\nreversal",ha="center",va="top",fontsize=6,color=cols[k])

ax.set_ylim(hi,lo); ax.set_xlim(-0.08,1.0); ax.set_xticks([])
ax.set_ylabel("Age (Ma)")
ax.set_title("Purussaurus (ceiling) + palynology (floor) bracket the TAM column to 12.7-13.8 Ma\n"
             "Best-fit magnetostratigraphic reversal (~12.9-13.0 Ma) => C5AAr",fontsize=10)
plt.tight_layout()
p=os.path.join(OUT,"fig_purussaurus_bracket.png"); fig.savefig(p,dpi=200,bbox_inches="tight"); plt.close()
print("\nwrote",p,"and surviving_chrons.csv")
