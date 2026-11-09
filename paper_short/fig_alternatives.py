#!/usr/bin/env python3
"""Fig: which reversed chron did SC record? Candidate placements vs independent constraints."""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

GPTS = [
 ("C5An.1n",12.049,12.174,"N"), ("",12.174,12.272,"R"),
 ("C5An.2n",12.272,12.474,"N"), ("C5Ar.1r",12.474,12.735,"R"),
 ("",12.735,12.770,"N"), ("",12.770,12.829,"R"),
 ("",12.829,12.887,"N"), ("C5Ar.3r",12.887,13.032,"R"),
 ("C5AAn",13.032,13.183,"N"), ("C5AAr",13.183,13.363,"R"),
 ("C5ABn",13.363,13.608,"N"), ("C5ABr",13.608,13.739,"R"),
 ("C5ACn",13.739,14.070,"N"), ("C5ACr",14.070,14.163,"R"),
 ("C5ADn",14.163,14.609,"N"),
]
CAND = [
 ("C5Ar.2r", 12.770,12.829, 12.62,12.98, 729,123,"ex_palyno","section top 12.62 Ma — younger than the pollen floor (12.7 Ma)"),
 ("C5Ar.3r", 12.887,13.032, 12.53,13.41, 297,129,"ex_palyno","section top 12.53 Ma — younger than the pollen floor (12.7 Ma)"),
 ("C5AAr",  13.183,13.363, 12.74,13.84, 239, 96,"survives","fits pollen window and the $\\it{Purussaurus}$ range (base 13.84 $\\approx$ 13.8 Ma — marginal)"),
 ("C5ABr",  13.608,13.739, 13.29,14.08, 328, 43,"ex_fossil","best shape fit, but section base 14.08 Ma predates the $\\it{Purussaurus}$ range ($\\leq$13.8 Ma)"),
 ("C5ACr",  14.070,14.163, 13.84,14.41, 462,168,"ex_both","section base 14.41 Ma — too old for pollen (14.2 Ma) and $\\it{Purussaurus}$ (13.8 Ma)"),
]
PALYNO=(12.7,14.2); FOSSIL_CEIL=13.8
C_PAL="#2166ac"; C_FOS="#d55e00"; C_OK="#1b7837"; C_BAD="#b2182b"

fig,ax=plt.subplots(figsize=(11.2,7.6))
X0,X1=12.30,14.72
ax.set_xlim(X0,X1); ax.set_ylim(-1.80,9.4); ax.axis("off")

ax.axvspan(PALYNO[0],FOSSIL_CEIL,color=C_OK,alpha=0.07,zorder=0)
for x,c,lab in [(PALYNO[0],C_PAL,"pollen floor\n12.7 Ma"),
                (PALYNO[1],C_PAL,"pollen ceiling\n14.2 Ma"),
                (FOSSIL_CEIL,C_FOS,"$\\it{Purussaurus}$ ceiling\n13.8 Ma")]:
    ax.axvline(x,color=c,lw=1.3,ls=(0,(5,3)),ymin=0.06,ymax=0.90,zorder=2)
    ax.text(x,8.72,lab,color=c,fontsize=8.3,ha="center",va="bottom",weight="bold")

ys,yh=7.55,0.60
for nm,t,b,pol in GPTS:
    fc="0.15" if pol=="N" else "white"
    ax.add_patch(mpatches.Rectangle((t,ys),b-t,yh,facecolor=fc,edgecolor="0.15",lw=0.7,zorder=3))
    if nm and pol=="N" and (b-t)>0.17:
        ax.text((t+b)/2,ys+yh/2,nm,fontsize=7.0,ha="center",va="center",color="white",zorder=4)
    elif nm and pol=="N":
        ax.text((t+b)/2,ys+yh+0.14,nm,fontsize=6.8,ha="center",va="bottom",color="0.15")
    elif nm and pol=="R":
        ax.text((t+b)/2,ys-0.26,nm,fontsize=7.2,ha="center",va="top",color="0.35")
ax.text(X0-0.015,ys+yh/2,"GPTS\n(ATNTS2012)",fontsize=8.2,ha="right",va="center",weight="bold")
for a in np.arange(12.4,14.61,0.2):
    ax.plot([a,a],[ys+yh,ys+yh+0.06],color="0.5",lw=0.7)
    ax.text(a,ys+yh+0.36,f"{a:.1f}",fontsize=7,color="0.5",ha="center")
ax.text(X1,ys+yh+0.36,"Ma (older $\\rightarrow$)",fontsize=7,color="0.5",ha="left")

rh=0.42
for i,(nm,ct,cb,st,sb,rate,mis,verdict,reason) in enumerate(CAND):
    y=6.00-1.30*i
    ok = verdict=="survives"
    if ok:
        ax.add_patch(mpatches.Rectangle((X0+0.01,y-0.50),X1-X0-0.02,1.14,facecolor=C_OK,
                     alpha=0.09,edgecolor=C_OK,lw=1.1,zorder=1))
    for (xa,xb,pol) in [(st,ct,"N"),(ct,cb,"R"),(cb,sb,"N")]:
        fc="0.15" if pol=="N" else "white"
        ax.add_patch(mpatches.Rectangle((xa,y),xb-xa,rh,facecolor=fc,edgecolor="0.15",lw=0.9,zorder=4))
    viols=[]
    if st<PALYNO[0]: viols.append((st,PALYNO[0],C_PAL))
    if sb>PALYNO[1]: viols.append((PALYNO[1],sb,C_PAL))
    if sb>FOSSIL_CEIL: viols.append((FOSSIL_CEIL,min(sb,X1-0.01),C_FOS))
    for xa,xb,c in viols:
        ax.add_patch(mpatches.Rectangle((xa,y-0.07),xb-xa,rh+0.14,facecolor="none",
                     edgecolor=c,hatch="////",lw=1.4,zorder=5))
    mark = "✓" if ok else "✗"
    mc = C_OK if ok else C_BAD
    ax.text(X0+0.03,y+rh+0.13,f"{mark}  R zone = {nm}",fontsize=9.6,color=mc,weight="bold",va="bottom")
    ax.text(X0+0.03,y-0.13,f"implied rate {rate} cm/Myr $\\cdot$ N–R–N shape misfit {mis} cm"
            +("  (best shape)" if mis==43 else ""),fontsize=7.4,color="0.35",va="top")
    ax.text(X1-0.03,y+rh+0.13,("SURVIVES — " if ok else "excluded — ")+reason,
            fontsize=7.9,color=mc,ha="right",va="bottom",weight="bold" if ok else "normal")

ax.text(X0+0.01,-0.62,
 "Each row stretches the measured Santa Carina section (Normal–Reversed–Normal; 103 / 43 / 71 cm) so that its reversed zone occupies one candidate\n"
 "reversed chron of the GPTS. Independent evidence — the pollen zone window (T-15/T-16) and the known $\\it{Purussaurus}$ range — then excludes every\n"
 "candidate but C5AAr. Reversed chrons younger than C5Ar.2r, or much longer than the reversed zone (e.g. C5Ar.1r), fail these tests outright.",
 fontsize=8.6,va="top",color="0.2")
ly=-1.52
ax.add_patch(mpatches.Rectangle((X0+0.02,ly),0.12,0.22,facecolor="0.15",edgecolor="0.15")); ax.text(X0+0.18,ly+0.11,"normal polarity",fontsize=7.8,va="center")
ax.add_patch(mpatches.Rectangle((X0+0.78,ly),0.12,0.22,facecolor="white",edgecolor="0.15")); ax.text(X0+0.94,ly+0.11,"reversed polarity",fontsize=7.8,va="center")
ax.add_patch(mpatches.Rectangle((X0+1.60,ly),0.12,0.22,facecolor="none",edgecolor=C_BAD,hatch="////",lw=1.2)); ax.text(X0+1.76,ly+0.11,"violates a constraint",fontsize=7.8,va="center")

ax.set_title("Which reversed chron did Santa Carina record?  Candidate GPTS placements vs. independent age constraints",
             fontsize=11.5,weight="bold",pad=12)
fig.tight_layout()
fig.savefig(__file__.rsplit("/",1)[0]+"/figs/fig_alternatives.png",dpi=250,bbox_inches="tight")
print("written figs/fig_alternatives.png")
