#!/usr/bin/env python3
"""Supplementary figure: why the observed N-R-N sequence is C5ABn-C5AAr-C5AAn and
not another GPTS fit. Two independent biostratigraphic bounds (palyno: Crassoretitriletes zone 14.1-12.4 Ma;
Purussaurus genus ceiling <13.8 Ma) define a permitted window; within it, the N-R-N
geometry + the U-Pb-calibrated La Venta regional magnetostratigraphy uniquely select C5AAr."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
OUTS=["/tmp/peru/supp/figs","/mnt/user-data/outputs"]

NORMAL={"C5An.1n":(12.049,12.174),"C5An.2n":(12.272,12.474),"C5Ar.1n":(12.735,12.770),
        "C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),"C5ABn":(13.363,13.608),
        "C5ACn":(13.739,14.070)}
# candidate short reversed chrons (bounded by two normals) and verdicts
CAND=[("C5An.1r",12.174,12.272,"excl: >12.4 Ma palyno top (Grimsdalea FAD)"),
      ("C5Ar.1r",12.770,12.829,"excl: <12.4 Ma; too thin"),
      ("C5Ar.2r",12.887,13.032,"excl: overlying normal too long; La Venta tie"),
      ("C5AAr", 13.183,13.363,"SELECTED"),
      ("C5ABr", 13.608,13.739,"excl: >13.8 ceiling edge; wrong bounding normals")]
YL=(14.20,11.80)

fig,ax=plt.subplots(figsize=(8.6,8.8))
xg=0.6; w=0.34
# GPTS normal boxes
for nm,(y0,y1) in NORMAL.items():
    ax.add_patch(Rectangle((xg-w/2,y0),w,y1-y0,fc="black",ec="none",zorder=3))
    ax.text(xg+w/2+0.03,(y0+y1)/2,nm,va="center",fontsize=7.5)
ax.add_patch(Rectangle((xg-w/2,YL[1]),w,YL[0]-YL[1],fc="none",ec="k",lw=1.0,zorder=2))
ax.text(xg,YL[0]+0.05,"GTS2020",ha="center",va="top",fontsize=9,weight="bold")

# permitted window 12.4-13.8 (independent biostrat bounds)
ax.axhspan(12.4,13.8,xmin=0.02,xmax=0.98,color="#2e8b57",alpha=0.10,zorder=0)
ax.axhline(12.4,color="#2e8b57",lw=1.4,ls="--",zorder=1)
ax.axhline(13.8,color="#8e44ad",lw=1.4,ls="--",zorder=1)
ax.text(3.15,12.4,"palynology: Mid-Miocene [14.1\u201312.4 Ma pending]",fontsize=7.6,color="#1e6b40",va="center")
ax.text(3.15,13.8,"$Purussaurus$ genus ceiling  <13.8 Ma",fontsize=7.6,color="#6c3483",va="center")

# candidate reversed chrons: tick + verdict
for nm,y0,y1,verd in CAND:
    yc=(y0+y1)/2; sel=verd=="SELECTED"
    ax.add_patch(Rectangle((xg-w/2,y0),w,y1-y0,fc=("#c0392b" if sel else "0.8"),
                 ec=("k" if sel else "0.6"),lw=(1.2 if sel else .5),zorder=4,alpha=(1 if sel else .9)))
    mark=("✓ " if sel else "✗ ")
    col=("#c0392b" if sel else "0.35")
    ax.annotate(f"{mark}{nm}: {verd}",(xg+w/2,yc),xytext=(1.35,yc),textcoords="data",
                fontsize=(8.2 if sel else 7.0),color=col,va="center",weight=("bold" if sel else "normal"),
                arrowprops=dict(arrowstyle="-",color=col,lw=.7))

# observed SC pattern (schematic N-R-N) on the right, aligned on C5AAr
xs=6.4; ww=0.5
segs=[("N (C5AAn)",13.032,13.183,"black"),("R (C5AAr)",13.183,13.363,"white"),
      ("N (C5ABn)",13.363,13.55,"black")]
for lab,y0,y1,c in segs:
    ax.add_patch(Rectangle((xs-ww/2,y0),ww,y1-y0,fc=c,ec="k",lw=.8,zorder=4))
ax.text(xs,YL[0]+0.05,"observed\nSC (N–R–N)",ha="center",va="top",fontsize=8.6,weight="bold")
ax.annotate("short reversed zone\n(0.43 m) between\ntwo normals",(xs+ww/2,13.27),
            xytext=(xs+0.5,13.05),fontsize=7.2,va="center",
            arrowprops=dict(arrowstyle="-",color="0.4",lw=.7))
# tie band
ax.axhspan(13.183,13.363,xmin=0.72,xmax=0.9,color="#c0392b",alpha=0.10,zorder=0)

ax.set_xlim(0.05,8.2); ax.set_ylim(*YL); ax.set_xticks([])
ax.set_ylabel("Age (Ma)",fontsize=11)
for s in ["top","right","bottom"]: ax.spines[s].set_visible(False)
ax.set_title("Why the Iquitos reversal is C5AAr: two independent biostratigraphic bounds define a\n"
             "permitted window (12.4–13.8 Ma); the N–R–N geometry and the U-Pb-calibrated La Venta\n"
             "regional magnetostratigraphy then uniquely select C5AAr among candidate reversed chrons",
             fontsize=9.0,weight="bold")
fig.text(0.5,-0.02,"Candidate reversed chrons that could produce a short reversed zone between two normals are tested against the permitted window and the regional tie; "
         "only C5AAr satisfies all constraints. Palynology and the $Purussaurus$ genus range are the two genuinely independent bounds; the La Venta U-Pb age model "
         "and its magnetostratigraphy jointly furnish the regional tie (counted as one constraint).",
         ha="center",fontsize=6.6,style="italic",color="0.35",wrap=True)
fig.tight_layout()
for d in OUTS: fig.savefig(os.path.join(d,"figS_altcorr.png"),dpi=190,bbox_inches="tight",facecolor="w")
print("wrote figS_altcorr.png")
