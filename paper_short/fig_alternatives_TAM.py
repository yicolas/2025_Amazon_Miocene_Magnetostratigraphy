#!/usr/bin/env python3
"""TAM version of the alternatives figure: what age is the Purussaurus horizon
under each candidate chron? TAM's R is truncated above (Group 1 lignite,
paramagnetic, unrecoverable), so TAM is drawn at its slowest allowed rate
(R zone just fills the chron); faster rates push the horizon toward the chron base.
Two outputs: compacted (as measured) and decompacted (mud x2.6, lignite/peat x5)."""
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
# name, T, B, verdict, reason (verdicts inherited from the joint SC+TAM exclusion)
CAND = [
 ("C5Ar.2r",12.770,12.829,"ex","excluded at Santa Carina — pollen floor (12.7 Ma)"),
 ("C5Ar.3r",12.887,13.032,"ex","excluded at Santa Carina — pollen floor (12.7 Ma)"),
 ("C5AAr", 13.183,13.363,"ok","SURVIVES"),
 ("C5ABr", 13.608,13.739,"ex","excluded — Santa Carina base (14.08 Ma) predates the $\\it{Purussaurus}$ range"),
 ("C5ACr", 14.070,14.163,"ex","excluded — horizon older than the $\\it{Purussaurus}$ range (and pollen ceiling)"),
]
TAM_R, TAM_NLO, GAP = 50.0, 48.5, 17.0
PALYNO=(12.7,14.2); FOSSIL=13.8
C_PAL="#2166ac"; C_FOS="#d55e00"; C_OK="#1b7837"; C_BAD="#b2182b"; C_PUR="#6a3d9a"

for MODE in ("compacted","decompacted"):
    lig_stretch = 1.0 if MODE=="compacted" else 5.0/2.6   # lignite gap restored x5, mud rate x2.6
    rate_note = (lambda r: f"rate $\\geq${r:.0f} cm/Myr (as measured)") if MODE=="compacted" \
           else (lambda r: f"rate $\\geq${r*2.6:.0f} cm/Myr (decompacted $\\times$2.6; lignite $\\times$5)")
    fig,ax=plt.subplots(figsize=(11.2,7.6))
    X0,X1=12.30,14.72
    ax.set_xlim(X0,X1); ax.set_ylim(-1.80,9.4); ax.axis("off")
    ax.axvspan(PALYNO[0],FOSSIL,ymin=0.17,ymax=0.97,color=C_OK,alpha=0.07,zorder=0)
    for x,c,lab in [(PALYNO[0],C_PAL,"pollen floor\n12.7 Ma"),
                    (PALYNO[1],C_PAL,"pollen ceiling\n14.2 Ma"),
                    (FOSSIL,C_FOS,"$\\it{Purussaurus}$ ceiling\n13.8 Ma")]:
        ax.axvline(x,color=c,lw=1.3,ls=(0,(5,3)),ymin=0.18,ymax=0.90,zorder=2)
        ax.text(x,8.72,lab,color=c,fontsize=8.3,ha="center",va="bottom",weight="bold")
    ys,yh=7.55,0.60
    for nm,t,b,pol in GPTS:
        if b<X0 or t>X1: continue
        fc="0.15" if pol=="N" else "white"
        ax.add_patch(mpatches.Rectangle((max(t,X0),ys),min(b,X1)-max(t,X0),yh,facecolor=fc,edgecolor="0.15",lw=0.7,zorder=3))
        if nm and pol=="N" and (b-t)>0.17 and t>X0:
            ax.text((t+b)/2,ys+yh/2,nm,fontsize=7.0,ha="center",va="center",color="white",zorder=4)
        elif nm and pol=="N" and t>X0:
            ax.text((t+b)/2,ys+yh+0.14,nm,fontsize=6.8,ha="center",va="bottom",color="0.15")
        elif nm and pol=="R":
            ax.text((t+b)/2,ys-0.26,nm,fontsize=7.2,ha="center",va="top",color="0.35")
    ax.text(X0-0.015,ys+yh/2,"GPTS\n(ATNTS2012)",fontsize=8.2,ha="right",va="center",weight="bold")
    for a in np.arange(12.4,14.61,0.2):
        ax.plot([a,a],[ys+yh,ys+yh+0.06],color="0.5",lw=0.7)
        ax.text(a,ys+yh+0.36,f"{a:.1f}",fontsize=7,color="0.5",ha="center")
    ax.text(X1,ys+yh+0.36,"Ma (older $\\rightarrow$)",fontsize=7,color="0.5",ha="left")

    rh=0.42
    for i,(nm,T,B,verdict,reason) in enumerate(CAND):
        y=6.05-1.42*i
        D=B-T; r_min=TAM_R/D
        ok = verdict=="ok"
        if ok:
            ax.add_patch(mpatches.Rectangle((X0+0.01,y-0.50),X1-X0-0.02,1.14,facecolor=C_OK,
                         alpha=0.09,edgecolor=C_OK,lw=1.1,zorder=1))
        # R zone (white) fills the chron at the slowest allowed rate; N (black) below
        n_span=(TAM_NLO/TAM_R)*D
        ax.add_patch(mpatches.Rectangle((T,y),D,rh,facecolor="white",edgecolor="0.15",lw=0.9,zorder=4))
        ax.add_patch(mpatches.Rectangle((B,y),n_span,rh,facecolor="0.15",edgecolor="0.15",lw=0.9,zorder=4))
        # record lost above (Group 1 lignite): dotted open box
        lost=0.12
        ax.add_patch(mpatches.Rectangle((T-lost,y),lost,rh,facecolor="none",edgecolor="0.45",
                     ls=(0,(2,2)),lw=1.1,zorder=4))
        if i==0:
            ax.text(T-lost-0.03,y+rh/2,"record lost above\n(lignite, paramagnetic)",fontsize=6.8,
                    color="0.45",ha="right",va="center")
        # Purussaurus horizon bracket: [T - GAP*stretch/r_min, B]
        h_y = T - (GAP*lig_stretch)/r_min
        ax.plot([h_y,B],[y-0.22,y-0.22],color=C_PUR,lw=2.2,solid_capstyle="butt",zorder=6)
        ax.plot([h_y,h_y],[y-0.30,y-0.14],color=C_PUR,lw=1.6); ax.plot([B,B],[y-0.30,y-0.14],color=C_PUR,lw=1.6)
        ax.scatter([(h_y+B)/2],[y-0.22],marker="*",s=110,color=C_PUR,zorder=7)
        ax.text(B+0.02,y-0.22,f"$\\it{{Purussaurus}}$ {h_y:.2f}–{B:.2f} Ma",fontsize=7.6,
                color=C_PUR,va="center",ha="left")
        mark = "✓" if ok else "✗"; mc = C_OK if ok else C_BAD
        ax.text(X0+0.03,y+rh+0.13,f"{mark}  R zone = {nm}",fontsize=9.6,color=mc,weight="bold",va="bottom")
        ax.text(X0+0.03,y-0.44,rate_note(r_min),fontsize=7.4,color="0.35",va="top")
        txt = ("SURVIVES — the horizon dates to "+f"{h_y:.2f}–{B:.2f} Ma") if ok else reason
        ax.text(X1-0.03,y+rh+0.13,txt,fontsize=7.9,color=mc,ha="right",va="bottom",
                weight="bold" if ok else "normal")

    ax.text(X0+0.01,-0.62,
     "Tamshiyacu preserves only R-over-N: the upper normal was lost in the paramagnetic $\\it{Purussaurus}$-bearing lignite (Group 1), so each row pins\n"
     "the N$\\rightarrow$R boundary at the candidate chron base and draws the section at the slowest rate that fits the measured reversed zone (50 cm) inside\n"
     "the chron; faster rates shift the horizon (star) toward the chron base. Exclusions carry over from the Santa Carina N–R–N test (Fig. 3).",
     fontsize=8.6,va="top",color="0.2")
    ly=-1.52
    ax.add_patch(mpatches.Rectangle((X0+0.02,ly),0.12,0.22,facecolor="0.15",edgecolor="0.15")); ax.text(X0+0.18,ly+0.11,"normal polarity",fontsize=7.8,va="center")
    ax.add_patch(mpatches.Rectangle((X0+0.78,ly),0.12,0.22,facecolor="white",edgecolor="0.15")); ax.text(X0+0.94,ly+0.11,"reversed polarity",fontsize=7.8,va="center")
    ax.scatter([X0+1.62],[ly+0.11],marker="*",s=110,color=C_PUR); ax.text(X0+1.72,ly+0.11,"$\\it{Purussaurus}$ cranial horizon (age bracket)",fontsize=7.8,va="center")
    ttl_mode = "as measured (compacted)" if MODE=="compacted" else "decompacted (mud $\\times$2.6, lignite $\\times$5)"
    ax.set_title(f"Tamshiyacu: the $\\it{{Purussaurus}}$ horizon under each candidate chron — {ttl_mode}",
                 fontsize=11.5,weight="bold",pad=12)
    fig.tight_layout()
    out=__file__.rsplit("/",1)[0]+f"/figs/fig_alternatives_TAM_{MODE}.png"
    fig.savefig(out,dpi=250,bbox_inches="tight"); plt.close(fig)
    print("written",out.rsplit("/",1)[1])
