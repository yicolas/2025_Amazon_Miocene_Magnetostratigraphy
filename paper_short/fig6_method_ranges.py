#!/usr/bin/env python3
"""Fig 6 (working): which chronometers reach which low-latitude Amazonian units.
C14/OSL/Be-10 date the terraces and young cover; only detrital-zircon U-Pb and
magnetostratigraphy reach the Miocene Pebas. Ranges approximate (draft)."""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({"font.family":"serif","font.size":9})
LEAD_BLUE="#4F81BD"; DARK="#345A8A"; RED="#c0392b"; GREY="0.45"

# method, min age (a), max age (a), color, note
METHODS=[
 ("Radiocarbon ($^{14}$C)",            3e2, 5.5e4, "#6baed6", "organics; ~55 ka ceiling"),
 ("OSL (quartz / feldspar pIRIR)",     1e2, 2e5,  "#74c476", "to ~0.2 Ma; pIRIR to ~0.7 Ma"),
 ("$^{10}$Be surface exposure",        1e3, 5e6,  "#fd8d3c", "low-latitude production low"),
 ("$^{26}$Al/$^{10}$Be burial",        2e5, 5e6,  "#e6550d", "burial of quartz-bearing fills"),
 ("Magnetostratigraphy (+ biostrat anchor)", 7.8e5, 3e7, "#756bb1", "needs reversals + an anchor"),
 ("Detrital-zircon U–Pb (MDA)",       1e6, 3e7,  "#31a354", "any age; ~1–2% precision"),
]
PEBAS=(13.18e6,13.36e6)

fig,ax=plt.subplots(figsize=(9.2,4.6))
X0,X1=1e2,3e7
ax.set_xscale("log"); ax.set_xlim(X0,X1); ax.set_ylim(-1.0,len(METHODS)+1.1)

ax.axvspan(*PEBAS,color=RED,alpha=0.30,zorder=1)
ax.text(PEBAS[1]*1.25,len(METHODS)+0.62,"Pebas Fm (this study)\n13.18–13.36 Ma",fontsize=7.6,
        ha="left",va="center",color=RED,weight="bold")

for i,(nm,a0,a1,c,note) in enumerate(METHODS):
    y=len(METHODS)-1-i
    reach = a1>=PEBAS[0]
    ax.barh(y,a1-a0,left=a0,height=0.52,color=c,alpha=0.9 if reach else 0.55,
            edgecolor="0.2",lw=0.6,zorder=3)
    ax.text(a0*0.88,y,nm,ha="right",va="center",fontsize=8.6,
            color=DARK if reach else "0.25",weight="bold" if reach else "normal")
    ax.text(a0*1.1, y-0.40, note, fontsize=7.0, color=GREY, va="top", ha="left", style="italic")
    if not reach:
        ax.plot([a1*1.5],[y],marker="x",color=RED,ms=7,mew=2,zorder=4)

ax.axvline(np.sqrt(PEBAS[0]*PEBAS[1]),color=RED,lw=1.0,ls=(0,(4,3)),zorder=2)
ax.set_yticks([])
ax.set_xlabel("Age (years before present, log scale)")
xt=[1e2,1e3,1e4,1e5,1e6,1e7]
ax.set_xticks(xt); ax.set_xticklabels(["100 a","1 ka","10 ka","100 ka","1 Ma","10 Ma"])
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
ax.set_title("Geochronometer age ranges vs. the middle Miocene Pebas target — why U–Pb + magnetostratigraphy",
             fontsize=10.5,color=DARK,weight="bold",loc="left",pad=10)
fig.text(0.06,0.0,"Working figure — method limits approximate; red $\\times$ marks where a method expires before the Pebas target (red band).",
        fontsize=7.4,color=GREY,va="bottom")
fig.tight_layout()
fig.savefig(__file__.rsplit("/",1)[0]+"/figs/fig6_method_ranges.png",dpi=250,bbox_inches="tight")
print("written figs/fig6_method_ranges.png")
