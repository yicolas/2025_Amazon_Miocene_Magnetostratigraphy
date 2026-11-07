import os
#!/usr/bin/env python3
"""Figure 4 (combined, 3 panels):
(a) Age-depth across chron C5AAr for both Iquitos sections (TAM, SC) vs marine
    reference gradients — the primary measured rates.
(b) Accumulation-rate compilation (forest plot) across the Pebas system.
(c) Sadler plot (accumulation rate vs measurement timespan) placing the Iquitos
    bay margin among global depositional settings; labels = setting / Author date.

Layout: (a) and (b) on the top row; (c) spans the full bottom row so its
two-line source labels stay legible."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from fig_sadler import draw_sadler
OUT=os.environ.get("PEBAS_OUT","/mnt/user-data/outputs")

# ---- shared age model (chron C5AAr) ----
C5AAr=(13.183,13.363); C5AAn=(13.032,13.183)
DATUM=13.363; TOP=13.183; dt=DATUM-TOP           # 0.180 Myr
SC_TH=0.43; TAM_TH=0.50
SC_R=SC_TH/dt; TAM_R=TAM_TH/dt                    # compacted m/Myr
C_TAM="#c0392b"; C_SC="#e67e22"; C_1237="#2166ac"; C_926="#1b7837"; C_1238="#8e6c1f"
C_BAS="#7f5539"; C_ORB="#b8860b"

fig=plt.figure(figsize=(13.6,12.8))
gs=gridspec.GridSpec(2,2,height_ratios=[1.0,1.30],hspace=0.30,wspace=0.26,
                     top=0.905,bottom=0.05,left=0.075,right=0.975)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,:])

# ================= (a) age-depth: TAM & SC =================
ax=axA; amin,amax=13.10,13.46
ax.axvspan(*C5AAr,color="#c0392b",alpha=0.08,zorder=0)
ax.axvspan(max(C5AAn[0],amin),C5AAn[1],color="0.88",alpha=0.6,zorder=0)
ax.text(np.mean(C5AAr),3.02,"C5AAr (reversed)",ha="center",fontsize=8,color="#a93226")
ax.text((C5AAn[1]+amin)/2,3.02,"C5AAn",ha="center",fontsize=8,color="0.45")
def grad(rate,c,lab):
    a=np.array([DATUM,TOP]); ax.plot(a,rate*(DATUM-a),color=c,lw=1.8,label=lab)
    ax.annotate(f"{rate:g}",(TOP,rate*dt),color=c,fontsize=8,weight="bold",va="center",ha="left",
                xytext=(3,0),textcoords="offset points")
grad(16.0,C_926,"ODP 926 Ceará Rise (~16, astro)")
grad(7.7,C_1237,"ODP 1237 Nazca Ridge (~7.7)")
a=np.array([DATUM,TOP])
ax.fill_between(a,SC_R*2.2*(DATUM-a),SC_R*3.0*(DATUM-a),color=C_SC,alpha=0.15,zorder=2)
ax.plot(a,SC_R*(DATUM-a),color=C_SC,lw=2.6,zorder=3,label=f"SC (this study, {SC_R:.1f})")
ax.errorbar([DATUM,TOP],[0,SC_R*dt],xerr=0.008,yerr=0.03,fmt="o",ms=7,color=C_SC,mec="k",
            ecolor="k",elinewidth=0.8,capsize=2,zorder=5)
ax.plot(a,TAM_R*(DATUM-a),color=C_TAM,lw=2.2,ls="--",zorder=3,label=f"TAM (this study, ≥{TAM_R:.1f}, truncated)")
ax.plot(DATUM,0,"o",ms=7,color=C_TAM,mec="k",zorder=5)
ax.plot(TOP,TAM_R*dt,">",ms=8,color=C_TAM,mec="k",zorder=5)
ax.annotate("shared C5AAr base\n(N→R datum)",(DATUM,0),fontsize=7,color="0.25",
            xytext=(-4,12),textcoords="offset points",ha="right")
ax.axhline(0,color="k",lw=0.5,ls=":")
ax.set_xlim(amax,amin); ax.set_ylim(-0.15,3.25)
ax.set_xlabel("Age (Ma)"); ax.set_ylabel("Thickness above C5AAr base (m)")
ax.set_title("(a) Age–depth across C5AAr — TAM vs SC vs marine references",fontsize=9.3,weight="bold",loc="left")
ax.legend(fontsize=6.9,loc="upper left",frameon=False)
ax.text(13.15,1.7,"TAM ≈ SC: both\n~2–3 m/Myr (compacted),\n~6–9 decompacted",fontsize=6.6,color="#8a3b1e",
        ha="center",va="center")

# ================= (b) rate compilation =================
ax=axB
rows=[ # label, lo, hi, color, group
 ("TAM — this study (Iquitos bay margin)\ncompacted ≥2.8; decompacted ~6–9",2.8,9,C_TAM,"this study"),
 ("SC — this study (Iquitos bay margin)\ncompacted 2.4; decompacted ~5–7",2.4,7,C_SC,"this study"),
 ("ODP 1237 Nazca Ridge (SE Pacific)\nMix et al. 2003",7.7,7.7,C_1237,"marine ref"),
 ("ODP 926 Ceará Rise (equat. Atlantic)\nBackman & Raffi 1997",16,16,C_926,"marine ref"),
 ("ODP 1238 Carnegie Ridge (off Ecuador)\nMix et al. 2003",10,72,C_1238,"marine ref"),
 ("Pebas orbital (obliquity) cyclicity\nHoorn et al. 2022",100,100,C_ORB,"Pebas basin"),
 ("Pebas thickness/duration (350–1074 m)\nWesselingh et al. 2006",45,130,C_BAS,"Pebas basin"),
]
y=np.arange(len(rows))[::-1]
for yi,(lab,lo,hi,c,g) in zip(y,rows):
    if hi>lo: ax.plot([lo,hi],[yi,yi],color=c,lw=6.5,alpha=0.45,solid_capstyle="round")
    ax.plot([lo,hi],[yi,yi],color=c,lw=1.5)
    ax.scatter([np.sqrt(lo*hi)],[yi],color=c,s=30,zorder=5,ec="k",lw=0.4)
    ax.text(160,yi,lab,fontsize=6.5,va="center",ha="left",color=c)
ax.axvspan(2.4,9,color="#c0392b",alpha=0.06,zorder=0)
ax.set_xscale("log"); ax.set_xlim(1.8,3500); ax.set_ylim(-1.0,len(rows)-0.4)
ax.set_yticks([]); ax.set_xlabel("Linear sedimentation rate (m/Myr, log)")
ax.set_title("(b) Rate varies >1 order of magnitude across the Pebas system",fontsize=9.0,weight="bold",loc="left")
ax.spines[["top","right","left"]].set_visible(False)
ax.annotate("",xy=(2.4,-0.55),xytext=(130,-0.55),arrowprops=dict(arrowstyle="<->",color="0.5",lw=0.9))
ax.text(18,-0.80,"~1.5 orders of magnitude",fontsize=6.4,color="0.5",ha="center",va="top")

# ================= (c) Sadler plot =================
draw_sadler(axC, title=False, fs_lab=5.8)
axC.set_title("(c) Sadler plot — accumulation rate vs measurement timespan (labels: setting / source)",
              fontsize=9.6,weight="bold",loc="left")

fig.suptitle("Figure 4.  Middle-Miocene accumulation across the Pan-Amazonian Pebas system:\n"
             "measured Iquitos rates (a,b) and their timescale context among global settings (c)",
             fontsize=12.2,weight="bold",y=0.978)
fig.savefig(f"{OUT}/fig4_combined.png",dpi=185,bbox_inches="tight",facecolor="w")
fig.savefig(f"{OUT}/fig4_combined.pdf",bbox_inches="tight")
print("wrote fig4_combined.png  (SC=%.2f TAM=%.2f m/Myr compacted)"%(SC_R,TAM_R))
