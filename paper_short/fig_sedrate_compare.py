#!/usr/bin/env python3
"""Sedimentation-rate comparison across the middle Miocene.
(a) age-depth across the C5AA interval: Pebas vs ODP 1237 vs ODP 926 Ceara Rise.
(b) accumulation-rate context: this study vs open-marine references vs Pebas basin-scale
    estimates, showing the Iquitos margin is a starved sub-environment (gives the low rate credence)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
OUT="/mnt/user-data/outputs"
C5AAr=(13.183,13.363); NORM={"C5AAn":(13.032,13.183),"C5ABn":(13.363,13.608)}; DATUM=13.363

fig,(axL,axR)=plt.subplots(1,2,figsize=(12,5.2),gridspec_kw=dict(width_ratios=[1.15,1.0],wspace=0.28))

# ---- (a) age-depth ----
amin,amax=12.85,13.70
axL.axvspan(*C5AAr,color="#c0392b",alpha=0.10,zorder=0)
for nm,(t,b) in NORM.items(): axL.axvspan(max(t,amin),min(b,amax),color="0.86",alpha=0.5,zorder=0)
axL.text(np.mean(C5AAr),8.6,"C5AAr",ha="center",fontsize=8,color="#c0392b")
axL.text(13.10,0.25,"C5AAn",fontsize=7.5,color="0.4",ha="center"); axL.text(13.50,0.25,"C5ABn",fontsize=7.5,color="0.4",ha="center")
ages=np.linspace(amin,amax,50)
for nm,rate,c in [("Pebas SC/TAM (this study, decompacted ~7)",7.0,"#c0392b"),
                  ("ODP 1237 Nazca Ridge (~7.7)",7.7,"#2166ac"),
                  ("ODP 926 Ceara Rise (~16, astro-tuned)",16.2,"#1b7837")]:
    axL.plot(ages,rate*(DATUM-ages),color=c,lw=2,label=nm)
    axL.annotate(f"{rate:.0f}",(amin,rate*(DATUM-amin)),color=c,fontsize=8,va="center",ha="left",weight="bold")
axL.axhline(0,color="k",lw=0.6,ls=":"); axL.text(DATUM,-0.5,"C5AAr base (13.363 Ma)",fontsize=6.5,ha="center",va="top",color="0.3")
axL.set_xlim(amax,amin); axL.set_ylim(-1.3,9)
axL.set_xlabel("Age (Ma)"); axL.set_ylabel("Thickness above C5AAr base (m)")
axL.set_title("(a) Age–depth across C5AA",fontsize=9.5,weight="bold")
axL.legend(fontsize=6.8,loc="upper right",frameon=False); axL.spines[["top","right"]].set_visible(False)

# ---- (b) accumulation-rate context (log scale) ----
rows=[  # label, low, high (m/Myr), color, group
 ("Iquitos margin — this study\n(compacted 3; decompacted 6–9)",3,9,"#c0392b","site"),
 ("ODP 1237 Nazca Ridge (Pacific)",7.7,7.7,"#2166ac","marine"),
 ("ODP 926 Ceara Rise (Atlantic)",16,16,"#1b7837","marine"),
 ("Pebas orbital cyclicity\n(Hoorn et al. 2022)",100,100,"#7f5539","basin"),
 ("Pebas thickness/duration avg\n(Hoorn et al. 2010)",45,130,"#7f5539","basin"),
]
y=np.arange(len(rows))[::-1]
for yi,(lab,lo,hi,c,g) in zip(y,rows):
    if hi>lo: axR.plot([lo,hi],[yi,yi],color=c,lw=6,alpha=0.5,solid_capstyle="round")
    axR.plot([lo,hi],[yi,yi],color=c,lw=1.5)
    axR.scatter([np.sqrt(lo*hi)],[yi],color=c,s=28,zorder=5)
    axR.text(140,yi,lab,fontsize=6.6,va="center",ha="left",color=c)
axR.axvspan(3,9,color="#c0392b",alpha=0.08,zorder=0)
axR.set_xscale("log"); axR.set_xlim(2,3000); axR.set_ylim(-0.6,len(rows)-0.4); axR.set_yticks([])
axR.set_xlabel("Linear sedimentation rate (m/Myr, log)")
axR.set_title("(b) Accumulation-rate context\nIquitos site is starved vs the Pebas basin average",fontsize=9.5,weight="bold")
axR.spines[["top","right","left"]].set_visible(False)
axR.text(6,len(rows)-0.5,"our margin",fontsize=6.5,color="#c0392b",ha="center")
fig.suptitle("Middle-Miocene accumulation: the Iquitos Pebas margin is condensed",fontsize=11,weight="bold",y=1.01)
fig.savefig(f"{OUT}/fig_sedrate_compare.png",dpi=170,bbox_inches="tight",facecolor="w")
print("wrote fig_sedrate_compare.png")
