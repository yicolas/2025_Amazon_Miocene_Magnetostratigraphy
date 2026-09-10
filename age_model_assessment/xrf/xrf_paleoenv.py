#!/usr/bin/env python3
"""XRF paleoenvironment: downcore redox/provenance proxies for the SC composite
(the long record). Fe, S (pyrite), Ca (carbonate/shell), Fe/Ti (authigenic Fe
enrichment), K/Ti (clay/detrital), Mn (redox). Overlays the SC polarity column.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BASE="/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy"
OUT="/tmp/peru/outxrf"; os.makedirs(OUT,exist_ok=True)
xrf=pd.read_csv(os.path.join(BASE,".notebook_append/xrf_tidy.csv"))

# SC groups in stratigraphic order + their measured lengths (cm) for compositing
SC_ORDER=[("SC_1C_1B_1A_2_3C",123.5),("SC_3B_3A_4D_4C_4B_4A",136.5),
          ("SC_5_6_7C_7B_7A",113.5),("SC_8",50)]
off=0.0; parts=[]
for g,L in SC_ORDER:
    d=xrf[(xrf.group==g)&(xrf.valid==1)].copy()
    if d.empty: off+=L; continue
    d["cdepth"]=(d.pos_cm-d.pos_cm.min())+off
    parts.append(d); off+=L
d=pd.concat(parts,ignore_index=True).sort_values("cdepth")
d["FeTi"]=d.Fe/d.Ti.replace(0,np.nan); d["KTi"]=d.K/d.Ti.replace(0,np.nan)

# SC polarity zones (from synthesis) on the same composite axis (approx)
ZONES=[(0,103,"N"),(106,149,"R"),(150,173,"N"),(179,227,"N"),(245,262,"N")]

def smooth(x,y,w=7):
    y=np.asarray(y,float); m=np.isfinite(y)
    return x[m],pd.Series(y[m]).rolling(w,center=True,min_periods=1).median().values

fig,axes=plt.subplots(1,6,figsize=(14,8),sharey=True)
panels=[("Fe","#c0392b","Fe (magnetic+authigenic)"),
        ("S","#d4ac0d","S (pyrite)"),
        ("Ca","#2980b9","Ca (carbonate/shell)"),
        ("FeTi","#7d3c98","Fe/Ti (authigenic Fe)"),
        ("KTi","#148f77","K/Ti (clay/detrital)"),
        ("Mn","#7f8c8d","Mn (redox)")]
for ax,(col,c,lab) in zip(axes,panels):
    xs,ys=smooth(d.cdepth.values,d[col].values)
    ax.plot(d[col],d.cdepth,".",ms=1.5,color=c,alpha=.25)
    ax.plot(ys,xs,"-",color=c,lw=1.3)
    ax.set_xlabel(lab,fontsize=8)
    if col in("Fe","S","Ca"): ax.set_xscale("log")
axes[0].set_ylim(d.cdepth.max()+10,-10)   # 0 (top/young) at top, deep/old at bottom
axes[0].set_ylabel("SC composite scan depth (cm)  — young (top) to old (bottom)")
# polarity column on the far right (add an axis)
axp=fig.add_axes([0.92,0.11,0.03,0.77])
for a,b,p in ZONES:
    axp.axhspan(a,b,color=("k" if p=="N" else "w"),ec="k")
    axp.text(0.5,(a+b)/2,p,ha="center",va="center",color=("w" if p=="N" else "k"),fontweight="bold",fontsize=8)
axp.set_ylim(axes[0].get_ylim()); axp.set_xlim(0,1); axp.set_xticks([]); axp.set_title("pol",fontsize=8)
fig.suptitle("SC elemental stratigraphy — redox/provenance proxies (Pebas mega-wetland). "
             "High Fe–S–Ca = pyritic, reducing, shell-bearing (marine-influenced) intervals.",y=0.98,fontsize=10)
plt.tight_layout(rect=[0,0,0.9,1])
plt.savefig(os.path.join(OUT,"fig_xrf_paleoenv.png"),dpi=200,bbox_inches="tight"); plt.close()

# quick stats for the text
print("SC valid XRF n=",len(d))
for col in ["Fe","S","Ca","FeTi","KTi","Mn"]:
    print(f"  {col:5} median={np.nanmedian(d[col]):.1f}  p10={np.nanpercentile(d[col],10):.1f}  p90={np.nanpercentile(d[col],90):.1f}")
# correlation Fe-S, Fe-Ca, S-Ca (pooled all valid)
av=xrf[xrf.valid==1]
def r(a,b):
    la,lb=np.log10(av[a].replace(0,np.nan)),np.log10(av[b].replace(0,np.nan)); m=np.isfinite(la)&np.isfinite(lb)
    return np.corrcoef(la[m],lb[m])[0,1]
print(f"\npooled log-log r: Fe-S={r('Fe','S'):+.2f} Fe-Ca={r('Fe','Ca'):+.2f} S-Ca={r('S','Ca'):+.2f} Ca-Sr={r('Ca','Sr'):+.2f}")
print("wrote fig_xrf_paleoenv.png")
