#!/usr/bin/env python3
"""VGP diagnosis figure + local/regional/GPTS barcode comparison + La Venta U-Pb table.
Regional data from Mora-Rojas et al. 2023 (Geodiversitas 45(6)); GPTS = GTS2020/ATNTS2012.
"""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE="/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy"
OUT="/tmp/peru/outlaventa"; os.makedirs(OUT,exist_ok=True)

# ---------------------------------------------------------------- VGP figure
vdir=os.path.join(BASE,"OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed/dir/vgp_corrected")
vf=[f for f in glob.glob(os.path.join(vdir,"VGP*.txt"))]
frames=[]
for f in vf:
    df=pd.read_csv(f,sep="\t")
    df.columns=[c.strip() for c in df.columns]
    df["site"]="SC" if os.path.basename(f).split(" - ")[1].startswith("SC") else "TAM"
    frames.append(df)
V=pd.concat(frames,ignore_index=True)
vlat=[c for c in V.columns if "VGP_Lat" in c][0]; vlon=[c for c in V.columns if "VGP_Long" in c][0]

def vgp_lat(site_lat,inc,dec):
    slat,I,D=np.radians(site_lat),np.radians(inc),np.radians(dec)
    p=np.arctan2(2.0,np.tan(I))
    return np.degrees(np.arcsin(np.sin(slat)*np.cos(p)+np.cos(slat)*np.sin(p)*np.cos(D)))

fig,ax=plt.subplots(1,2,figsize=(11,5))
# (a) observed VGP scatter: latitude vs longitude
for s,c in [("SC","#2980b9"),("TAM","#c0392b")]:
    d=V[V.site==s]
    ax[0].scatter(d[vlon],d[vlat],s=9,alpha=.5,color=c,label=f"{s} (n={len(d)})")
ax[0].axhline(90,color="0.4",ls=":"); ax[0].axhline(-90,color="0.4",ls=":")
ax[0].axhspan(80,90,color="0.85"); ax[0].axhspan(-90,-80,color="0.85")
ax[0].text(-175,84,"where poles SHOULD sit (normal)",fontsize=7,color="0.3")
ax[0].set_xlabel("VGP longitude (deg)"); ax[0].set_ylabel("VGP latitude (deg)")
ax[0].set_ylim(-95,95); ax[0].set_xlim(-180,180); ax[0].legend(fontsize=8)
ax[0].set_title("Observed VGPs: scatter across ALL latitudes/longitudes\n(not clustered at the poles)")
# (b) why: declination sensitivity at shallow inclination
D=np.linspace(0,180,181)
ax[1].plot(D,[vgp_lat(-3.4,-8,d) for d in D],lw=2,color="#2980b9",label="shallow  I=-8 (~4S normal)")
ax[1].plot(D,[vgp_lat(-3.4,-60,d) for d in D],lw=2,color="#7f8c8d",label="steep  I=-60")
ax[1].axhline(90,color="0.5",ls=":"); ax[1].axhline(-90,color="0.5",ls=":")
ax[1].set_xlabel("assumed declination (deg) — UNKNOWN for unoriented cores")
ax[1].set_ylabel("computed VGP latitude (deg)"); ax[1].legend(fontsize=8)
ax[1].set_title("Why: at shallow I, VGP latitude swings +89->-83\nas the (arbitrary) declination varies")
fig.suptitle("VGP scatter is an artifact of unoriented cores at low paleolatitude — "
             "polarity must be read from inclination sign, not VGPs",y=1.03,fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_vgp_diagnosis.png"),dpi=200,bbox_inches="tight"); plt.close()
print("wrote fig_vgp_diagnosis.png  (VGP_lat median SC=%.0f TAM=%.0f)"
      %(V[V.site=='SC'][vlat].median(),V[V.site=='TAM'][vlat].median()))

# ---------------------------------------------------------------- La Venta U-Pb table
upb=pd.DataFrame([
 dict(sample="TVV-01",strat_m=932,median_Ma=11.46,CI95_low=10.52,CI95_high=12.09,type="detrital (MDA)"),
 dict(sample="SEL-02",strat_m=503,median_Ma=12.30,CI95_low=11.79,CI95_high=12.74,type="tuff (YSP 12.8+-0.1)"),
 dict(sample="TVV-04",strat_m=474,median_Ma=12.31,CI95_low=11.80,CI95_high=12.77,type="tuff (YSP 12.8+-0.1)"),
 dict(sample="44011", strat_m=218.8,median_Ma=13.12,CI95_low=12.69,CI95_high=13.47,type="tuff (YSP 13.9+-0.3)"),
 dict(sample="44017", strat_m=63.5,median_Ma=16.00,CI95_low=15.94,CI95_high=16.00,type="tuff (YSP 17.6+-0.4)"),
])
upb.to_csv(os.path.join(OUT,"laventa_upb_ages.csv"),index=False)

# ---------------------------------------------------------------- barcode comparison
NORM={"C5An.1n":(12.049,12.174),"C5An.2n":(12.272,12.474),"C5Ar.1n":(12.735,12.770),
 "C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),"C5ABn":(13.363,13.608),"C5ACn":(13.739,14.070)}
def draw_gpts(ax,lo,hi,x0,x1,label=True):
    ax.add_patch(plt.Rectangle((x0,lo),x1-x0,hi-lo,fc="white",ec="k",lw=.8))
    for nm,(t,b) in NORM.items():
        if b>lo and t<hi:
            ax.add_patch(plt.Rectangle((x0,max(t,lo)),x1-x0,min(b,hi)-max(t,lo),fc="k"))
            if label and lo<=(t+b)/2<=hi: ax.text(x1+0.01,(t+b)/2,nm,va="center",fontsize=6)

LO,HI=11.8,14.2
# Honda occupied interval C5ABn(base 13.608) -> C5An.1n(top 12.049)
HONDA_TOP,HONDA_BASE=12.049,13.608
MONKEY=(12.829,13.363)   # C5AAr..C5Ar.2n = Monkey Beds / P. neivensis interval
LAVENTAN=(11.8,13.5)
PALYNO_FLOOR=12.7

# Peru preferred correlation: SC reversal = C5AAr; place columns in age.
# SC zones (cm): N0-103, R106-149, N150-173,179-227,245-262. Tie R to C5AAr (13.183-13.363).
def sc_age(cm, r_top_cm=106, r_bot_cm=149, r_top_age=13.183, r_bot_age=13.363):
    sed=(r_bot_cm-r_top_cm)/(r_bot_age-r_top_age)  # cm/Ma
    return r_top_age + (cm-r_top_cm)/sed
SC=[(0,103,"N"),(106,149,"R"),(150,173,"N"),(179,227,"N"),(245,262,"N")]
TAMz=[(0,50,"R"),(51,72,"N"),(78,99,"N")]
# TAM: tie its R (0-50) to C5AAr too (top of TAM section = top of R ~ C5AAr top)
def tam_age(cm, r_top_cm=0, r_bot_cm=50, r_top_age=13.20, r_bot_age=13.363):
    sed=(r_bot_cm-r_top_cm)/(r_bot_age-r_top_age)
    return r_top_age+(cm-r_top_cm)/sed

fig,ax=plt.subplots(figsize=(10,8))
# Col1 GPTS
draw_gpts(ax,LO,HI,0.02,0.14)
ax.text(0.08,HI+0.03,"GPTS\n(GTS2020)",ha="center",va="bottom",fontsize=8)
# Col2 La Venta Honda: greyed full, occupied solid, Monkey band, U-Pb points
x0,x1=0.30,0.42
ax.add_patch(plt.Rectangle((x0,LO),x1-x0,HI-LO,fc="0.93",ec="k",lw=.8))
for nm,(t,b) in NORM.items():          # occupied chrons drawn black only within Honda interval
    tt,bb=max(t,HONDA_TOP),min(b,HONDA_BASE)
    if bb>tt: ax.add_patch(plt.Rectangle((x0,tt),x1-x0,bb-tt,fc="k"))
ax.add_patch(plt.Rectangle((x0,HONDA_TOP),x1-x0,HONDA_BASE-HONDA_TOP,fc="none",ec="#e67e22",lw=1.6))
ax.text(0.36,HI+0.03,"La Venta\n(Honda Gp)",ha="center",va="bottom",fontsize=8)
# Monkey/Purussaurus band
ax.add_patch(plt.Rectangle((x0-0.005,MONKEY[0]),(x1-x0)+0.01,MONKEY[1]-MONKEY[0],fc="none",ec="#16a085",lw=2,ls="--"))
ax.text(x1+0.09,(MONKEY[0]+MONKEY[1])/2,"Monkey Beds /\nP. neivensis\n(C5AAr–C5Ar.2n)",va="center",fontsize=6.5,color="#0e6655")
# U-Pb ages (median + CI) as points at a fixed x
xu=0.47
for _,r in upb.iterrows():
    if LO<=r.median_Ma<=HI:
        ax.plot([xu,xu],[r.CI95_low,r.CI95_high],color="#8e44ad",lw=1)
        ax.plot(xu,r.median_Ma,"o",color="#8e44ad",ms=4)
        ax.text(xu+0.015,r.median_Ma,r["sample"],fontsize=5.5,va="center",color="#6c3483")
ax.text(xu,HI+0.03,"U-Pb\nzircon",ha="center",va="bottom",fontsize=7,color="#6c3483")
# Col3 SC (Peru), Col4 TAM (Peru)
def draw_peru(ax,zones,agefn,x0,x1,title):
    for a,b,p in zones:
        aa,bb=agefn(a),agefn(b)
        ax.add_patch(plt.Rectangle((x0,min(aa,bb)),x1-x0,abs(bb-aa),fc=("k" if p=="N" else "white"),ec="k",lw=.6))
    ax.text((x0+x1)/2,HI+0.03,title,ha="center",va="bottom",fontsize=8)
draw_peru(ax,SC,sc_age,0.62,0.72,"SC\n(Peru)")
draw_peru(ax,TAMz,tam_age,0.80,0.90,"TAM\n(Peru)")
# constraint bands across all
ax.axhline(PALYNO_FLOOR,color="#c0392b",ls="--",lw=1.2)
ax.text(0.98,PALYNO_FLOOR,"palyno T15/T16\nfloor 12.7 Ma",fontsize=6,color="#c0392b",va="center",ha="right")
ax.set_ylim(HI,LO); ax.set_xlim(0,1.05); ax.set_xticks([]); ax.set_ylabel("Age (Ma)")
ax.set_title("Local (Peru TAM/SC) vs regional (La Venta Honda Gp, U-Pb-calibrated) vs GPTS\n"
             "Preferred tie: Peru reversal = C5AAr = La Venta Monkey-Beds / Purussaurus interval (~13.2–13.4 Ma)",
             fontsize=9.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_barcode_comparison.png"),dpi=200,bbox_inches="tight"); plt.close()
print("wrote fig_barcode_comparison.png, laventa_upb_ages.csv")
print(upb.to_string(index=False))
