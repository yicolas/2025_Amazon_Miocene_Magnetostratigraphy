#!/usr/bin/env python3
"""
Peru Amazon Miocene — internal correlation + multi-GPTS age-model assessment
============================================================================
Consumes the polarity synthesis produced earlier (polarity_zones.csv) plus the
18-candidate cross_parameter_comparison.csv, and tests the observed polarity
against three timescale calibrations and the independent Amazon biostratigraphy.

Outputs (./out2):
  gpts_multiscale.csv        chron ages in GTS2020/GTS2012/CK95
  age_scenarios.csv          each candidate model -> implied ages + chron misfit
  fig_internal_corr.png      TAM vs SC polarity on sediment depth
  fig_multiscale_gpts.png    observed polarity vs 3 timescales + biostrat band
  fig_nonuniqueness.png      the single reversal tie is non-unique
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "out2"); os.makedirs(OUT, exist_ok=True)
ZONES = pd.read_csv(os.path.join(HERE, "out", "polarity_zones.csv"))
XP    = pd.read_csv("/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy/"
                    "cross_parameter_comparison.csv")

# ----------------------------------------------------------------- GPTS tables
# normal subchrons; (top, base) Ma. GTS2020==GTS2012 in this window (both ATNTS2012).
GPTS = {
 "GTS2020": {"C5n.2n":(9.984,11.056),"C5An.1n":(12.049,12.174),"C5An.2n":(12.272,12.474),
   "C5Ar.1n":(12.735,12.770),"C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),
   "C5ABn":(13.363,13.608),"C5ACn":(13.739,14.070),"C5ADn":(14.163,14.609)},
 "GTS2012": {"C5n.2n":(9.984,11.056),"C5An.1n":(12.049,12.174),"C5An.2n":(12.272,12.474),
   "C5Ar.1n":(12.735,12.770),"C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),
   "C5ABn":(13.363,13.608),"C5ACn":(13.739,14.070),"C5ADn":(14.163,14.609)},
 "CK95":    {"C5n.2n":(9.920,10.949),"C5An.1n":(11.935,12.078),"C5An.2n":(12.184,12.401),
   "C5Ar.1n":(12.678,12.708),"C5Ar.2n":(12.775,12.819),"C5AAn":(12.991,13.139),
   "C5ABn":(13.302,13.510),"C5ACn":(13.703,14.076),"C5ADn":(14.178,14.612)},
}
pd.DataFrame([{"chron":c,
   "GTS2020_top":GPTS["GTS2020"][c][0],"GTS2020_base":GPTS["GTS2020"][c][1],
   "GTS2012_top":GPTS["GTS2012"][c][0],"GTS2012_base":GPTS["GTS2012"][c][1],
   "CK95_top":GPTS["CK95"][c][0],"CK95_base":GPTS["CK95"][c][1]}
   for c in GPTS["GTS2020"]]).to_csv(os.path.join(OUT,"gpts_multiscale.csv"),index=False)

# reversed intervals (gaps between normals) for each scale, in the study window
def rev_intervals(scale, lo=11.8, hi=14.8):
    norms = sorted(GPTS[scale].values())
    revs=[]
    for (t0,b0),(t1,b1) in zip(norms[:-1],norms[1:]):
        if b0<hi and t1>lo: revs.append((b0,t1))   # base of younger-normal .. top of older-normal
    return revs

# biostrat constraint (from research): T15/T16 boundary 12.7 Ma; window T-15(14.2-12.7)+base T-16;
# marine incursion ME-4 = 14.6-11.62 Ma. Pebas has ~no datable ash.
BIOSTRAT = {"T15_T16_boundary":12.7, "T15":(12.7,14.2), "T16_base":(7.1,12.7), "ME4":(11.62,14.6)}

# ----------------------------------------------------------------- age scenarios
# depth->age for a linear model given (top_age, sed cm/Ma) and a site's zone depths.
def zone_ages(site, top_age, sed):
    z = ZONES[ZONES.site==site].copy()
    z0 = z.top_cm.min()
    z["top_age"]  = top_age + (z.top_cm   - z0)/sed
    z["bot_age"]  = top_age + (z.bottom_cm- z0)/sed
    return z

def nearest_chron_boundary(age, scale):
    bnds=[]
    for c,(t,b) in GPTS[scale].items():
        bnds.append((t,f"{c}.top")); bnds.append((b,f"{c}.base"))
    a=min(bnds,key=lambda x:abs(x[0]-age))
    return a[1], round(a[0]-age,3)

# candidate models = the inclination + declination rows from cross_parameter_comparison
rows=[]
for _,r in XP.iterrows():
    site=r.Site; top=r.Start_Ma; sed=r.Sed_Rate_cm_Ma
    z=zone_ages(site, top, sed)
    # for each site's reversed zones, find misfit to nearest boundary of each scale
    revs=z[z.polarity=="R"]
    for _,zz in revs.iterrows():
        for scale in ("GTS2020","CK95"):
            bt,dt=nearest_chron_boundary(zz.top_age,scale)
            bb,db=nearest_chron_boundary(zz.bot_age,scale)
            rows.append(dict(param=r.Model,rank=r.Rank,site=site,score=r.Score_Value,
                corr=r.Correlation,top_Ma=r.Start_Ma,end_Ma=r.End_Ma,sed=sed,scale=scale,
                Rzone_top_cm=zz.top_cm,Rzone_bot_cm=zz.bottom_cm,
                Rzone_top_age=round(zz.top_age,3),Rzone_bot_age=round(zz.bot_age,3),
                near_top=bt,dtop_kyr=int(round(dt*1000)),near_bot=bb,dbot_kyr=int(round(db*1000)),
                in_ME4= (BIOSTRAT["ME4"][0]<=zz.top_age<=BIOSTRAT["ME4"][1])))
scen=pd.DataFrame(rows)
scen.to_csv(os.path.join(OUT,"age_scenarios.csv"),index=False)

# ----------------------------------------------------------------- console
print("="*74)
print("INTERNAL POLARITY (from polarity_zones.csv)")
print("="*74)
for s in ("TAM","SC"):
    z=ZONES[ZONES.site==s]
    print(f"{s}: "+" ".join(f"{p}({int(a)}-{int(b)}cm)" for a,b,p in
          zip(z.top_cm,z.bottom_cm,z.polarity)))
print("\nNumber of Reversed zones: TAM =",(ZONES[(ZONES.site=='TAM')&(ZONES.polarity=='R')]).shape[0],
      " SC =",(ZONES[(ZONES.site=='SC')&(ZONES.polarity=='R')]).shape[0])

print("\n"+"="*74)
print("BEST INCLINATION MODELS (trusted param) — reversed-zone tie & misfit")
print("="*74)
inc=scen[(scen.param=="Inclination")]
for scale in ("GTS2020","CK95"):
    print(f"\n-- {scale} --")
    sub=inc[inc.scale==scale].sort_values(["site","rank"])
    for _,r in sub.iterrows():
        print(f"  {r['site']} rk{int(r['rank'])} sed={r['sed']:>4.0f} top={r['top_Ma']:.3f}  "
              f"R@{int(r['Rzone_top_cm']):>3}-{int(r['Rzone_bot_cm']):<3}cm -> "
              f"{r['Rzone_top_age']:.3f}-{r['Rzone_bot_age']:.3f}Ma  "
              f"near {r['near_top']}({r['dtop_kyr']:+d}kyr)/{r['near_bot']}({r['dbot_kyr']:+d}kyr)  "
              f"ME4={'Y' if r['in_ME4'] else 'n'}")

print("\n"+"="*74)
print("SED-RATE SPREAD across the 18 candidate models (cm/Ma)")
print("="*74)
for (site,param),g in XP.groupby(["Site","Model"]):
    print(f"  {site:3} {param:12} sed {g.Sed_Rate_cm_Ma.min():.0f}-{g.Sed_Rate_cm_Ma.max():.0f}  "
          f"start {g.Start_Ma.min():.3f}-{g.Start_Ma.max():.3f}  score {g.Score_Value.max():.3f}")

# ----------------------------------------------------------------- figures
def fig_internal():
    fig,axes=plt.subplots(1,2,figsize=(6,9),sharey=False)
    for ax,s in zip(axes,("TAM","SC")):
        z=ZONES[ZONES.site==s]
        for _,r in z.iterrows():
            ax.axhspan(r.top_cm,r.bottom_cm,color=("k" if r.polarity=="N" else "w"),ec="k")
            ax.text(0.5,(r.top_cm+r.bottom_cm)/2,r.polarity,ha="center",va="center",
                    color=("w" if r.polarity=="N" else "k"),fontweight="bold")
        ax.set_xlim(0,1);ax.set_xticks([]);ax.invert_yaxis()
        ax.set_title(f"Site {s}\n{z.shape[0]} zones, "
                     f"{int(z.bottom_cm.max())} cm");ax.set_ylabel("sediment depth (cm)")
    # tie line: the single well-constrained reversal in each
    fig.suptitle("Internal correlation — one shared, well-constrained reversal\n"
                 "(TAM R = TAM-4;  SC R = SC-4C/4B). Flanking normals correlate.",y=1.02,fontsize=10)
    plt.tight_layout();p=os.path.join(OUT,"fig_internal_corr.png")
    fig.savefig(p,dpi=200,bbox_inches="tight");plt.close();return p

def fig_multiscale():
    # observed polarity placed with the two best coeval inclination models
    models={"TAM":(12.935,778.0),"SC":(12.945,786.0)}
    fig,axes=plt.subplots(1,5,figsize=(13,8),sharey=True)
    lo,hi=11.9,14.7
    scales=["GTS2020","GTS2012","CK95"]
    for i,scale in enumerate(scales):
        ax=axes[i];ax.set_title(scale,fontsize=10);ax.set_ylim(hi,lo);ax.set_xlim(0,1);ax.set_xticks([])
        ax.axhspan(lo,hi,color="w",ec="k")
        for c,(t,b) in GPTS[scale].items():
            if b>lo and t<hi:
                ax.axhspan(t,b,color="k")
                if lo<=(t+b)/2<=hi: ax.text(1.05,(t+b)/2,c,va="center",fontsize=6)
        if i==0: ax.set_ylabel("Age (Ma)")
    for j,site in enumerate(("TAM","SC")):
        ax=axes[3+j];top,sed=models[site];z=zone_ages(site,top,sed)
        for _,r in z.iterrows():
            ax.axhspan(min(r.top_age,r.bot_age),max(r.top_age,r.bot_age),
                       color=("k" if r.polarity=="N" else "w"),ec="k")
        ax.set_title(f"{site}\n{top}Ma,{sed:.0f}cm/Ma",fontsize=9)
        ax.set_xlim(0,1);ax.set_xticks([])
    # biostrat band: T15/T16 boundary + ME4 window on every panel
    for ax in axes:
        ax.axhline(BIOSTRAT["T15_T16_boundary"],color="#c0392b",lw=1.3,ls="--")
        ax.axhspan(BIOSTRAT["ME4"][0],BIOSTRAT["ME4"][1],color="#3498db",alpha=0.06)
    axes[0].text(0.02,BIOSTRAT["T15_T16_boundary"]-0.03,"T15/T16 (12.7 Ma)",
                 color="#c0392b",fontsize=6)
    fig.suptitle("Observed polarity vs three timescales (black=Normal). Red dash = palyno "
                 "T15/T16 boundary; blue = ME-4 marine-incursion window.",y=1.02,fontsize=10)
    plt.tight_layout();p=os.path.join(OUT,"fig_multiscale_gpts.png")
    fig.savefig(p,dpi=200,bbox_inches="tight");plt.close();return p

def fig_nonuniqueness():
    # A single R-flanked-by-N pattern can tie to several reversed chrons. Show candidates.
    revs=rev_intervals("GTS2020")
    names=[]
    norms=sorted(GPTS["GTS2020"].items(),key=lambda kv:kv[1])
    for i in range(len(norms)-1):
        names.append(f"{norms[i][0]}r")   # reversed interval below each normal
    fig,ax=plt.subplots(figsize=(7,5))
    for i,((b0,t1),nm) in enumerate(zip(revs,names)):
        ax.barh(i,t1-b0,left=b0,color="0.8",ec="k")
        ax.text((b0+t1)/2,i,nm,ha="center",va="center",fontsize=8)
    ax.set_yticks(range(len(revs)));ax.set_yticklabels([f"cand {i+1}" for i in range(len(revs))])
    ax.set_xlabel("Age (Ma)");ax.invert_xaxis()
    ax.axvspan(BIOSTRAT["T15"][0],BIOSTRAT["T15"][1],color="#2ecc71",alpha=0.12,label="palyno T-15 window")
    ax.axvline(BIOSTRAT["T15_T16_boundary"],color="#c0392b",ls="--",label="T15/T16 (12.7 Ma)")
    ax.set_title("Non-uniqueness: one reversed zone flanked by normals could be any of these\n"
                 "mid-Miocene reversed chrons. Biostrat narrows, doesn't uniquely fix.",fontsize=9)
    ax.legend(fontsize=7,loc="lower left")
    plt.tight_layout();p=os.path.join(OUT,"fig_nonuniqueness.png")
    fig.savefig(p,dpi=200,bbox_inches="tight");plt.close();return p

for p in (fig_internal(),fig_multiscale(),fig_nonuniqueness()):
    print("wrote",os.path.basename(p))
print("wrote gpts_multiscale.csv, age_scenarios.csv")
print("DONE")
