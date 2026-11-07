#!/usr/bin/env python3
"""XRF raw-vs-cleaned diagnostic for Fe & S (SC and TAM), plotted vs field height.
Cleaning steps (conservative): (1) drop low-total-count points (cps<CPS_MIN) and
poor spectral fits (MSE>MSE_MAX); (2) break the line at inter-segment gaps (>GAP cm,
no bridging); (3) 3-point median despike within each contiguous run."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE="/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy"
PS=os.path.join(BASE,"polarity_synthesis"); OUT="/tmp/peru/fig_build"
CPS_MIN=80000; MSE_MAX=3.5; GAP=2.0            # <-- tunable cleaning thresholds
C_FE="#2a9d8f"; C_S="#e76f51"

PZ=pd.read_csv(os.path.join(PS,"polarity_zones.csv")); POL={s:d for s,d in PZ.groupby("site")}
S=pd.read_csv(os.path.join(PS,"samples_corrected.csv"))
ORDER={"TAM":["TAM_3A_4_5E_5D_5C","TAM_5B_5A_6_7"],
       "SC":["SC_1C_1B_1A_2_3C","SC_3B_3A_4D_4C_4B_4A","SC_5_6_7C_7B_7A","SC_8"]}
MARKpmag={"TAM":(0.00,2.54),"SC":(1.60,5.00)}; XTRA={"TAM":{"TAM_1_2_3B":(2.54,3.78)}}
def META(site):
    d=S[S.group.isin(ORDER[site])]; base=0.0; m={}
    for g in ORDER[site]:
        th=d[d.group==g].sed_depth_cm.max(); m[g]=(base,th); base+=th
    return m
def toH(site,cm):
    z=POL[site]; lo=z.top_cm.min(); hi=z.bottom_cm.max(); a,b=MARKpmag[site]
    return b-(cm-lo)/(hi-lo)*(b-a)

X=pd.read_csv(os.path.join(BASE,".notebook_append/xrf_tidy.csv"))
def register(site):
    m=META(site); out=[]
    Xv=X[(X.group.isin(ORDER[site]))].copy()
    for g in ORDER[site]:
        dg=Xv[Xv.group==g].copy()
        if not len(dg): continue
        b,th=m[g]; gmax=dg[dg.valid==1].pos_cm.max()
        dg["H"]=toH(site,b+(dg.pos_cm/gmax)*th)
        out.append(dg)
    for g,(fb,ft) in XTRA.get(site,{}).items():
        dg=X[X.group==g].copy()
        if not len(dg): continue
        gmax=dg[dg.valid==1].pos_cm.max(); dg["H"]=ft-(dg.pos_cm/gmax)*(ft-fb)
        out.append(dg)
    return pd.concat(out)

def clean(df, elem):
    """Return (raw_df, cleaned series with NaN breaks, removed_mask) sorted by H desc."""
    d=df[df.valid==1].copy().sort_values(["group","pos_cm"])
    d["removed"]=(d.cps<CPS_MIN)|(d.MSE>MSE_MAX)
    # build cleaned with gap breaks + median despike, per group
    rows=[]
    for g,dg in d.groupby("group"):
        dg=dg.sort_values("pos_cm").reset_index(drop=True)
        val=dg[elem].where(~dg.removed).astype(float)
        # despike: 3-pt rolling median on kept values
        val=val.rolling(3,center=True,min_periods=1).median()
        prev=None
        for i,r in dg.iterrows():
            if prev is not None and (r.pos_cm-prev)>GAP:
                rows.append(dict(H=np.nan, v=np.nan))     # break the line across a gap
            rows.append(dict(H=r.H, v=(np.nan if r.removed else val.iloc[i])))
            prev=r.pos_cm
    cl=pd.DataFrame(rows)
    return d, cl

fig,axes=plt.subplots(2,4,figsize=(15,9),sharey="row")
for row,site in enumerate(["SC","TAM"]):
    reg=register(site)
    for col,(elem,color) in enumerate([("Fe",C_FE),("S",C_S)]):
        araw=axes[row,2*col]; acln=axes[row,2*col+1]
        d,cl=clean(reg,elem)
        # RAW panel: all valid points, line bridges everything (current behaviour)
        dr=d.sort_values("H",ascending=False)
        araw.plot(dr[elem],dr.H,color="0.5",lw=0.5)
        rem=d[d.removed]
        araw.scatter(rem[elem],rem.H,s=10,c="red",marker="x",lw=0.7,zorder=5,label=f"flagged n={len(rem)}")
        araw.set_title(f"{site} {elem} — RAW (valid flag only)",fontsize=9)
        araw.legend(fontsize=6,loc="lower right")
        # CLEANED panel
        acln.plot(cl.v,cl.H,color=color,lw=0.7)
        acln.set_title(f"{site} {elem} — CLEANED",fontsize=9)
        for a in (araw,acln):
            a.set_xlabel(f"{elem} (cts)",fontsize=7); a.tick_params(labelsize=6); a.grid(alpha=.2)
    axes[row,0].set_ylabel(f"{site}\nfield height (m)",fontsize=9)
fig.suptitle(f"XRF Fe & S: raw vs cleaned  (drop cps<{CPS_MIN:,} or MSE>{MSE_MAX}; break gaps>{GAP} cm; 3-pt median despike)",
             fontsize=11,weight="bold")
plt.tight_layout(rect=[0,0,1,0.97])
fig.savefig(os.path.join(OUT,"xrf_diagnostic.png"),dpi=150,bbox_inches="tight",facecolor="w")
# report
for site in ["SC","TAM"]:
    reg=register(site)
    for elem in ["Fe","S"]:
        d,_=clean(reg,elem)
        print(f"{site} {elem}: valid={len(d)} flagged(cps/MSE)={d.removed.sum()}")
print("wrote xrf_diagnostic.png")
