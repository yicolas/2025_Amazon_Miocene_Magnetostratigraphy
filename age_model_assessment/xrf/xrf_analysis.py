#!/usr/bin/env python3
"""XRF vs NRM-intensity: test lithologic control of remanence intensity.
Registers valid XRF (Fe, Ti, Ca, K...) to NRM intensity per group by shared
within-group depth (XRF pos_mm/10 == paleomag pos_cm), then correlates.
Outputs (./outxrf): xrf_registered.csv, xrf_correlations.csv, fig_xrf_*.png
"""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE="/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy"
OUT=os.path.join("/tmp/peru","outxrf"); os.makedirs(OUT,exist_ok=True)
xrf=pd.read_csv(os.path.join(BASE,".notebook_append/xrf_tidy.csv"))

# map group -> its NRM .int file
INT={
 "TAM_1_2_3B":"TAM_1_2_3B_NRM_01192026.NRM.int",
 "TAM_3A_4_5E_5D_5C":"TAM_3A_4_5E_5D_5C_NRM_01192026.NRM.int",
 "TAM_5B_5A_6_7":"TAM_5B_5A_6_7_NRM_01202026.NRM.int",
 "SC_1C_1B_1A_2_3C":"SC_1C_1B_1A_2_3C_NRM_01212026.NRM.int",
 "SC_3B_3A_4D_4C_4B_4A":"SC_3B_3A_4D_4C_4B_4A_NRM_01202026.NRM.int",
 "SC_5_6_7C_7B_7A":"SC_5_6_7C_7B_7A_NRM_01212026.NRM.int",
 "SC_8":"SC_8_NRM_01212026.NRM.int"}
def load_int(fn):
    p=os.path.join(BASE,"OSU_uchannel_magnetometer_NRM_data/int",fn)
    df=pd.read_csv(p,sep="\t")
    df.columns=[c.strip() for c in df.columns]
    df=df.rename(columns={"Position":"pos_cm"})
    df["pos_cm"]=pd.to_numeric(df["pos_cm"],errors="coerce")
    nrm_col = "Int 0" if "Int 0" in df.columns else None
    if nrm_col is None:
        return None   # no NRM step (e.g. TAM Group 1 two-step file) -> skip
    df["NRM"]=pd.to_numeric(df[nrm_col],errors="coerce")
    df["I300"]=pd.to_numeric(df["Int 300"],errors="coerce") if "Int 300" in df.columns else np.nan
    return df[["pos_cm","NRM","I300"]].dropna(subset=["pos_cm"])

def pear(x,y):
    m=np.isfinite(x)&np.isfinite(y)
    if m.sum()<5: return np.nan,m.sum()
    return np.corrcoef(x[m],y[m])[0,1],int(m.sum())

rows=[]; reg_all=[]
GROUPS=list(INT.keys())
for g in GROUPS:
    x=xrf[(xrf.group==g)&(xrf.valid==1)].copy()
    x=x[(x.Fe>0)]
    if x.empty: continue
    it=load_int(INT[g])
    if it is None:
        print(f"  skip {g}: no NRM (Int 0) step in intensity file"); continue
    it=it[it.NRM>0]
    # interpolate XRF (as f(pos_cm)) onto NRM positions, within overlap
    x=x.sort_values("pos_cm")
    lo=max(x.pos_cm.min(),it.pos_cm.min()); hi=min(x.pos_cm.max(),it.pos_cm.max())
    itc=it[(it.pos_cm>=lo)&(it.pos_cm<=hi)].copy()
    def interp(col): return np.interp(itc.pos_cm, x.pos_cm, x[col])
    itc["Fe"]=interp("Fe"); itc["Ti"]=interp("Ti"); itc["Ca"]=interp("Ca")
    itc["K"]=interp("K"); itc["Mn"]=interp("Mn"); itc["cps"]=interp("cps")
    itc["FeTi"]=itc.Fe/itc.Ti.replace(0,np.nan)
    itc["group"]=g
    reg_all.append(itc)
    lNRM=np.log10(itc.NRM.values)
    r_fe,n=pear(lNRM,np.log10(itc.Fe.replace(0,np.nan).values))
    r_feti,_=pear(lNRM,np.log10(itc.FeTi.values))
    r_ti,_=pear(lNRM,np.log10(itc.Ti.replace(0,np.nan).values))
    r_ca,_=pear(lNRM,np.log10(itc.Ca.replace(0,np.nan).values))
    r_cps,_=pear(lNRM,np.log10(itc.cps.replace(0,np.nan).values))
    rows.append(dict(group=g,n=n,r_logFe=round(r_fe,2),r_logFeTi=round(r_feti,2),
                     r_logTi=round(r_ti,2),r_logCa=round(r_ca,2),r_logCPS=round(r_cps,2)))
reg=pd.concat(reg_all,ignore_index=True)
reg.to_csv(os.path.join(OUT,"xrf_registered.csv"),index=False)
corr=pd.DataFrame(rows)

# pooled (within-group z-scored to remove per-group offsets)
def z(a):
    a=np.asarray(a,float); return (a-np.nanmean(a))/np.nanstd(a)
pooled=[]
for g,gd in reg.groupby("group"):
    pooled.append(pd.DataFrame({"zlNRM":z(np.log10(gd.NRM)),
        "zlFe":z(np.log10(gd.Fe.replace(0,np.nan))),
        "zlFeTi":z(np.log10(gd.FeTi))}))
P=pd.concat(pooled,ignore_index=True)
r_pool_fe,n_fe=pear(P.zlNRM.values,P.zlFe.values)
r_pool_feti,_=pear(P.zlNRM.values,P.zlFeTi.values)
corr.loc[len(corr)]=dict(group="POOLED(z)",n=n_fe,r_logFe=round(r_pool_fe,2),
    r_logFeTi=round(r_pool_feti,2),r_logTi=np.nan,r_logCa=np.nan,r_logCPS=np.nan)
# --- diagnostics: Fe vs S (pyrite) and Fe vs Ca (carbonate), valid XRF only ---
diag=[]
for g,gd in xrf[xrf.valid==1].groupby("group"):
    fe=np.log10(gd.Fe.replace(0,np.nan)).values
    rS,_=pear(fe,np.log10(gd.S.replace(0,np.nan)).values)
    rCa,_=pear(fe,np.log10(gd.Ca.replace(0,np.nan)).values)
    diag.append(dict(group=g,r_logFe_logS=round(rS,2),r_logFe_logCa=round(rCa,2)))
diag=pd.DataFrame(diag)
corr=corr.merge(diag,on="group",how="left")
corr.to_csv(os.path.join(OUT,"xrf_correlations.csv"),index=False)
print("XRF correlations (Pearson r on log10):")
print(corr.to_string(index=False))

# --- pyrite/carbonate figure: Fe vs S and Fe vs Ca (pooled, colored by group) ---
figp,axp=plt.subplots(1,2,figsize=(11,5)); cm=plt.cm.tab10
av=xrf[xrf.valid==1]
for i,(g,gd) in enumerate(av.groupby("group")):
    axp[0].scatter(np.log10(gd.S.replace(0,np.nan)),np.log10(gd.Fe.replace(0,np.nan)),s=6,alpha=.4,color=cm(i),label=g)
    axp[1].scatter(np.log10(gd.Ca.replace(0,np.nan)),np.log10(gd.Fe.replace(0,np.nan)),s=6,alpha=.4,color=cm(i))
rS_all,_=pear(np.log10(av.Fe.replace(0,np.nan)).values,np.log10(av.S.replace(0,np.nan)).values)
rCa_all,_=pear(np.log10(av.Fe.replace(0,np.nan)).values,np.log10(av.Ca.replace(0,np.nan)).values)
axp[0].set_xlabel("log$_{10}$ S"); axp[0].set_ylabel("log$_{10}$ Fe"); axp[0].set_title(f"Fe vs S (pyrite): pooled r={rS_all:+.2f}"); axp[0].legend(fontsize=6,ncol=2)
axp[1].set_xlabel("log$_{10}$ Ca"); axp[1].set_ylabel("log$_{10}$ Fe"); axp[1].set_title(f"Fe vs Ca (carbonate): pooled r={rCa_all:+.2f}")
plt.suptitle("Most iron is authigenic (pyrite / carbonate-associated), not detrital magnetite\n"
             "-> XRF total Fe is a poor proxy for the remanence carrier",y=1.03,fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_xrf_pyrite.png"),dpi=200,bbox_inches="tight"); plt.close()

# ---- figure 1: pooled cross-plot log Fe vs log NRM
fig,ax=plt.subplots(1,2,figsize=(11,5))
cmap=plt.cm.tab10
for i,(g,gd) in enumerate(reg.groupby("group")):
    ax[0].scatter(np.log10(gd.Fe.replace(0,np.nan)),np.log10(gd.NRM),s=8,alpha=.5,color=cmap(i),label=g)
ax[0].set_xlabel("log$_{10}$ Fe (XRF peak area)"); ax[0].set_ylabel("log$_{10}$ NRM intensity")
ax[0].set_title(f"Fe vs NRM intensity (pooled z-scored r={r_pool_fe:+.2f}, n={n_fe})")
ax[0].legend(fontsize=6,ncol=2)
ax[1].scatter(P.zlFe,P.zlNRM,s=8,alpha=.4,color="0.3")
b=np.polyfit(P.zlFe[np.isfinite(P.zlFe)&np.isfinite(P.zlNRM)],P.zlNRM[np.isfinite(P.zlFe)&np.isfinite(P.zlNRM)],1)
xx=np.linspace(-3,3,10); ax[1].plot(xx,b[0]*xx+b[1],"r-",lw=1.5)
ax[1].set_xlabel("z(log Fe), within-group"); ax[1].set_ylabel("z(log NRM)")
ax[1].set_title(f"Within-group standardized: r={r_pool_fe:+.2f}")
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_xrf_Fe_vs_NRM.png"),dpi=200,bbox_inches="tight"); plt.close()

# ---- figure 2: downhole panels per group (NRM, Fe, Fe/Ti, Ca)
gs=[g for g in GROUPS if g in reg.group.values]
fig,axes=plt.subplots(1,len(gs),figsize=(2.4*len(gs),8),sharey=False)
for ax,g in zip(axes,gs):
    gd=reg[reg.group==g]
    ax.plot(z(np.log10(gd.NRM)),gd.pos_cm,color="k",lw=1,label="log NRM")
    ax.plot(z(np.log10(gd.Fe.replace(0,np.nan))),gd.pos_cm,color="#c0392b",lw=1,label="log Fe")
    ax.plot(z(np.log10(gd.FeTi)),gd.pos_cm,color="#2980b9",lw=0.9,label="log Fe/Ti")
    ax.invert_yaxis(); ax.set_title(g.replace("_","\n"),fontsize=7); ax.set_xlabel("z-score")
    if ax is axes[0]: ax.set_ylabel("within-group depth (cm)"); ax.legend(fontsize=6)
plt.suptitle("Downhole (z-scored): NRM intensity vs Fe and Fe/Ti — minima coincide = lithologic control",y=1.02,fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_xrf_downhole.png"),dpi=200,bbox_inches="tight"); plt.close()
print("\nwrote xrf_registered.csv, xrf_correlations.csv, fig_xrf_Fe_vs_NRM.png, fig_xrf_downhole.png")
