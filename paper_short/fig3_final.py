#!/usr/bin/env python3
"""FINAL Figure 3 (v2):
   TOP  : dense down-section columns (TAM + SC) — NRM, MDF, GRM, inc, dec, MAD, polarity
   BOTTOM row (3 panels):
     A  merged AF-demag decay curves (TAM + SC overlaid)
     B  GRM ratio vs MDF crossplot (greigite screen)
     C  representative Zijderveld (best clean specimen)
Exports flat PNG + layered SVG/PDF + per-panel PNGs."""
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd
import fig3_rockmag_dense as R      # dense columns: draw(), DAT, POL, OFF, LAN, BLOCK
import fig3_bottom_candidates as B  # SER, SCOL, mdf_grm
plt.close("all")
OUT=R.OUT; LAYERS=os.path.join(OUT,"fig3_layers"); os.makedirs(LAYERS,exist_ok=True)

# ---- choose the single best representative Zijderveld specimen ----
def best_specimen():
    cand=[]
    for site in ("TAM","SC"):
        for (g,pos),df in B.SER[site].items():
            if len(df)<6 or df.M.iloc[0]<=0: continue
            m,r=B.mdf_grm(df); nrm=df.M.iloc[0]; frac_end=df.M.iloc[-1]/nrm
            if not (np.isfinite(m) and np.isfinite(r)): continue
            if not (0.85<r<1.15): continue          # exclude greigite (GRM) specimens
            if m<15 or frac_end>0.45: continue        # want stable + well-demagnetized
            # straightness: fraction of variance along principal direction
            V=df[["X","Y","Z"]].values; V=V-V.mean(0)
            if len(V)<4: continue
            s=np.linalg.svd(V,compute_uv=False); lin=1-s[1]/s[0]
            cand.append((lin*np.log10(nrm*1e9+1), site, (g,pos), m, r, nrm))
    cand.sort(key=lambda t:-t[0]); return cand[0] if cand else None
BEST=best_specimen()

# ================= figure =================
fig=plt.figure(figsize=(15,11.8))
gs=GridSpec(2,4,height_ratios=[3.0,1.15],hspace=0.16,wspace=0.30,
            left=0.055,right=0.985,top=0.945,bottom=0.075)

# ---------- TOP: dense columns ----------
axtop=fig.add_subplot(gs[0,:]); axtop.set_gid("panel_columns")
R.draw(axtop,"TAM",0.0); R.draw(axtop,"SC",R.BLOCK)
allH=pd.concat([R.DAT['TAM'].H,R.DAT['SC'].H])
axtop.set_ylim(allH.min()-0.4,allH.max()+0.9); axtop.set_xlim(-0.5,R.BLOCK+R.LAN["POL"][1]+0.3)
axtop.set_ylabel("Height in reversal-aligned frame (m)",fontsize=10)
axtop.set_xticks([]); [axtop.spines[s].set_visible(False) for s in ("top","right","bottom")]
leg=[Line2D([0],[0],color="0.6",lw=1,label="raw"),Line2D([0],[0],color="k",lw=1.4,label="smoothed"),
     Line2D([0],[0],color="#1b7837",lw=1,label="MDF"),Line2D([0],[0],color="#c0392b",lw=1,label="GRM ratio (>1 greigite)"),
     Line2D([0],[0],color="#6a51a3",lw=1,label="MAD"),Line2D([0],[0],marker="o",color="w",mfc="#2166ac",ms=6,label="inc<0 (normal)"),
     Line2D([0],[0],marker="o",color="w",mfc="#b2182b",ms=6,label="inc>0 (reversed)"),
     Patch(fc="k",ec="k",label="N"),Patch(fc="w",ec="k",label="R"),Patch(fc="white",ec="0.6",hatch="////",label="no data")]
axtop.legend(handles=leg,loc="upper center",ncol=5,fontsize=7.5,frameon=False,bbox_to_anchor=(0.5,1.005))

# ---------- A: merged decay curves ----------
axA=fig.add_subplot(gs[1,0]); axA.set_gid("panel_A_decay")
grid=np.linspace(0,100,51)
for site in ("TAM","SC"):
    stk=[]
    for (g,pos),df in B.SER[site].items():
        if len(df)<4 or df.M.iloc[0]<=0: continue
        stk.append(np.interp(grid,df.AF,df.M.values/df.M.values[0],left=np.nan,right=np.nan))
    stk=np.array(stk); med=np.nanmedian(stk,0); q1=np.nanpercentile(stk,25,0); q3=np.nanpercentile(stk,75,0)
    axA.fill_between(grid,q1,q3,color=B.SCOL[site],alpha=0.18,lw=0)
    axA.plot(grid,med,color=B.SCOL[site],lw=2.4,label=f"{site} (n={len(stk)})")
axA.axhline(0.5,color="0.6",lw=0.6,ls=":"); axA.text(2,0.53,"MDF",fontsize=7,color="0.45")
axA.axvspan(70,100,color="#c0392b",alpha=0.06); axA.text(85,0.66,"GRM upturn\n(greigite)",fontsize=6.5,color="#c0392b",ha="center",va="center")
axA.set_xlim(0,100); axA.set_ylim(0,1.35); axA.set_xlabel("AF peak field (mT)",fontsize=9)
axA.set_ylabel("M / M$_0$",fontsize=9); axA.legend(fontsize=8,loc="upper right",frameon=False)
axA.set_title("A · AF-demag decay (median ± IQR)",fontsize=9,weight="bold")
axA.spines[["top","right"]].set_visible(False)

# ---------- B: GRM ratio vs MDF ----------
axB=fig.add_subplot(gs[1,1]); axB.set_gid("panel_B_grm")
axB.axhspan(1.0,3.0,color="#c0392b",alpha=0.07,zorder=0); axB.axhline(1.0,color="#c0392b",lw=0.8,ls="--")
axB.text(58,1.04,"greigite (GRM>1)",fontsize=7,color="#c0392b",ha="right",va="bottom")
for site in ("TAM","SC"):
    xs=[];ys=[]
    for (g,pos),df in B.SER[site].items():
        m,r=B.mdf_grm(df)
        if np.isfinite(m) and np.isfinite(r): xs.append(m);ys.append(r)
    axB.scatter(xs,ys,s=12,color=B.SCOL[site],alpha=0.5,lw=0,label=f"{site} (n={len(xs)})")
axB.set_xlim(0,60); axB.set_ylim(0,2.6); axB.set_xlabel("MDF (mT)",fontsize=9)
axB.set_ylabel("GRM ratio  I$_{100}$/I$_{70}$",fontsize=9); axB.legend(fontsize=8,frameon=False,loc="upper left")
axB.set_title("B · GRM screen (greigite)",fontsize=9,weight="bold"); axB.spines[["top","right"]].set_visible(False)

# ---------- C: representative Zijderveld ----------
axC=fig.add_subplot(gs[1,2]); axC.set_gid("panel_C_zij")
if BEST:
    _,site,(g,pos),mdf,grm,nrm=BEST; df=B.SER[site][(g,pos)]
    X,Y,Z=df.X.values*1e9,df.Y.values*1e9,df.Z.values*1e9
    axC.plot(Y,X,"-o",color="k",ms=4,lw=1.0,mfc="k",label="horizontal (N vs E)")
    axC.plot(Y,-Z,"-s",color="0.55",ms=4,lw=1.0,mfc="none",label="vertical (Up vs E)")
    axC.plot(Y[0],X[0],"o",color="#c0392b",ms=7,zorder=6)  # NRM start
    lim=np.nanmax(np.abs([X,Y,Z]))*1.15
    axC.axhline(0,color="0.7",lw=0.6); axC.axvline(0,color="0.7",lw=0.6)
    axC.set_xlim(-lim,lim); axC.set_ylim(-lim,lim); axC.set_aspect("equal")
    axC.set_xlabel("E (nA m²)",fontsize=8); axC.set_ylabel("N / Up (nA m²)",fontsize=8)
    axC.legend(fontsize=6.6,loc="lower right",frameon=False)
    axC.set_title(f"C · Zijderveld — {site} pos {pos:.0f}\nMDF {mdf:.0f} mT, GRM {grm:.2f} (ChRM to origin)",
                  fontsize=8,weight="bold")
    axC.tick_params(labelsize=6.5)

# ---------- D: SC inclination distribution -> paleolatitude (orientation-independent) ----------
axD=fig.add_subplot(gs[1,3]); axD.set_gid("panel_D_inc")
dS=R.S[R.S.group.isin(R.ORDER['SC'])].dropna(subset=["inc"])
incs=dS.inc.values.astype(float); pol=dS.polarity.values
site_lat=-3.36
I_gad=np.degrees(np.arctan(2*np.tan(np.radians(abs(site_lat)))))      # GAD-expected |I| (~6.7 deg)
binw=3.0; bins=np.arange(-42,42+binw,binw)
axD.hist(incs[pol=="N"],bins=bins,orientation="horizontal",color="#2166ac",alpha=0.80,label="normal (inc<0)")
axD.hist(incs[pol=="R"],bins=bins,orientation="horizontal",color="#b2182b",alpha=0.80,label="reversed (inc>0)")
axD.axhline(0,color="0.6",lw=0.6)
axD.axhline(-I_gad,color="k",ls="--",lw=0.8); axD.axhline(I_gad,color="k",ls="--",lw=0.8)
mean_absinc=np.nanmean(np.abs(incs))
paleolat=np.degrees(np.arctan(np.tan(np.radians(mean_absinc))/2.0))
axD.set_ylim(-42,42); axD.set_ylabel("ChRM inclination (deg)",fontsize=7.5)
axD.set_xlabel("count",fontsize=7.5); axD.tick_params(labelsize=6.5)
# right axis: paleolatitude equivalent of the (signed) inclination via tan I = 2 tan(lat)
def _i2l(i):
    i=np.asarray(i,float); return np.sign(i)*np.degrees(np.arctan(np.tan(np.radians(np.abs(i)))/2.0))
def _l2i(l):
    l=np.asarray(l,float); return np.sign(l)*np.degrees(np.arctan(2*np.tan(np.radians(np.abs(l)))))
try:
    axr=axD.secondary_yaxis('right', functions=(_i2l,_l2i))
    axr.set_ylabel("paleolatitude (deg)",fontsize=7.5); axr.tick_params(labelsize=6.5)
except Exception as e:
    print("secondary axis skipped:",e)
axD.set_title(f"D · SC inclination → paleolatitude\nmean |I|={mean_absinc:.0f}° → λ={paleolat:.0f}°  (GAD |I|={I_gad:.0f}° at {abs(site_lat):.1f}°S)",
              fontsize=7.3,weight="bold")
axD.legend(fontsize=5.8,frameon=False,loc="lower right")
axD.text(0.5,-0.24,"dashed = GAD-expected inclination; declination not used (orientation-independent)",
         transform=axD.transAxes,ha="center",fontsize=5.6,color="0.45")

fig.suptitle("Figure 3 · Paleomagnetic & rock-magnetic data — TAM & SC",fontsize=13,weight="bold",y=0.985)

# ---------- exports ----------
fig.savefig(os.path.join(OUT,"fig3_final.png"),dpi=190,facecolor="w")
fig.savefig(os.path.join(OUT,"fig3_final.svg"),transparent=True)
fig.savefig(os.path.join(OUT,"fig3_final.pdf"),transparent=True)
def save_layer(keep):
    for a in fig.axes: a.set_visible(a.get_gid()==keep)
    fig.patch.set_alpha(0); fig.savefig(os.path.join(LAYERS,f"{keep}.png"),transparent=True,dpi=190)
for a in fig.axes:
    if a.get_gid(): save_layer(a.get_gid())
for a in fig.axes: a.set_visible(True)
print("wrote fig3_final.png/.svg/.pdf ; best Zijderveld =",BEST[1:4] if BEST else None)
