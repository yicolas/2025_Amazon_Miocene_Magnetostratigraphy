import os
#!/usr/bin/env python3
"""Candidate BOTTOM panels for the 3-panel Fig 3. Renders four options as
separate PNGs so we can pick which to keep:
  A  af_demag_decay.png   normalized AF-demag decay curves (M/M0 vs field), per site
  B  grm_mdf_crossplot.png GRM ratio (I100/I70) vs MDF, greigite field shaded
  C  stereonets.png        equal-area ChRM directions, TAM & SC, Fisher mean + a95
  D  zijderveld.png        orthogonal vector plots for representative specimens
Derived from the u-channel demag DAT series + samples_corrected."""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

BASE=os.environ.get("PEBAS_BASE","/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy")
PS=os.path.join(BASE,"polarity_synthesis")
DATD=os.path.join(BASE,"OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed")
OUT=os.environ.get("PEBAS_OUT","/mnt/user-data/outputs")
S=pd.read_csv(os.path.join(PS,"samples_corrected.csv"))
ORDER={"TAM":["TAM_3A_4_5E_5D_5C","TAM_5B_5A_6_7"],
       "SC":["SC_1C_1B_1A_2_3C","SC_3B_3A_4D_4C_4B_4A","SC_5_6_7C_7B_7A","SC_8"]}
C_TAM="#1f6f8b"; C_SC="#b5651d"; SCOL={"TAM":C_TAM,"SC":C_SC}

def datfile(group):
    toks=group.split("_"); pref=[toks[0]]+(["SE"] if toks[0]=="SC" else [])
    for c in glob.glob(os.path.join(DATD,"*.DAT")):
        b=os.path.basename(c)
        if any(p in b for p in pref) and toks[1] in b: return c
    return None

def demag(group):
    """per-position demag series: dict pos -> DataFrame(AF,M,X,Y,Z,Inc,Dec)."""
    f=datfile(group)
    if not f: return {}
    d=pd.read_csv(f,sep="\t",header=0)
    d=d[d["Data type"].astype(str).str.contains("SAMPLE",na=False)].copy()
    d["AF"]=d[["AF X","AF Y","AF Z"]].max(axis=1)/10.0
    out={}
    for pos,dp in d.groupby("Position"):
        dp=dp.sort_values("AF")
        out[(group,pos)]=pd.DataFrame(dict(AF=dp["AF"].values,M=dp["Intensity"].values.astype(float),
            X=dp["X intensity"].values.astype(float),Y=dp["Y intensity"].values.astype(float),
            Z=dp["Z intensity"].values.astype(float),Inc=dp["Inclination: Unrotated"].values.astype(float),
            Dec=dp["Declination: Unrotated"].values.astype(float)))
    return out
SER={s:{k:v for g in ORDER[s] for k,v in demag(g).items()} for s in ("TAM","SC")}

def mdf_grm(df):
    af=df.AF.values; M=df.M.values
    if len(af)<3 or M[0]<=0: return np.nan,np.nan
    frac=M/M[0]; mdf=np.nan
    for k in range(1,len(af)):
        if frac[k]<=0.5: mdf=np.interp(0.5,[frac[k],frac[k-1]],[af[k],af[k-1]]); break
    at=lambda x:M[np.argmin(np.abs(af-x))]
    grm=at(100)/at(70) if at(70)>0 else np.nan
    return mdf,grm

# ======================= A: AF-demag decay curves =======================
figA,axes=plt.subplots(1,2,figsize=(11,4.6))
for ax,site in zip(axes,("TAM","SC")):
    grid=np.linspace(0,100,51); stack=[]
    for (g,pos),df in SER[site].items():
        if len(df)<4 or df.M.iloc[0]<=0: continue
        frac=df.M.values/df.M.values[0]
        ax.plot(df.AF,frac,color=SCOL[site],lw=0.4,alpha=0.10)
        stack.append(np.interp(grid,df.AF,frac,left=np.nan,right=np.nan))
    stk=np.array(stack); med=np.nanmedian(stk,0); q1=np.nanpercentile(stk,25,0); q3=np.nanpercentile(stk,75,0)
    ax.fill_betweenx(med,grid*0,grid*0,alpha=0)  # noop keep API
    ax.fill_between(grid,q1,q3,color=SCOL[site],alpha=0.20,lw=0)
    ax.plot(grid,med,color=SCOL[site],lw=2.2)
    ax.axhline(0.5,color="0.6",lw=0.6,ls=":"); ax.text(2,0.52,"MDF",fontsize=7,color="0.4")
    ax.set_title(f"{site}   (n={len(SER[site])})",fontsize=10,weight="bold",color=SCOL[site])
    ax.set_xlabel("AF peak field (mT)"); ax.set_xlim(0,100); ax.set_ylim(0,1.35)
    ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("M / M$_0$")
figA.suptitle("A · Normalized AF-demagnetization decay curves",fontsize=12,weight="bold")
figA.text(0.5,0.005,"Upturn above ~70 mT = gyroremanent (GRM) growth, the greigite fingerprint.",
          ha="center",fontsize=8,color="0.35")
figA.tight_layout(rect=[0,0.03,1,0.96]); figA.savefig(f"{OUT}/af_demag_decay.png",dpi=180,facecolor="w"); plt.close(figA)

# ======================= B: GRM-MDF crossplot =======================
figB,ax=plt.subplots(figsize=(6.2,5.4))
ax.axhspan(1.0,3.0,color="#c0392b",alpha=0.07,zorder=0)
ax.axhline(1.0,color="#c0392b",lw=0.8,ls="--")
ax.text(2,1.02,"greigite field  (GRM ratio > 1)",fontsize=8,color="#c0392b",va="bottom")
for site in ("TAM","SC"):
    xs=[];ys=[]
    for (g,pos),df in SER[site].items():
        m,r=mdf_grm(df)
        if np.isfinite(m) and np.isfinite(r): xs.append(m);ys.append(r)
    ax.scatter(xs,ys,s=14,color=SCOL[site],alpha=0.55,lw=0,label=f"{site} (n={len(xs)})")
ax.set_xlabel("MDF (mT)"); ax.set_ylabel("GRM ratio  I$_{100}$/I$_{70}$")
ax.set_xlim(0,60); ax.set_ylim(0,2.6); ax.legend(frameon=False,fontsize=9)
ax.spines[["top","right"]].set_visible(False)
figB.suptitle("Figure S1 · GRM ratio vs MDF (greigite screen)",fontsize=12,weight="bold")
figB.tight_layout(); figB.savefig(f"{OUT}/grm_mdf_crossplot.png",dpi=180,facecolor="w"); plt.close(figB)

# ======================= C: equal-area stereonets =======================
def eqarea(D,I):
    p=np.abs(I); r=np.sqrt(2)*np.sin(np.radians((90-p)/2))
    return r*np.sin(np.radians(D)), r*np.cos(np.radians(D))
def fisher(dec,inc):
    d=np.radians(dec); i=np.radians(inc)
    x=np.cos(i)*np.cos(d); y=np.cos(i)*np.sin(d); z=np.sin(i)
    R=np.sqrt(x.sum()**2+y.sum()**2+z.sum()**2); N=len(dec)
    xm,ym,zm=x.sum()/R,y.sum()/R,z.sum()/R
    Dm=np.degrees(np.arctan2(ym,xm))%360; Im=np.degrees(np.arcsin(zm))
    k=(N-1)/(N-R) if N>R else np.nan
    a95=np.degrees(np.arccos(1-(N-R)/R*(20**(1/(N-1))-1))) if N>R else np.nan
    return Dm,Im,k,a95,N
figC,axes=plt.subplots(1,2,figsize=(11,5.6))
for ax,site in zip(axes,("TAM","SC")):
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Circle((0,0),1,fc="none",ec="k",lw=1.2))
    for a in range(0,360,10):
        r0=0.98 if a%90 else 0.95
        ax.plot([r0*np.sin(np.radians(a)),np.sin(np.radians(a))],[r0*np.cos(np.radians(a)),np.cos(np.radians(a))],color="k",lw=0.6)
    for lab,a in [("N",0),("E",90),("S",180),("W",270)]:
        ax.text(1.08*np.sin(np.radians(a)),1.08*np.cos(np.radians(a)),lab,ha="center",va="center",fontsize=9)
    d=S[S.group.isin(ORDER[site])].dropna(subset=["inc","dec"])
    for pol,col,fill in [("N","#2166ac",True),("R","#b2182b",False)]:
        sub=d[d.polarity==pol]
        if not len(sub): continue
        xx,yy=eqarea(sub.dec.values,sub.inc.values)
        ax.scatter(xx,yy,s=16,facecolors=(col if fill else "none"),edgecolors=col,lw=0.8,alpha=0.8)
        Dm,Im,k,a95,N=fisher(sub.dec.values,sub.inc.values)
        mx,my=eqarea(np.array([Dm]),np.array([Im]))
        ax.scatter(mx,my,s=120,marker="*",facecolors=col,edgecolors="k",lw=0.8,zorder=6)
        # a95 small circle
        th=np.linspace(0,360,120); cx=[];cy=[]
        for t in th:
            # rotate a95 cone around mean (approx in projection)
            pass
        ax.text(-1.15,-1.02-(0.0 if pol=="N" else 0.14),
                f"{'N' if pol=='N' else 'R'}: D={Dm:.0f}° I={Im:.0f}° α95={a95:.1f}° (n={N})",
                fontsize=7.5,color=col,ha="left")
    ax.set_title(f"{site}",fontsize=11,weight="bold",color=SCOL[site])
    ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.3,1.2)
figC.suptitle("Figure S2 · Equal-area ChRM directions (filled = lower hemi / normal, open = upper / reversed)",fontsize=11,weight="bold")
figC.text(0.5,0.02,"Note: u-channel samples are not azimuthally oriented, so declination is RELATIVE; polarity is assigned on inclination sign. "
          "Clusters are therefore not expected to be antipodal in declination.",ha="center",fontsize=7.5,color="0.4")
figC.tight_layout(rect=[0,0.04,1,0.95]); figC.savefig(f"{OUT}/stereonets.png",dpi=180,facecolor="w"); plt.close(figC)

# ======================= D: Zijderveld =======================
def pick(site,n=3):
    """choose representative positions: strongest normal, strongest reversed, one greigite."""
    items=list(SER[site].items()); scored=[]
    for k,df in items:
        if len(df)<5 or df.M.iloc[0]<=0: continue
        m,r=mdf_grm(df); scored.append((k,df.Inc.iloc[0],r,df.M.iloc[0]))
    if not scored: return []
    norm=sorted([s for s in scored if s[1]<0],key=lambda z:-z[3])
    rev=sorted([s for s in scored if s[1]>0],key=lambda z:-z[3])
    gre=sorted([s for s in scored if np.isfinite(s[2])],key=lambda z:-z[2])
    out=[]
    if norm: out.append((norm[0][0],"normal"))
    if rev: out.append((rev[0][0],"reversed"))
    if gre and gre[0][0] not in [o[0] for o in out]: out.append((gre[0][0],"greigite GRM"))
    return out[:n]
figD,axes=plt.subplots(2,3,figsize=(12,7.6))
for row,site in enumerate(("TAM","SC")):
    reps=pick(site)
    for col in range(3):
        ax=axes[row][col]
        if col>=len(reps): ax.axis("off"); continue
        (g,pos),tag=reps[col]; df=SER[site][(g,pos)]
        X,Y,Z=df.X.values*1e9,df.Y.values*1e9,df.Z.values*1e9  # nAm^2 scale
        # horizontal (solid): E on x, N on y ; vertical (open): E on x, Up(-Z) on y
        ax.plot(Y,X,"-o",color="k",ms=3,lw=0.9,mfc="k",label="horiz (N vs E)")
        ax.plot(Y,-Z,"-s",color="0.5",ms=3,lw=0.9,mfc="none",label="vert (Up vs E)")
        lim=np.nanmax(np.abs([X,Y,Z]))*1.15
        ax.axhline(0,color="0.7",lw=0.6); ax.axvline(0,color="0.7",lw=0.6)
        ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect("equal")
        ax.set_title(f"{site} pos {pos:.0f} · {tag}",fontsize=8.5,weight="bold",color=SCOL[site])
        ax.tick_params(labelsize=6)
        if row==1 and col==0: ax.set_xlabel("E (nAm²)",fontsize=7); ax.set_ylabel("N / Up",fontsize=7)
axes[0][2].legend(fontsize=7,frameon=False,loc="upper right")
figD.suptitle("Figure S3 · Zijderveld orthogonal vector plots (representative specimens)",fontsize=12,weight="bold")
figD.tight_layout(rect=[0,0,1,0.96]); figD.savefig(f"{OUT}/zijderveld.png",dpi=180,facecolor="w"); plt.close(figD)

print("wrote af_demag_decay.png, grm_mdf_crossplot.png, stereonets.png, zijderveld.png")
for s in ("TAM","SC"): print(s,"series",len(SER[s]))
