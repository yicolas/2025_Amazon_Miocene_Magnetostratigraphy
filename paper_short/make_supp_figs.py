import os
#!/usr/bin/env python3
"""Supplement figures S1 (GRM/MDF crossplot), S2 (reversal stereonet), S3 (Zijderveld)."""
import glob, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.environ.get("PEBAS_BASE") or glob.glob("/mnt/user-data/uploads/**/2025_Amazon_Miocene_Magnetostratigraphy", recursive=True)[0]
INTDIR = os.path.join(BASE, "OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed/int")
DATDIR = os.path.join(BASE, "OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed")
SYN    = os.path.join(BASE, "polarity_synthesis/samples_corrected.csv")
OUT    = os.environ.get("PEBAS_OUT","/tmp/peru/supp/figs")
os.makedirs(OUT, exist_ok=True)

STEPS_MT = np.array([0,10,20,25,30,35,40,50,70,90,100], float)  # header 0..1000 = *0.1 mT

def site_of(name):
    n = name.upper()
    if n.startswith("TAM"): return "TAM"
    return "SC"          # SE_* is an SC typo group

# ---------------------------------------------------------------- S1: GRM / MDF
rows=[]
for f in sorted(glob.glob(os.path.join(INTDIR,"*.int"))):
    df = pd.read_csv(f, sep="\t")
    df.columns=[c.strip() for c in df.columns]
    intcols=[c for c in df.columns if c.strip().startswith("Int")]
    if len(intcols)>len(STEPS_MT):            # stray trailing empty column
        intcols=intcols[:len(STEPS_MT)]
    if len(intcols)!=len(STEPS_MT):
        print("  skip (cols=%d):"%len(intcols), os.path.basename(f)); continue
    run=os.path.basename(f).split("_NRM")[0]
    for _,r in df.iterrows():
        I=np.array([float(r[c]) for c in intcols])
        if I[0]<=0 or np.any(~np.isfinite(I)): continue
        Inorm=I/I[0]
        # MDF: first field where normalized intensity <= 0.5 (linear interp)
        mdf=np.nan
        for k in range(1,len(Inorm)):
            if Inorm[k]<=0.5:
                x0,x1=STEPS_MT[k-1],STEPS_MT[k]; y0,y1=Inorm[k-1],Inorm[k]
                mdf=x0+(0.5-y0)*(x1-x0)/(y1-y0) if y1!=y0 else x1
                break
        if not np.isfinite(mdf): mdf=STEPS_MT[-1]
        # GRM index: max intensity above 50 mT relative to the running minimum up to 50 mT
        hi=I[STEPS_MT>=50]; base=I[STEPS_MT<=50].min()
        grm=hi.max()/base if base>0 else np.nan
        rows.append(dict(run=run, site=site_of(run), mdf=mdf, grm=grm))
S1=pd.DataFrame(rows).replace([np.inf,-np.inf],np.nan).dropna()

fig,ax=plt.subplots(figsize=(6.6,5.2))
# data-driven classification matching the text's physical criteria:
#   greigite / GRM-affected  = high coercivity (MDF>=50 mT) OR high-field intensity rise (GRM index>=1.5)
#   magnetite (soft)         = otherwise (MDF ~10-35 mT, monotone decay)
S1["greig"]=(S1.mdf>=50)|(S1.grm>=1.5)
ax.axvspan(50,105,color="#b2182b",alpha=.05,zorder=0)
g0=S1[~S1.greig]; g1=S1[S1.greig]
ax.scatter(g0.mdf,g0.grm,s=16,alpha=.5,edgecolor="none",color="#2166ac",
           label=f"magnetite, soft (n={len(g0)})")
ax.scatter(g1.mdf,g1.grm,s=18,alpha=.6,edgecolor="none",color="#b2182b",
           label=f"greigite / GRM-affected (n={len(g1)})")
ax.axhline(1.5,ls=":",c="grey",lw=1); ax.text(6,1.55,"GRM-index = 1.5",fontsize=6.5,color="grey")
ax.text(77,ax.get_ylim()[1]*0.02,"high-coercivity field (MDF $\\geq$ 50 mT)",fontsize=7,color="#b2182b",ha="center")
ax.set_xlabel("Median destructive field, MDF (mT)")
ax.set_ylabel("GRM index  (max intensity$_{>50\\,mT}$ / min intensity$_{\\leq50\\,mT}$)")
ax.set_title("Figure S1. Coercivity vs. gyroremanence (GRM) behaviour",fontsize=10,loc="left")
ax.legend(fontsize=8,frameon=False,loc="upper right")
ax.grid(alpha=.2)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"figS1_grm_mdf.png"),dpi=200); plt.close(fig)
print("S1:",len(S1),"positions; GRM>1:",int((S1.grm>1).sum()))

# ---------------------------------------------------------------- S2: reversal stereonet
syn=pd.read_csv(SYN)
def eqarea(dec,inc):
    inc=np.asarray(inc,float); dec=np.deg2rad(np.asarray(dec,float))
    R=np.sqrt(2)*np.sin(np.deg2rad(90-np.abs(inc))/2)
    return R*np.sin(dec), R*np.cos(dec)
# single u-channel run spanning the datum reversal: relative declination is internally
# consistent within one run, so the normal and reversed populations are meaningfully antipodal.
GRP="SC_3B_3A_4D_4C_4B_4A"
sg=syn[syn.group==GRP]
fig,ax=plt.subplots(figsize=(5.8,6.4))
th=np.linspace(0,2*np.pi,361)
ax.plot(np.sin(th),np.cos(th),color="k",lw=1.3)
for t in range(0,360,10):
    tr=np.deg2rad(t); L=0.025
    ax.plot([np.sin(tr)*(1-L),np.sin(tr)],[np.cos(tr)*(1-L),np.cos(tr)],color="k",lw=.6)
ax.plot([0,0],[-1,1],color="k",lw=.5,alpha=.25); ax.plot([-1,1],[0,0],color="k",lw=.5,alpha=.25)
norm=sg[sg.inc<0]; rev=sg[sg.inc>=0]
xn,yn=eqarea(norm.dec,norm.inc); xr,yr=eqarea(rev.dec,rev.inc)
# mean directions
def meandir(d,i):
    d=np.deg2rad(d); i=np.deg2rad(i)
    x=np.mean(np.cos(i)*np.cos(d)); y=np.mean(np.cos(i)*np.sin(d)); z=np.mean(np.sin(i))
    D=np.rad2deg(np.arctan2(y,x))%360; I=np.rad2deg(np.arctan2(z,np.hypot(x,y)))
    return D,I
Dn,In=meandir(norm.dec.values,norm.inc.values); Dr,Ir=meandir(rev.dec.values,rev.inc.values)
mxn,myn=eqarea([Dn],[In]); mxr,myr=eqarea([Dr],[Ir])
ax.scatter(xn,yn,facecolors="none",edgecolors="#2166ac",s=26,lw=1.0,label=f"Normal (inc<0, up)  n={len(norm)}")
ax.scatter(xr,yr,facecolors="#b2182b",edgecolors="k",s=26,lw=.3,alpha=.85,label=f"Reversed (inc>0, down)  n={len(rev)}")
ax.scatter(mxn,myn,marker="*",s=260,facecolor="#08306b",edgecolor="k",lw=.8,zorder=5)
ax.scatter(mxr,myr,marker="*",s=260,facecolor="#67000d",edgecolor="k",lw=.8,zorder=5)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.12,1.12); ax.set_ylim(-1.15,1.18)
ax.text(0,1.08,"N (relative declination)",ha="center",fontsize=8.5)
ax.set_title(f"Figure S2. Reversal in equal-area projection (run {GRP.replace('_',' ')})",
             fontsize=9.5,loc="left")
ax.legend(fontsize=8,frameon=False,loc="upper left",bbox_to_anchor=(-0.02,1.0))
ax.text(0,-1.12,"Stars = Fisher mean directions. Within one unoriented u-channel run relative declination\n"
        "is internally consistent, so the ~180$^\\circ$ swing between populations is a meaningful reversal;\n"
        "absolute declination is not defined. Polarity is read from inclination sign.",
        ha="center",fontsize=6.8,style="italic",color="dimgray")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"figS2_stereonet.png"),dpi=200,bbox_inches="tight"); plt.close(fig)
print("S2 (%s): normal %d reversed %d  meanN=(%.0f,%.1f) meanR=(%.0f,%.1f)"%(GRP,len(norm),len(rev),Dn,In,Dr,Ir))

# ---------------------------------------------------------------- S3: Zijderveld
def load_dat(path):
    df=pd.read_csv(path,sep="\t")
    df.columns=[c.strip() for c in df.columns]
    # AF level (mT) from AF X/Y/Z (they degauss all three to same peak); take max
    for c in ["AF X","AF Y","AF Z"]:
        df[c]=pd.to_numeric(df[c],errors="coerce")
    df["af_mt"]=df[["AF X","AF Y","AF Z"]].max(axis=1)/10.0   # *0.1 mT
    df["pos"]=pd.to_numeric(df["Position"],errors="coerce")
    for c in ["X intensity","Y intensity","Z intensity"]:
        df[c]=pd.to_numeric(df[c],errors="coerce")
    df["dtype"]=df["Data type"].astype(str).str.strip()
    return df

def series(sub):
    sub=sub.sort_values("af_mt")
    x=sub["X intensity"].values; y=sub["Y intensity"].values; z=sub["Z intensity"].values
    mt=sub["af_mt"].values
    N=np.sqrt(x**2+y**2+z**2)
    return mt,x,y,z,N

def zij(ax,sub,title,color):
    mt,x,y,z,N=series(sub)
    s=1.0/np.nanmax(N); x,y,z=x*s,y*s,z*s
    ax.plot(y,x,"-o",color=color,ms=4,lw=1,label="horizontal (E,N)")
    ax.plot(y,z,"-s",color="grey",ms=3.5,lw=1,mfc="white",label="vertical (E,Down)")
    ax.axhline(0,color="k",lw=.6); ax.axvline(0,color="k",lw=.6)
    ax.set_aspect("equal")
    lim=max(np.nanmax(np.abs(np.concatenate([x,y,z])))*1.15,0.05)
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    ax.set_title(title,fontsize=8.5)
    ax.set_xlabel("E",fontsize=7); ax.set_ylabel("N  /  Down",fontsize=7)
    ax.tick_params(labelsize=6)
    for xi,yi,m in zip(y,x,mt):
        if m in (0,20,50,100):
            ax.annotate(f"{int(m)}",(xi,yi),fontsize=5.5,color="dimgray",xytext=(2,2),textcoords="offset points")

def metrics(df):
    """per position: fraction of magnitude removed, GRM rise, ChRM inclination (30-70 mT)."""
    out=[]
    for p,sub in df.groupby("pos"):
        if sub["af_mt"].nunique()<9: continue
        mt,x,y,z,N=series(sub)
        if not np.isfinite(N[0]) or N[0]<=0: continue
        frem=1-N[-1]/np.nanmax(N)
        base=N[mt<=50].min(); rise=(N[mt>50].max()/base) if base>0 else np.nan
        m=(mt>=30)&(mt<=70)
        inc=np.rad2deg(np.arctan2(np.mean(z[m]),np.hypot(np.mean(x[m]),np.mean(y[m]))))
        out.append(dict(pos=p,frem=frem,rise=rise,inc=inc,nrm=N[0]))
    return pd.DataFrame(out)

sc1=load_dat(os.path.join(DATDIR,"SC_1C_1B_1A_2_3C_NRM_01212026.uchecked.DAT")); sc1=sc1[sc1.dtype!="LEADER"]
scr=load_dat(os.path.join(DATDIR,"SE_3B_3A_4D_4C_4B_4A_NRM_01202026.uchecked.DAT")); scr=scr[scr.dtype!="LEADER"]  # 'SE' = SC typo in filename
tam=load_dat(os.path.join(DATDIR,"TAM_3A_4_5E_5D_5C_NRM_01192026.uchecked.DAT")); tam=tam[tam.dtype!="LEADER"]

m1=metrics(sc1); m2=metrics(tam); m3=metrics(scr)
# panel A: cleanest soft normal magnetite (high frem, low rise, negative inc)
a=m1[(m1.rise<1.3)&(m1.inc<0)].sort_values("frem",ascending=False).iloc[0]
# panel B: greigite/GRM in TAM grp3 (largest high-field rise)
b=m2.sort_values("rise",ascending=False).iloc[0]
# panel C: cleanest reversed-polarity position (positive inc), high frem
cc=m3[(m3.inc>0)&(m3.rise<1.5)].sort_values("frem",ascending=False).iloc[0]

fig,axes=plt.subplots(1,3,figsize=(12.6,4.6))
zij(axes[0],sc1[sc1.pos==a.pos],f"SC normal (magnetite)  pos {a.pos:.0f} cm\n{a.frem*100:.0f}% removed by 100 mT, inc {a.inc:+.0f}°","#2166ac")
zij(axes[1],tam[tam.pos==b.pos],f"TAM Grp 3 greigite (GRM)  pos {b.pos:.0f} cm\nintensity rises ×{b.rise:.1f} above 50 mT","#b2182b")
zij(axes[2],scr[scr.pos==cc.pos],f"SC reversed  pos {cc.pos:.0f} cm\n{cc.frem*100:.0f}% removed, inc {cc.inc:+.0f}°","#4d9221")
p1,p2,p3=a.pos,b.pos,cc.pos
axes[0].legend(fontsize=6.5,frameon=False,loc="lower right")
fig.suptitle("Figure S3. Representative orthogonal (Zijderveld) demagnetization diagrams  —  intensities normalized; labels = AF field (mT)",
             fontsize=9.5,x=0.5)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(os.path.join(OUT,"figS3_zijderveld.png"),dpi=200); plt.close(fig)
print("S3 positions:",p1,p2,p3)
print("DONE ->",OUT)
