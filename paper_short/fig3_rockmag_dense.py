#!/usr/bin/env python3
"""Data-dense down-section paleomagnetic / rock-magnetic columns for TAM & SC.
Columns per site (vs reversal-aligned height): NRM intensity (raw+smoothed, log),
susceptibility, MDF, GRM ratio (I100/I70), inclination (raw+median), declination, MAD, polarity.
Derived from the u-channel DAT demag series + samples_corrected + OSU susceptibility runs."""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

BASE=os.environ.get("PEBAS_BASE","/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy")
PS=os.path.join(BASE,"polarity_synthesis"); DATD=os.path.join(BASE,"OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed")
OUT=os.environ.get("PEBAS_OUT","/mnt/user-data/outputs")

PZ=pd.read_csv(os.path.join(PS,"polarity_zones.csv")); POL={s:d for s,d in PZ.groupby("site")}
S=pd.read_csv(os.path.join(PS,"samples_corrected.csv"))
ORDER={"TAM":["TAM_3A_4_5E_5D_5C","TAM_5B_5A_6_7"],
       "SC":["SC_1C_1B_1A_2_3C","SC_3B_3A_4D_4C_4B_4A","SC_5_6_7C_7B_7A","SC_8"]}
MARKp={"TAM":(0.00,2.54),"SC":(1.60,5.00)}
def META(site):
    d=S[S.group.isin(ORDER[site])]; base=0.0; m={}
    for g in ORDER[site]:
        th=d[d.group==g].sed_depth_cm.max(); m[g]=(base,th); base+=th
    return m
def toH(site,cm):
    z=POL[site]; lo=z.top_cm.min(); hi=z.bottom_cm.max(); a,b=MARKp[site]
    return b-(cm-lo)/(hi-lo)*(b-a)
def rev_base(site):
    z=POL[site]; return toH(site,z[z.polarity=="R"].bottom_cm.max())
datum=rev_base("TAM"); offSC=datum-rev_base("SC"); OFF={"TAM":0.0,"SC":offSC}

# ---- DAT-derived per-position NRM intensity, MDF, GRM ratio ----
def dat_for(group):
    f=[x for x in glob.glob(os.path.join(DATD,"*.DAT")) if group.split("_")[0] in x and group.replace("SC","SE").split("_")[1] in x]
    # robust: match by tokens
    cand=glob.glob(os.path.join(DATD,"*.DAT"))
    best=None
    for c in cand:
        b=os.path.basename(c)
        if group[:6].replace("SC","S") in b or group in b or b.startswith(group.split("_")[0]):
            best=c
    return best
def per_pos(group):
    # find DAT file by matching leading tokens; tolerate the SE/SC prefix typo in exports
    cand=glob.glob(os.path.join(DATD,"*.DAT")); f=None
    toks=group.split("_")
    pref=[toks[0]]+(["SE"] if toks[0]=="SC" else [])   # SC_3B_3A_4D_4C_4B_4A DAT is misnamed SE_...
    for c in cand:
        b=os.path.basename(c)
        if any(p in b for p in pref) and toks[1] in b: f=c;break
    if f is None: return None
    d=pd.read_csv(f,sep="\t",header=0)
    d=d[d["Data type"].astype(str).str.contains("SAMPLE",na=False)].copy()
    d["AF"]=d[["AF X","AF Y","AF Z"]].max(axis=1)/10.0   # Oe->mT
    rows=[]
    for pos,dp in d.groupby("Position"):
        dp=dp.sort_values("AF"); af=dp["AF"].values; I=dp["Intensity"].values.astype(float)
        if len(af)<3 or I[0]<=0: continue
        nrm=I[0]; frac=I/nrm
        mdf=np.nan
        for k in range(1,len(af)):
            if frac[k]<=0.5: mdf=np.interp(0.5,[frac[k],frac[k-1]],[af[k],af[k-1]]); break
        def at(x):
            j=np.argmin(np.abs(af-x)); return I[j]
        grm=at(100)/at(70) if at(70)>0 else np.nan
        rows.append(dict(Position=pos,NRM=nrm,MDF=mdf,GRM=grm))
    return pd.DataFrame(rows)

def build(site):
    m=META(site); frames=[]
    for g in ORDER[site]:
        sg=S[S.group==g].copy().sort_values("uchannel_pos_cm")
        pp=per_pos(g)
        if pp is not None and len(pp):
            pps=pp.sort_values("Position")
            sg["NRM"]=np.interp(sg.uchannel_pos_cm,pps.Position,pps.NRM,left=np.nan,right=np.nan)
            sg["MDF"]=np.interp(sg.uchannel_pos_cm,pps.Position,pps.MDF,left=np.nan,right=np.nan)
            sg["GRM"]=np.interp(sg.uchannel_pos_cm,pps.Position,pps.GRM,left=np.nan,right=np.nan)
        else:
            sg["NRM"]=np.nan; sg["MDF"]=np.nan; sg["GRM"]=np.nan
        b,th=m[g]; sg["comp_cm"]=b+sg.sed_depth_cm; sg["H"]=toH(site,sg.comp_cm)+OFF[site]
        frames.append(sg)
    return pd.concat(frames).sort_values("H")
DAT={s:build(s) for s in ("TAM","SC")}

# ================= figure =================
LAN={"NRM":(0.0,0.9),"MDF":(1.05,1.75),"GRM":(1.9,2.6),
     "INC":(2.75,3.95),"DEC":(4.1,5.1),"MAD":(5.25,5.95),"POL":(6.1,6.4)}
BLOCK=7.5
NRM_LOG=(-7.0,-5.0)   # fixed, shared A/m log range so TAM & SC are comparable
# per-lane x-axis definition: key, value-lo, value-hi, tick values, tick labels
AXDEF=[("NRM",NRM_LOG[0],NRM_LOG[1],[-7,-6,-5],["10⁻⁷","10⁻⁶","10⁻⁵"]),
       ("MDF",0,60,[0,30,60],["0","30","60"]),
       ("GRM",0.5,1.6,[0.5,1.0,1.6],["0.5","1","1.6"]),
       ("INC",-45,45,[-40,0,40],["−40","0","40"]),
       ("DEC",0,360,[0,180,360],["0","180","360"]),
       ("MAD",0,20,[0,10,20],["0","10","20"])]
def smooth(s): return pd.Series(s).rolling(5,center=True,min_periods=1).median()
def sc(v,lane,lo,hi): return lane[0]+np.clip((v-lo)/(hi-lo),0,1)*(lane[1]-lane[0])
def draw(ax,site,xs):
    d=DAT[site]; z=POL[site]; off=OFF[site]
    zh=[toH(site,c)+off for c in list(z.top_cm)+list(z.bottom_cm)]; ylo,yhi=min(zh),max(zh)
    # NRM intensity (log, fixed shared range)
    ln=np.log10(d.NRM.replace(0,np.nan)); lo,hi=NRM_LOG
    ax.plot(xs+sc(ln,LAN["NRM"],lo,hi),d.H,color="0.6",lw=0.5)
    ax.plot(xs+sc(ln.rolling(7,center=True,min_periods=1).median(),LAN["NRM"],lo,hi),d.H,color="k",lw=1.1)
    # MDF (faint raw + smoothed)
    mlo,mhi=0,60
    ax.plot(xs+sc(d.MDF,LAN["MDF"],mlo,mhi),d.H,color="#1b7837",lw=0.4,alpha=0.35)
    ax.plot(xs+sc(smooth(d.MDF),LAN["MDF"],mlo,mhi),d.H,color="#1b7837",lw=1.1)
    # GRM ratio (>1 greigite): highlight (faint raw + smoothed)
    g1=sc(1.0,LAN["GRM"],0.5,1.6)
    ax.plot([xs+g1]*2,[ylo,yhi],color="0.7",lw=0.5,ls=":")
    ax.plot(xs+sc(d.GRM,LAN["GRM"],0.5,1.6),d.H,color="#c0392b",lw=0.4,alpha=0.35)
    ax.plot(xs+sc(smooth(d.GRM),LAN["GRM"],0.5,1.6),d.H,color="#c0392b",lw=1.1)
    # inclination
    exp=np.degrees(np.arctan(2*np.tan(np.radians(3.7))))
    for xv,c,st in [(0,"#c0392b","-"),(-exp,"0.7",":"),(exp,"0.7",":")]:
        ax.plot([xs+sc(xv,LAN["INC"],-45,45)]*2,[ylo,yhi],color=c,lw=0.6,ls=st)
    ax.scatter(xs+sc(d.inc.clip(-45,45),LAN["INC"],-45,45),d.H,s=3,c=np.where(d.inc<0,"#2166ac","#b2182b"),lw=0)
    ax.plot(xs+sc(d.inc_smooth.clip(-45,45),LAN["INC"],-45,45),d.H,color="k",lw=1.0)
    # declination
    ax.scatter(xs+sc(d.dec%360,LAN["DEC"],0,360),d.H,s=3,c="0.35",lw=0)
    # MAD (faint raw + smoothed)
    ax.plot([xs+sc(5,LAN["MAD"],0,20)]*2,[ylo,yhi],color="0.7",lw=0.5,ls=":")
    ax.plot(xs+sc(d.mad,LAN["MAD"],0,20),d.H,color="#6a51a3",lw=0.4,alpha=0.35)
    ax.plot(xs+sc(smooth(d.mad),LAN["MAD"],0,20),d.H,color="#6a51a3",lw=1.1)
    # polarity
    ax.add_patch(Rectangle((xs+LAN["POL"][0],ylo),LAN["POL"][1]-LAN["POL"][0],yhi-ylo,fc="white",ec="0.6",lw=.4,hatch="////",zorder=2))
    for _,r in z.iterrows():
        y0,y1=sorted((toH(site,r.top_cm)+off,toH(site,r.bottom_cm)+off))
        ax.add_patch(Rectangle((xs+LAN["POL"][0],y0),LAN["POL"][1]-LAN["POL"][0],y1-y0,fc=("k" if r.polarity=="N" else "w"),ec="k",lw=.4,zorder=3))
    # faint value gridlines behind data (all measurement lanes)
    for k,lo,hi,ticks,labs in AXDEF:
        for t in ticks:
            xt=xs+sc(t,LAN[k],lo,hi)
            ax.plot([xt,xt],[ylo,yhi],color="0.90",lw=0.4,zorder=0)
    # headers
    heads=[("NRM","NRM\n(A/m)"),("MDF","MDF\n(mT)"),("GRM","GRM\nI₁₀₀/I₇₀"),
           ("INC","inc (°)"),("DEC","dec (°)"),("MAD","MAD\n(°)"),("POL","pol")]
    for k,lab in heads:
        ax.text(xs+sum(LAN[k])/2,yhi+0.10,lab,ha="center",va="bottom",fontsize=6.5)
    # per-lane bottom x-axis: baseline + tick marks + labels
    ab=ylo-0.05                     # axis baseline
    for k,lo,hi,ticks,labs in AXDEF:
        x0,x1=xs+LAN[k][0],xs+LAN[k][1]
        ax.plot([x0,x1],[ab,ab],color="k",lw=0.7,zorder=4)
        for t,lb in zip(ticks,labs):
            xt=xs+sc(t,LAN[k],lo,hi)
            ax.plot([xt,xt],[ab,ab-0.055],color="k",lw=0.7,zorder=4)
            ax.text(xt,ab-0.10,lb,ha="center",va="top",fontsize=6.3)
    ax.text(xs+(LAN['NRM'][0]+LAN['POL'][1])/2,yhi+0.55,{"TAM":"Tamshiyacu (TAM)","SC":"Santa Carina (SC)"}[site],
            ha="center",fontsize=11,weight="bold")

fig,ax=plt.subplots(figsize=(15,9))
draw(ax,"TAM",0.0); draw(ax,"SC",BLOCK)
ymin=min(0,S[S.group.isin(ORDER['SC'])].sed_depth_cm.max()*0)  # placeholder
allH=pd.concat([DAT['TAM'].H,DAT['SC'].H]); ax.set_ylim(allH.min()-0.4,allH.max()+0.9)
ax.set_xlim(-0.5,BLOCK+LAN["POL"][1]+0.3)
ax.set_ylabel("Height in reversal-aligned frame (m)",fontsize=10)
ax.set_xticks([]); [ax.spines[s].set_visible(False) for s in ("top","right","bottom")]
leg=[Line2D([0],[0],color="0.6",lw=1,label="raw"),Line2D([0],[0],color="k",lw=1.4,label="smoothed"),
     Line2D([0],[0],color="#1b7837",lw=1,label="MDF"),
     Line2D([0],[0],color="#c0392b",lw=1,label="GRM ratio (>1 greigite)"),Line2D([0],[0],color="#6a51a3",lw=1,label="MAD"),
     Line2D([0],[0],marker="o",color="w",mfc="#2166ac",ms=6,label="inc<0 (normal)"),
     Line2D([0],[0],marker="o",color="w",mfc="#b2182b",ms=6,label="inc>0 (reversed)"),
     Patch(fc="k",ec="k",label="N"),Patch(fc="w",ec="k",label="R"),Patch(fc="white",ec="0.6",hatch="////",label="no data")]
fig.legend(handles=leg,loc="lower center",ncol=6,fontsize=7.5,frameon=False,bbox_to_anchor=(0.5,-0.02))
fig.suptitle("Paleomagnetic & rock-magnetic data down-section — TAM & SC",fontsize=12.5,weight="bold",y=0.98)
fig.savefig(os.path.join(OUT,"fig3_rockmag_dense.png"),dpi=190,bbox_inches="tight",facecolor="w")
print("wrote fig3_rockmag_dense.png; TAM n=%d SC n=%d"%(len(DAT['TAM']),len(DAT['SC'])))
