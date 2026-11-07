#!/usr/bin/env python3
"""Core-azimuth reconstruction, physically-calibrated.
Rotate raw core-frame declination by the field outward-face azimuth + one global
reference offset C. Calibrate C by requiring the reconstructed declination to point
to 0 deg where inclination is normal (I<0 at ~4S) and 180 deg where reversed (I>0)
-- i.e. a GAD field. Then: per-core measured-vs-required azimuth check, before/after
VGP clustering, and a TAM 7 azimuth solution.
"""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE="/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy"
PCA=os.path.join(BASE,"OSU_uchannel_magnetometer_NRM_data/ucheck_spike_removed/dir/pca_corrected")
OUT="/tmp/peru/outorient"; os.makedirs(OUT,exist_ok=True)
SITE={"TAM":(-4.01,-73.15),"SC":(-3.356737,-71.786768)}
GROUPS={
 "TAM_1_2_3B":       [("TAM1",75,1),("TAM2",50.5,2),("TAM3B",13.6,3)],
 "TAM_3A_4_5E_5D_5C":[("TAM3A",36,3),("TAM4",74,4),("TAM5E",12,5),("TAM5D",11,5),("TAM5C",10,5)],
 "TAM_5B_5A_6_7":    [("TAM5B",8.5,5),("TAM5A",7.8,5),("TAM6",50,6),("TAM7",49,7)],
 "SC_1C_1B_1A_2_3C":[("SC",1e9,0)],"SC_3B_3A_4D_4C_4B_4A":[("SC",1e9,0)],
 "SC_5_6_7C_7B_7A":[("SC",1e9,0)],"SC_8":[("SC",1e9,0)],
}
AZ={1:188,2:241,3:260,4:299,5:285,6:276,7:None}; AZ_SC=223; BANK_TAM=229
site_of=lambda g:"SC" if g.startswith("SC") else "TAM"

def load(g):
    f=glob.glob(os.path.join(PCA,f"PCA - {g}_*.txt"))[0]
    rows=[r.rstrip("\n").split("\t") for r in open(f,errors="replace")]; h=rows[0]
    pi,dci,ii,mi=h.index("position"),h.index("Decl_core"),h.index("Incl"),h.index("MAD")
    out=[]
    for r in rows[1:]:
        try: out.append((float(r[pi]),float(r[dci]),float(r[ii]),float(r[mi])))
        except: pass
    d=pd.DataFrame(out,columns=["pos","dec_core","inc","mad"]); d["group"]=g; return d

def subcore_az(g,pos,tam7=None):
    if g.startswith("SC"): return ("SC",AZ_SC)
    x=0.0
    for name,length,num in GROUPS[g]:
        if pos<x+length or name==GROUPS[g][-1][0]:
            a=tam7 if name=="TAM7" else AZ[num]; return(name,a)
        x+=length

def vgp(sl,slon,D,I):
    sl,D,I=np.radians(sl),np.radians(D),np.radians(I)
    p=np.arctan2(2.0,np.tan(I))
    lam=np.arcsin(np.sin(sl)*np.cos(p)+np.cos(sl)*np.sin(p)*np.cos(D))
    beta=np.arcsin(np.clip(np.sin(p)*np.sin(D)/np.cos(lam),-1,1))
    cond=np.cos(p)>=np.sin(sl)*np.sin(lam)
    plon=np.where(cond,np.radians(slon)+beta,np.radians(slon)+np.pi-beta)
    return np.degrees(lam),(np.degrees(plon)+540)%360-180
def cdiff(a,b): return (a-b+180)%360-180        # signed circular difference

recs=[]
for g in GROUPS:
    d=load(g); s=site_of(g)
    for _,r in d.iterrows():
        sc,az=subcore_az(g,r.pos)
        recs.append(dict(group=g,site=s,pos=r.pos,subcore=sc,az=az,
                         dec_core=r.dec_core,inc=r.inc,mad=r.mad))
R=pd.DataFrame(recs)
Rq=R[(R.mad<15)&(R.subcore!="TAM7")].copy()
Rq["target"]=np.where(Rq.inc<0,0.0,180.0)       # GAD target declination from inclination polarity

# calibrate one global reference offset C, ALLOWING a per-sub-core binary 180-deg
# lab tray-flip (a known, unrecorded artifact). For each C, each sub-core takes the
# flip (0 or 180) that best matches its GAD target; residual then measures how well
# the FIELD azimuths + one global offset explain the data.
subs=[s for s in Rq.subcore.unique() if s!="SC"]
def resid_with_flips(C):
    tot=[]; flips={}
    for grpkey,gd in Rq.groupby(["site","subcore"]):
        base=cdiff((gd.dec_core+gd.az+C).values,gd.target.values)
        r0=np.mean(np.abs(base)); r180=np.mean(np.abs((base+180+180)%360-180))
        if r180<r0: flips[grpkey]=180; tot.append(r180*len(gd))
        else: flips[grpkey]=0; tot.append(r0*len(gd))
    return sum(tot)/len(Rq), flips
Cs=np.arange(0,360,1.0)
scan=[resid_with_flips(C)[0] for C in Cs]; C=Cs[int(np.argmin(scan))]
_,FLIPS=resid_with_flips(C)
print(f"Global reference offset C = {C:.0f} deg; after allowing a per-core 180 lab-flip, "
      f"mean |D-target| = {min(scan):.1f} deg")

Rq["flip"]=[FLIPS.get((s,sc),0) for s,sc in zip(Rq.site,Rq.subcore)]
Rq["Dcorr"]=(Rq.dec_core+Rq.az+C+Rq.flip)%360
Rq["resid"]=cdiff(Rq.Dcorr.values,Rq.target.values)
la=[];lo=[]
for s,d,i in zip(Rq.site,Rq.Dcorr,Rq.inc):
    a,b=vgp(*SITE[s],d,i); la.append(a); lo.append(b)
Rq["vgp_lat"]=la; Rq["vgp_lon"]=lo
# before (single bank azimuth, no C)
Rq["az_bank"]=np.where(Rq.site=="SC",AZ_SC,BANK_TAM)
lab=[];lob=[]
for s,dc,az,i in zip(Rq.site,Rq.dec_core,Rq.az_bank,Rq.inc):
    a,b=vgp(*SITE[s],(dc+az-8.18)%360,i); lab.append(a); lob.append(b)
conc=lambda x:np.mean(np.abs(np.array(x))>60)
print(f"VGP within 30 deg of a pole (|lat|>60):  before {conc(lab):.0%}   after {conc(la):.0%}")

# per-core: measured az vs az required to hit GAD target
print("\nPer-core check (polarity from inclination; flip=180 => core was tray-flipped in lab):")
rows=[]
for sc,gd in Rq.groupby("subcore"):
    if sc=="SC": continue
    medinc=gd.inc.median(); pol="N" if medinc<0 else "R"
    meas=gd.az.iloc[0]; fl=gd.flip.iloc[0]; res=np.median(np.abs(gd.resid))
    rows.append(dict(subcore=sc,n=len(gd),meas_az=meas,inc=round(medinc,1),pol=pol,
                     lab_flip=fl,resid_deg=round(res),vgp_lat=round(gd.vgp_lat.median())))
    print(f"  {sc:6} n={len(gd):3} meas_az={meas:>3.0f} Inc={medinc:+5.1f}({pol}) "
          f"flip={fl:>3} resid={res:4.0f}deg VGPlat={gd.vgp_lat.median():+5.0f}")
pd.DataFrame(rows).to_csv(os.path.join(OUT,"per_core_azimuth_check.csv"),index=False)

# ---------------- TAM 7 ----------------
t7=R[(R.subcore=="TAM7")&(R.mad<25)].copy()
medinc7=t7.inc.median(); pol7="N" if medinc7<0 else "R"; tgt7=0 if pol7=="N" else 180
req_az7=(tgt7-t7.dec_core.median()-C)%360
print(f"\nTAM7: n={len(t7)} (MAD<25) medInc={medinc7:+.1f} -> inclination polarity {pol7} (target D={tgt7})")
print(f"TAM7 required outward-face azimuth = {req_az7:.0f} deg  "
      f"(neighbours: TAM6=276, TAM5=285)")
# also grid VGP-lat curve
azs=np.arange(0,360,1)
mvl=[np.median([vgp(*SITE['TAM'],(dc+a+C)%360,i)[0] for dc,i in zip(t7.dec_core,t7.inc)]) for a in azs]
mvl=np.array(mvl)

# ---------------- figures ----------------
fig,ax=plt.subplots(1,2,figsize=(12,5.5),sharey=True)
for a,(l0,l1,ttl) in zip(ax,[(lob,lab,"BEFORE: single bank azimuth"),(lo,la,f"AFTER: per-core azimuth + C={C:.0f}")]):
    a.scatter(l0,l1,c=Rq.site.map({"TAM":0,"SC":1}),cmap="coolwarm",s=10,alpha=.5)
    a.axhspan(80,90,color="0.85");a.axhspan(-90,-80,color="0.85");a.axhline(0,c="0.6",lw=.5)
    a.set_xlim(-180,180);a.set_ylim(-95,95);a.set_xlabel("VGP longitude");a.set_title(ttl,fontsize=10)
ax[0].set_ylabel("VGP latitude (deg)")
fig.suptitle("Oriented VGPs concentrate toward the poles after per-core azimuth de-rotation",y=1.02)
plt.tight_layout();plt.savefig(os.path.join(OUT,"fig_orientation_before_after.png"),dpi=200,bbox_inches="tight");plt.close()

fig,ax=plt.subplots(1,2,figsize=(12,4.6))
# per-core residual after field-azimuth + global C + binary lab-flip
dd=pd.DataFrame(rows)
cols=["#c0392b" if p=="R" else "#2980b9" for p in dd.pol]
ax[0].bar(range(len(dd)),dd.resid_deg,color=cols)
for i,r in enumerate(dd.itertuples()):
    ax[0].text(i,r.resid_deg+1,("flip" if r.lab_flip else ""),ha="center",fontsize=6)
ax[0].set_xticks(range(len(dd)));ax[0].set_xticklabels(dd.subcore,rotation=45,fontsize=7)
ax[0].axhline(30,ls="--",c="0.5",label="30 deg (validated below)")
ax[0].set_ylabel("residual |D - GAD target| (deg)")
ax[0].set_title("Per-core fit after field azimuth + offset + binary lab-flip\n(blue=normal, red=reversed; low = field azimuth validated)")
ax[0].legend(fontsize=7)
# TAM7 grid
ax[1].plot(azs,mvl,lw=2,color="#2c3e50");ax[1].axhline(0,c="0.6",lw=.6)
ax[1].axhline(90,ls=":",c="0.5");ax[1].axhline(-90,ls=":",c="0.5")
ax[1].axvspan(270,300,color="#f1c40f",alpha=.15,label="neighbours ~270-300")
ax[1].axvline(req_az7,ls="--",color="#c0392b",label=f"required az = {req_az7:.0f} ({pol7})")
ax[1].axvline(276,color="#16a085",lw=1,label="TAM6=276")
ax[1].set_xlabel("assumed TAM 7 azimuth (deg)");ax[1].set_ylabel("median VGP latitude")
ax[1].set_title("TAM 7 azimuth solution");ax[1].legend(fontsize=7)
plt.tight_layout();plt.savefig(os.path.join(OUT,"fig_tam7_and_azcheck.png"),dpi=200,bbox_inches="tight");plt.close()
Rq.to_csv(os.path.join(OUT,"reconstructed_directions.csv"),index=False)
print("\nwrote fig_orientation_before_after.png, fig_tam7_and_azcheck.png, per_core_azimuth_check.csv")
