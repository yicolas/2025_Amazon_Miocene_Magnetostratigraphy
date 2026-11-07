#!/usr/bin/env python3
"""Equal-area (Schmidt) projection plots of the reconstructed directions + VGPs,
in the notebook's style (equal_area_xy from the notebook cell 35).
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/tmp/peru/outorient"
R=pd.read_csv(os.path.join(OUT,"reconstructed_directions.csv"))

def equal_area_xy(dec,inc):                      # notebook convention
    d=np.radians(np.asarray(dec,float)); i=np.radians(np.asarray(inc,float))
    Rr=np.sqrt(2.0)*np.sin((np.pi/2-np.abs(i))/2.0)
    return Rr*np.sin(d), Rr*np.cos(d)

def net(ax,title):
    th=np.linspace(0,2*np.pi,361); ax.plot(np.cos(th),np.sin(th),'k',lw=1.2)
    for r in (np.sqrt(2)*np.sin((np.pi/2-np.radians(a))/2) for a in (30,60)):
        ax.plot(r*np.cos(th),r*np.sin(th),color="0.8",lw=.6)
    ax.plot([0,0],[-1,1],color="0.85",lw=.5); ax.plot([-1,1],[0,0],color="0.85",lw=.5)
    ax.plot(0,0,'k+',ms=8)
    ax.text(0,1.04,"N",ha="center",fontsize=9); ax.text(1.04,0,"E",va="center",fontsize=9)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.15,1.15); ax.set_ylim(-1.15,1.15)
    ax.set_title(title,fontsize=10)

def plot_dirs(ax,df):
    dn=df[df.inc<0]; dr=df[df.inc>=0]          # up=normal(open), down=reversed(filled) at S hemi
    x,y=equal_area_xy(dn.Dcorr,dn.inc); ax.scatter(x,y,facecolors="none",edgecolors="#2980b9",s=14,lw=.7,label="Normal (inc<0, up)")
    x,y=equal_area_xy(dr.Dcorr,dr.inc); ax.scatter(x,y,color="#c0392b",s=14,label="Reversed (inc>0, down)")
    for D,I,mk,c in [(0,-8,'*','#1b4f72'),(180,8,'*','#7b241c')]:
        gx,gy=equal_area_xy(D,I); ax.scatter(gx,gy,marker=mk,s=200,color=c,edgecolors='k',zorder=5)

# ---- Figure 1: equal-area direction stereonets, before vs after, per site
fig,axes=plt.subplots(2,2,figsize=(10,10))
for row,site in enumerate(("TAM","SC")):
    d=R[R.site==site]
    # BEFORE: single bank azimuth (dec_core+bank+igrf), AFTER: reconstructed Dcorr
    bank=229 if site=="TAM" else 223
    dbefore=d.assign(Dcorr=(d.dec_core+bank-8.18)%360)
    net(axes[row,0],f"{site}  BEFORE (single bank azimuth)"); plot_dirs(axes[row,0],dbefore)
    net(axes[row,1],f"{site}  AFTER (per-core azimuth + lab-flip)"); plot_dirs(axes[row,1],d)
axes[0,1].legend(loc="upper right",fontsize=7,bbox_to_anchor=(1.15,1.1))
fig.suptitle("Equal-area stereonets of ChRM directions (geographic) — open=normal(up), filled=reversed(down); star=GAD",y=1.01,fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_stereonet_dirs.png"),dpi=200,bbox_inches="tight"); plt.close()

# ---- Figure 2: VGP polar equal-area (geographic N pole at centre)
def vgp_xy(vlat,vlon):
    colat=np.radians(90-np.asarray(vlat,float)); lon=np.radians(np.asarray(vlon,float))
    Rr=np.sqrt(2.0)*np.sin(colat/2.0)/np.sqrt(2.0)   # normalize so S pole -> r=1
    return Rr*np.sin(lon), Rr*np.cos(lon)
fig,axes=plt.subplots(1,2,figsize=(11,5.6))
for ax,site in zip(axes,("TAM","SC")):
    th=np.linspace(0,2*np.pi,361); ax.plot(np.cos(th),np.sin(th),'k',lw=1.2)       # equator
    ax.plot(0,0,'k+',ms=9); ax.text(0,1.05,"N pole",ha="center",fontsize=8); ax.text(0,-1.12,"S pole (edge)",ha="center",fontsize=7,color="0.4")
    for r in (0.5,):
        ax.plot(r*np.cos(th),r*np.sin(th),color="0.85",lw=.6)
    d=R[R.site==site]; dn=d[d.vgp_lat>0]; dr=d[d.vgp_lat<=0]
    x,y=vgp_xy(dn.vgp_lat,dn.vgp_lon); ax.scatter(x,y,facecolors="none",edgecolors="#2980b9",s=16,lw=.7,label="N-pole VGPs")
    x,y=vgp_xy(dr.vgp_lat,dr.vgp_lon); ax.scatter(x,y,color="#c0392b",s=16,label="S-pole VGPs")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.15,1.15); ax.set_ylim(-1.2,1.15)
    ax.set_title(f"Site {site} — VGPs (equal-area, N pole centre)",fontsize=10); ax.legend(fontsize=7,loc="lower right")
fig.suptitle("VGP equal-area projection (reconstructed): normal VGPs toward centre (N pole), reversed toward edge (S pole)",y=1.02,fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"fig_vgp_polar.png"),dpi=200,bbox_inches="tight"); plt.close()
print("wrote fig_stereonet_dirs.png, fig_vgp_polar.png")
