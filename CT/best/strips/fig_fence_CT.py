#!/usr/bin/env python3
"""Consolidated ODP-style fence with CT core-face lane.
Lanes per site: CT | lithology | inclination | polarity, reversal-aligned.
CT strips assigned by channel->group (see README_CT_strips.md); DRAFT registration."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Patch
from matplotlib.lines import Line2D
from PIL import Image

PS="/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy/polarity_synthesis"
CT="/mnt/user-data/outputs/ct_strips"; OUT="/tmp/peru/fig_build"

FAC={"lignite/coal":"#242424","organic clay":"#5c5138","mudstone (clay)":"#8fb0c0",
     "silty clay":"#a9c3d1","sandy clay":"#cbb48f","shell bed":"#e9d9a6"}
TAM_BEDS=pd.DataFrame([(0.00,0.20,"mudstone (clay)"),(0.20,0.26,"lignite/coal"),
 (0.26,2.06,"silty clay"),(2.06,2.66,"mudstone (clay)"),(2.66,2.76,"lignite/coal"),
 (2.76,3.40,"mudstone (clay)"),(3.40,3.60,"organic clay"),(3.60,3.78,"lignite/coal")],
 columns=["base_m","top_m","facies"])
SC_BEDS=pd.DataFrame([(0.00,0.45,"mudstone (clay)"),(0.45,0.65,"lignite/coal"),
 (0.65,0.90,"organic clay"),(0.90,1.80,"sandy clay"),(1.80,2.50,"mudstone (clay)"),
 (2.50,4.90,"silty clay"),(4.90,5.00,"lignite/coal")],columns=["base_m","top_m","facies"])
MARK={"TAM":dict(pmag=(0.00,2.54),nodata_top=(2.54,3.78),puru=2.71),
      "SC":dict(pmag=(1.60,5.00))}
PZ=pd.read_csv(os.path.join(PS,"polarity_zones.csv")); POL={s:d for s,d in PZ.groupby("site")}
S=pd.read_csv(os.path.join(PS,"samples_corrected.csv"))
ORD={"TAM":["TAM_3A_4_5E_5D_5C","TAM_5B_5A_6_7"],
     "SC":["SC_1C_1B_1A_2_3C","SC_3B_3A_4D_4C_4B_4A","SC_5_6_7C_7B_7A","SC_8"]}
def comp(site):
    d=S[S.group.isin(ORD[site])].copy(); base=0; out=[]
    for g in ORD[site]:
        dg=d[d.group==g].copy(); th=dg.sed_depth_cm.max(); dg["cc"]=base+dg.sed_depth_cm; base+=th; out.append(dg)
    return pd.concat(out)
DIR={s:comp(s) for s in ("TAM","SC")}
def toH(site,cm):
    z=POL[site]; lo=z.top_cm.min(); hi=z.bottom_cm.max(); a,b=MARK[site]["pmag"]; return b-(cm-lo)/(hi-lo)*(b-a)
def rev_base(site):
    z=POL[site]; return toH(site,z[z.polarity=="R"].bottom_cm.max())
datum=rev_base("TAM"); offSC=datum-rev_base("SC"); OFF={"TAM":0.0,"SC":offSC}

# ---- CT strip -> group -> field interval (DRAFT assignment; verify band/orient) ----
def load(fn):
    im=Image.open(os.path.join(CT,fn)).convert("L"); return np.asarray(im)
CTMAP={  # (base_m, top_m, strip file, flip_vertical)
 "TAM":[(0.00,1.13,"SE1_0014_IM14_b4_99cm.png",False),   # Group 3 (base)
        (1.13,2.54,"SE1_0014_IM14_b3_97cm.png",False),   # Group 2 (mid)
        (2.54,3.78,"SE1_0014_IM14_b2_96cm.png",False)],  # Group 1 (top; Purussaurus)
 "SC":[(4.35,5.00,"SE1_0014_IM14_b5_99cm.png",False),    # Group 4 (top)  [board1]
       (3.09,4.35,"SE2_0011_IM11_b1_59cm.png",False),    # Group 5        [board2]
       (1.87,3.09,"SE2_0011_IM11_b2_61cm.png",False),    # Group 6
       (1.65,1.87,"SE2_0011_IM11_b3_61cm.png",False)],   # Group 7
}

CTX=(-0.95,-0.55); LITHO=(0.0,0.42); INC=(0.55,1.72); POLX=(1.80,2.05)
INC0,INC1=-45,45
def xinc(v): return INC[0]+(v-INC0)/(INC1-INC0)*(INC[1]-INC[0])
LANE=3.2

def draw_site(ax,site,xs):
    off=OFF[site]; beds=(TAM_BEDS if site=="TAM" else SC_BEDS).copy()
    beds["b"]=beds.base_m+off; beds["t"]=beds.top_m+off; mk=MARK[site]
    ymin=beds.b.min(); ymax=beds.t.max()
    # CT lane
    for (b0,t0,fn,flip) in CTMAP[site]:
        img=load(fn)
        if flip: img=img[::-1]
        yb,yt=b0+off,t0+off
        ax.imshow(img,extent=[xs+CTX[0],xs+CTX[1],yb,yt],origin="upper",
                  aspect="auto",cmap="gray",zorder=3,vmin=0,vmax=255)
        ax.add_patch(Rectangle((xs+CTX[0],yb),CTX[1]-CTX[0],yt-yb,fill=False,ec="0.5",lw=0.4,zorder=4))
    # lithology
    for _,r in beds.iterrows():
        ax.add_patch(Rectangle((xs+LITHO[0],r.b),LITHO[1]-LITHO[0],r.t-r.b,facecolor=FAC[r.facies],ec="k",lw=0.5,zorder=3))
    # inclination
    exp=np.degrees(np.arctan(2*np.tan(np.radians(3.7))))
    for xv,st,c,lw in [(0,"-","#c0392b",1.0),(-exp,":","0.6",.7),(exp,":","0.6",.7)]:
        ax.plot([xs+xinc(xv)]*2,[ymin,ymax],color=c,lw=lw,ls=st,zorder=2)
    d=DIR[site].copy(); d["H"]=toH(site,d.cc)+off; d=d.sort_values("H")
    ax.plot(xs+xinc(d.inc.clip(INC0,INC1)),d.H,color="0.6",lw=0.6,zorder=3)
    ax.plot(xs+xinc(d.inc_smooth.clip(INC0,INC1)),d.H,color="k",lw=1.3,zorder=4)
    ax.scatter(xs+xinc(d.inc.clip(INC0,INC1)),d.H,s=4,c=np.where(d.inc<0,"#2166ac","#b2182b"),zorder=5,lw=0)
    # polarity
    z=POL[site]
    ax.add_patch(Rectangle((xs+POLX[0],ymin),POLX[1]-POLX[0],ymax-ymin,facecolor="#d9d9d9",ec="k",lw=.4,hatch="xxx",zorder=3))
    for _,r in z.iterrows():
        y0,y1=sorted((toH(site,r.top_cm)+off,toH(site,r.bottom_cm)+off))
        ax.add_patch(Rectangle((xs+POLX[0],y0),POLX[1]-POLX[0],y1-y0,facecolor=("k" if r.polarity=="N" else "w"),ec="k",lw=.4,zorder=4))
    gaps=[]
    if site=="TAM": gaps.append((mk["nodata_top"][0]+off,mk["nodata_top"][1]+off))
    for _,r in beds.iterrows():
        if "lignite" in r.facies: gaps.append((r.b,r.t))
    for g0,g1 in gaps:
        ax.add_patch(Rectangle((xs+POLX[0],g0),POLX[1]-POLX[0],g1-g0,facecolor="#d9d9d9",ec="k",lw=.4,hatch="xxx",zorder=6))
    for _,r in z.iterrows():
        if r.polarity=="R":
            yc=(toH(site,r.top_cm)+toH(site,r.bottom_cm))/2+off
            ax.annotate("C5AAr",(xs+POLX[1],yc),xytext=(xs+POLX[1]+0.06,yc),fontsize=7.5,color="#c0392b",va="center",weight="bold")
    if mk.get("puru") is not None:
        h=mk["puru"]+off; ax.plot(xs+0.21,h,marker=(5,1,0),ms=15,mfc="#8e44ad",mec="k",mew=.8,zorder=8)
        ax.annotate("Purussaurus",(xs+0.21,h),xytext=(xs+0.44,h+0.33),fontsize=7,color="#8e44ad",va="center",
                    arrowprops=dict(arrowstyle="-",color="#8e44ad",lw=.8))
    # headers
    yh=ymax+0.16
    ax.text(xs+sum(CTX)/2,yh,"CT",ha="center",fontsize=7.5)
    ax.text(xs+sum(LITHO)/2,yh,"lith.",ha="center",fontsize=7)
    ax.text(xs+sum(INC)/2,yh,"inclination (°)",ha="center",fontsize=7.5)
    ax.text(xs+sum(POLX)/2,yh,"pol.",ha="center",fontsize=7)
    ax.text(xs+sum(INC)/2,ymax+0.40,site,ha="center",fontsize=12,weight="bold")
    for v in (-40,-20,0,20,40):
        ax.text(xs+xinc(v),ymin-0.16,f"{v}",ha="center",fontsize=6)

fig,ax=plt.subplots(figsize=(12,9))
draw_site(ax,"TAM",0.0); draw_site(ax,"SC",LANE)
ax.axhline(datum,color="#c0392b",lw=1.4,ls="--",zorder=9)
ax.text(-1.05,datum,"N→R datum\n(C5AAr base)",fontsize=7.5,color="#c0392b",va="center",ha="right",weight="bold")
ymin=min(0,SC_BEDS.base_m.min()+offSC)-0.5; ymax=max(TAM_BEDS.top_m.max(),SC_BEDS.top_m.max()+offSC)+0.7
ax.set_ylim(ymin,ymax); ax.set_xlim(-1.6,LANE+POLX[1]+0.5)
ax.set_ylabel("Height above section base (m)",fontsize=10)
ax.set_yticks(np.arange(0,ymax,0.5)); ax.set_xticks([])
for sp in ("top","right","bottom"): ax.spines[sp].set_visible(False)
leg=[Patch(fc=c,ec="k",label=l) for l,c in FAC.items()]
leg+=[Line2D([0],[0],color="k",lw=1.4,label="ChRM inc (7-cm median)"),
      Line2D([0],[0],marker="o",color="w",mfc="#2166ac",ms=6,label="normal (inc<0)"),
      Line2D([0],[0],marker="o",color="w",mfc="#b2182b",ms=6,label="reversed (inc>0)"),
      Patch(fc="k",ec="k",label="N"),Patch(fc="w",ec="k",label="R"),
      Patch(fc="#d9d9d9",ec="k",hatch="xxx",label="no polarity data"),
      Line2D([0],[0],marker=(5,1,0),color="w",mfc="#8e44ad",mec="k",ms=11,label="Purussaurus"),
      Patch(fc="0.5",ec="k",label="CT core-face")]
fig.legend(handles=leg,loc="lower center",ncol=5,fontsize=7.5,frameon=False,bbox_to_anchor=(0.5,-0.05))
fig.suptitle("TAM & SC fence — CT · lithology · inclination · polarity, aligned on the shared C5AAr reversal",fontsize=12,weight="bold",y=0.98)
fig.text(0.5,0.005,"DRAFT CT registration: channel→group by board labels; CT strips scaled to each group's field interval; band pick & up/down orientation to verify.",
         ha="center",fontsize=7,style="italic",color="0.35")
fig.savefig(os.path.join(OUT,"fig_fence_CT.png"),dpi=200,bbox_inches="tight",facecolor="w")
print("wrote fig_fence_CT.png datum=%.2f offSC=%.2f"%(datum,offSC))
