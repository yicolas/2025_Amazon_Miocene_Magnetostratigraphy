import os
#!/usr/bin/env python3
"""Figure 2 - master magnetostratigraphic correlation (Pebas Fm, Iquitos; short & long manuscripts).
Left (metres, reversal-aligned): TAM & SC, each = [pystrat litho | feature lane (SVG symbols, Gingras 2002
parasequences) | inclination | polarity | ln(Fe/coh) | ln(S/coh)].
Right: La Venta composite section (0-625 m, Mora-Rojas et al. 2023; Monkey & Fish beds) | GPTS (ATNTS2012) |
La Venta U-Pb age model (Ma); C5AAr tied Peru -> La Venta -> GPTS by ribbons, La Venta reversals -> GPTS by lines.
Also computes the Monte-Carlo inferred age of the Tamshiyacu Purussaurus."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, Polygon as MplPolygon, Wedge, Ellipse
from matplotlib.lines import Line2D
import pystrat, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_symbols as ss   # SVG accessory symbols (./symbols/*.svg)

BASE=os.environ.get("PEBAS_BASE","/mnt/user-data/uploads/GitHub/2025_Amazon_Miocene_Magnetostratigraphy")
PS=os.path.join(BASE,"polarity_synthesis")
XRF=os.path.join(BASE,".notebook_append/xrf_tidy.csv")
OUT=os.environ.get("PEBAS_OUT","/tmp/peru/fig_build"); LIT="/tmp/peru/fig_build/litho"; os.makedirs(LIT,exist_ok=True)

# ---------------- lithology (pystrat) ----------------
FAC=pd.DataFrame([
 ("lignite/coal","#242424",0.24,0,"F3A"),("organic clay","#5c5138",0.32,0,"F1B"),
 ("mudstone (clay)","#8fb0c0",0.40,616,"F1A/B"),("silty clay","#a9c3d1",0.52,613,"F2C/D"),
 ("sandy clay","#cbb48f",0.66,607,"F2A/B"),("shell bed","#e9d9a6",0.58,627,"F3B")],
 columns=["label","color","width","swatch","fcode"])
STYLE=pystrat.Style(FAC.label.tolist(),FAC.color.tolist(),FAC.width.tolist(),swatch_values=FAC.swatch.tolist())
def leglabel(r): return f"{r.label} ({r.fcode})" if r.fcode else r.label
LWIDTH=0.72  # pystrat x-extent (grain-size width units) reserved in the litho lane
def render_litho(beds, tag):
    th=(beds.top_m-beds.base_m).tolist(); tot=beds.top_m.max()
    sec=pystrat.Section(th,beds.facies.tolist())
    ipu=2.0/LWIDTH  # inches per data-unit (equal x/y scale, same for both sites)
    fig,ax=plt.subplots(figsize=(LWIDTH*ipu,tot*ipu)); fig.subplots_adjust(0,0,1,1)
    sec.plot(STYLE,ax=ax,xticks=False); ax.set_xlim(0,LWIDTH); ax.set_ylim(0,tot); ax.axis("off")
    p=os.path.join(LIT,f"litho_{tag}.png"); fig.savefig(p,dpi=150,transparent=True); plt.close(fig)
    return p

TAM_BEDS=pd.read_csv(os.path.join(BASE,"stratigraphy/beds_TAM.csv")) if os.path.exists(os.path.join(BASE,"stratigraphy/beds_TAM.csv")) else pd.DataFrame([
 (0.00,0.20,"mudstone (clay)"),(0.20,0.26,"lignite/coal"),(0.26,2.06,"silty clay"),
 (2.06,2.66,"mudstone (clay)"),(2.66,2.76,"lignite/coal"),(2.76,3.40,"mudstone (clay)"),
 (3.40,3.60,"organic clay"),(3.60,3.78,"lignite/coal")],columns=["base_m","top_m","facies"])
SC_BEDS=pd.read_csv(os.path.join(BASE,"stratigraphy/beds_SC.csv")) if os.path.exists(os.path.join(BASE,"stratigraphy/beds_SC.csv")) else pd.DataFrame([
 (0.00,0.45,"mudstone (clay)"),(0.45,0.65,"lignite/coal"),(0.65,0.90,"organic clay"),
 (0.90,1.80,"sandy clay"),(1.80,2.50,"mudstone (clay)"),(2.50,4.90,"silty clay"),
 (4.90,5.00,"lignite/coal")],columns=["base_m","top_m","facies"])
BEDS={"TAM":TAM_BEDS,"SC":SC_BEDS}
LIMG={s:render_litho(BEDS[s],s) for s in ("TAM","SC")}

# ---- La Venta composite standard section (CS) of Mora-Rojas et al. 2023, digitised from suppl. s6 (lithology,
# lithofacies codes) and s15 (reinterpreted polarity) by age_model_assessment/laventa/build_laventa_cs.py.
# CS metres are the authors' metres: U-Pb samples (Table 4), fossil localities (s10) and StLs all use them,
# so the whole column - including the Monkey Beds - sits on one continuous scale.
LV_TOP=625.0
LV_FAC=pd.DataFrame([("claystone","#8f8a84",0.36,616,""),("mudstone","#a9a9a0",0.40,616,""),
                     ("siltstone","#c9c7b8",0.44,613,""),("mudstone_red","#c47a68",0.40,616,""),
                     ("sandstone","#e0c877",0.62,607,""),("conglomerate","#c69a5c",0.74,607,"")],
                    columns=["label","color","width","swatch","fcode"])
LV_GW={"clay":0.36,"silt":0.40,"very fine sand":0.48,"fine sand":0.54,"medium sand":0.60,"coarse sand":0.66,
       "very coarse sand":0.70,"granule":0.72,"pebble":0.74,"cobble":0.76}
_LVD=os.path.join(BASE,"age_model_assessment/laventa")
LV_CS=pd.read_csv(os.path.join(_LVD,"laventa_cs_beds.csv")).fillna({"codes":"","grading":""})
LV_CS=LV_CS[LV_CS.base_m<LV_TOP].copy(); LV_CS.loc[LV_CS.top_m>LV_TOP,"top_m"]=LV_TOP; LV_CS.loc[LV_CS.base_m<0,"base_m"]=0
# for display, fold beds <1 m into the thicker neighbour (they only render as outline slivers at this scale)
LV_CS=LV_CS.sort_values("base_m").reset_index(drop=True)
while ((LV_CS.top_m-LV_CS.base_m)<1.0).any():
    th=LV_CS.top_m-LV_CS.base_m; i=int(th.idxmin())
    j=i-1 if (i==len(LV_CS)-1 or (i>0 and th[i-1]>=th[i+1])) else i+1
    LV_CS.loc[j,"base_m"]=min(LV_CS.base_m[i],LV_CS.base_m[j]); LV_CS.loc[j,"top_m"]=max(LV_CS.top_m[i],LV_CS.top_m[j])
    for k in ("codes",): LV_CS.loc[j,k]=(str(LV_CS.loc[j,k])+" "+str(LV_CS.loc[i,k])).strip()
    LV_CS=LV_CS.drop(i).reset_index(drop=True)
LV_CS["facies"]=LV_CS.litho+"|"+LV_CS.grain.map(LV_GW).map("{:.2f}".format)
_fk=sorted(set(LV_CS.facies)); _fc=LV_FAC.set_index("label")
LV_STYLE=pystrat.Style(_fk,[_fc.color[k.split("|")[0]] for k in _fk],[float(k.split("|")[1]) for k in _fk],
                       swatch_values=[int(_fc.swatch[k.split("|")[0]]) for k in _fk])
LV_BEDS=LV_CS[["base_m","top_m","facies"]].reset_index(drop=True)
_P=pd.read_csv(os.path.join(_LVD,"laventa_cs_polarity.csv"))
LV_POL=[(r.base_m,min(r.top_m,LV_TOP),r.polarity,r.chron) for _,r in _P.iterrows() if r.base_m<LV_TOP]
LV_C5AAr=tuple(_P.loc[_P.chron=="C5AAr",["base_m","top_m"]].values[0])          # 153.5-176.6 m
LV_MONKEY=559.2                           # Los Mangos Monkey Unit locality, StL 9 (suppl. s10)
LV_MB=tuple(LV_CS.loc[LV_CS.codes.str.contains("_Red")&(LV_CS.base_m>550)&(LV_CS.base_m<570),["base_m","top_m"]].values[0])
LV_FISH=tuple(LV_CS.loc[LV_CS.fish==1,["base_m","top_m"]].values[0])
# La Venta U-Pb age model (Mora-Rojas 2023, table 4): name, strat_m, median, CI_young, CI_old
LV_UPB=[("44017",63.5,16.00,15.94,16.00),("Z1",120,15.20,13.18,16.00),("UR-01",154,13.26,12.72,13.77),
        ("UR-02",156.6,13.13,12.70,13.46),("JG-R89-1",205,13.13,12.70,13.46),("44011",218.8,13.13,12.69,13.47),
        ("JG-R90-3",242,13.13,12.68,13.45),("JG-R89-3",249,13.11,12.68,13.46),("Z2",303,12.40,11.81,12.92),
        ("JG-R90-1",329,12.40,11.82,12.91)]
def render_litho_lv():
    th=(LV_BEDS.top_m-LV_BEDS.base_m).tolist(); tot=LV_BEDS.top_m.max()
    sec=pystrat.Section(th,LV_BEDS.facies.tolist())
    f0,a0=plt.subplots(figsize=(2.0,9.0)); f0.subplots_adjust(0,0,1,1)
    sec.plot(LV_STYLE,ax=a0,xticks=False); a0.set_xlim(0,LWIDTH); a0.set_ylim(0,tot); a0.axis("off")
    p=os.path.join(LIT,"litho_LV.png"); f0.savefig(p,dpi=150,transparent=True); plt.close(f0); return p
LV_IMG=render_litho_lv()
# ---- Tamshiyacu sedimentology from Gingras et al. (2002, J. Sed. Res.; parasequences TmA, TmB) ----
# Gingras metres -> this section via the shared flooding-surface ties: lignite top / TmA base (0.82 -> 0.26),
# top of the F2B sand + ironstone-nodule horizon (2.15 -> 2.06), TmB lignite base (3.08 -> 2.66) and
# top (3.30 -> 2.76). Above 3.30 m the offset is held constant.
_GT=np.array([[0.82,0.26],[2.15,2.06],[3.08,2.66],[3.30,2.76]])
def g2t(g):
    g=np.asarray(g,float); return np.where(g<=_GT[-1,0],np.interp(g,_GT[:,0],_GT[:,1]),g-_GT[-1,0]+_GT[-1,1])
C_IRON="#a0522d"
# Santa Carina = Santa Julia locality (field sheet "SantaJulia_3"); Gingras' Santa Julia II log (SJA-SJC).
# Tie: SJA-capping lignite top / SJB basal TSE (3.25 m) -> 0.65 m here (top of the lower lignite, burrowed
# "organic clay with burrows" above it); SJC basal TSE (7.0 m) -> 4.90 m (base of the upper lignite).
def g2sc(g): return 0.65+(np.asarray(g,float)-3.25)*(4.90-0.65)/(7.0-3.25)
MARK={"TAM":dict(pmag=(0.00,2.54),nodata_top=(2.54,3.78),puru=2.71,
                 molluscs=[0.20,0.33,2.06,float(g2t(3.36))],   # + TmB basal shell bed (F3B, 3.30-3.42 m)
                 concretions=[],
                 iron_nodules=[float(g2t(2.15))],               # laterally continuous ironstone-nodule horizon
                 cu=[tuple(g2t([0.85,1.57])),                    # TmA: F1A mud -> F2C sand/mud interbeds
                     tuple(g2t([3.42,3.95]))],                   # TmB: laminated clay -> starved ripples thicken up
                 fu=[tuple(g2t([1.65,2.10]))],                   # post-storm, bioturbated silty sand (F2B)
                 firmground=[float(g2t(0.82)),float(g2t(3.08))], # Glossifungites (Thalassinoides) TSEs
                 hcs=[float(g2t(1.61))], ripples=[float(g2t(1.40)),float(g2t(3.60))],
                 bioturbation=[float(g2t(1.88)),float(g2t(2.60)),float(g2t(3.80))],
                 roots=[float(g2t(2.95))]),                     # roots from the TmB lignite into F1B
      "SC":dict(pmag=(1.60,5.00),
                molluscs=[0.30,0.60,1.90,2.30],concretions=[1.80,2.00,2.60,4.20],
                cu=[(0.0,0.45),                                  # top of SJA (CU sand below the lignite)
                    tuple(g2sc([5.80,6.95]))],                   # SJB top: silty -> locally sandy mud
                fu=[],
                firmground=[0.65,4.98],                          # SJB and SJC basal Glossifungites TSEs
                roots=[0.38],                                    # rhizoliths below the SJA-capping lignite
                ripples=[float(g2sc(3.45)),float(g2sc(4.60))],   # starved oscillation ripples + shell lenses
                hcs=[float(g2sc(4.05))],                         # low-angle cross-laminated bed (Arenicolites)
                bioturbation=[0.18,float(g2sc(3.9)),float(g2sc(5.2)),float(g2sc(6.3))])}
# parasequences (Gingras et al. 2002) in each section's own metres: (base, top, label)
PSEQ={"TAM":[(0.26,2.66,"TmA"),(2.66,3.78,"TmB")],
      "SC":[(0.0,0.65,"SJA"),(0.65,4.90,"SJB"),(4.90,5.00,"")]}
# grain-size classes mapped to the pystrat facies widths (rightward = coarser)
GRAINAX=[(0.40,"clay"),(0.52,"silt"),(0.66,"sand")]
C_SHELL="#7d6608"; C_CONC="#5d6d7e"; C_FE="#2a9d8f"; C_S="#e76f51"

# ---------------- polarity + directions ----------------
PZ=pd.read_csv(os.path.join(PS,"polarity_zones.csv")); POL={s:d for s,d in PZ.groupby("site")}
S=pd.read_csv(os.path.join(PS,"samples_corrected.csv"))
ORDER={"TAM":["TAM_3A_4_5E_5D_5C","TAM_5B_5A_6_7"],
       "SC":["SC_1C_1B_1A_2_3C","SC_3B_3A_4D_4C_4B_4A","SC_5_6_7C_7B_7A","SC_8"]}
def composite(site):
    d=S[S.group.isin(ORDER[site])].copy(); base=0.0; out=[]; meta={}
    for g in ORDER[site]:
        dg=d[d.group==g].copy(); th=dg.sed_depth_cm.max()
        dg["comp_cm"]=base+dg.sed_depth_cm; meta[g]=(base,th); base+=th; out.append(dg)
    return pd.concat(out), meta
DIR={}; META={}
for s in ("TAM","SC"): DIR[s],META[s]=composite(s)
def toH(site,cm):
    z=POL[site]; lo=z.top_cm.min(); hi=z.bottom_cm.max(); a,b=MARK[site]["pmag"]
    return b-(cm-lo)/(hi-lo)*(b-a)
def rev_base(site):
    z=POL[site]; return toH(site,z[z.polarity=="R"].bottom_cm.max())
def rev_top(site):
    z=POL[site]; return toH(site,z[z.polarity=="R"].top_cm.min())
datum=rev_base("TAM"); offSC=datum-rev_base("SC"); OFF={"TAM":0.0,"SC":offSC}

# ---------------- XRF -> field metres (per-group fractional registration) ----------------
# ---- XRF cleaning + log-ratio normalisation (Weltje & Tjallingii 2008) ----
# QC: keep surface-valid points with adequate total counts and good spectral fit;
# 3-pt median despike; normalise Fe,S to coherent scatter (removes matrix/count-rate effects).
CPS_MIN=80000; MSE_MAX=3.5; GAPCM=4.0
X=pd.read_csv(XRF)
X=X[(X.valid==1)&(X.cps>=CPS_MIN)&(X.MSE<=MSE_MAX)].copy()
# u-channels scanned for XRF but without palaeomagnetic data, placed in field metres: TAM_1_2_3B sits above
# the measured magnetics; the polarity column is extended over it as "no data" so XRF and pmag lanes match.
XTRA={"TAM":{"TAM_1_2_3B":(2.54,3.78)}}
def xrf_site(site):
    frames=[]
    def add(dg,hfun):
        dg=dg.sort_values("pos_cm").copy()
        for c in ("Fe","S","Mo_coh"): dg[c]=dg[c].rolling(3,center=True,min_periods=1).median()
        dg["H"]=hfun(dg.pos_cm.values)
        dg["lnFe"]=np.log((dg.Fe+1)/(dg.Mo_coh+1)); dg["lnS"]=np.log((dg.S+1)/(dg.Mo_coh+1))
        frames.append(dg[["group","pos_cm","H","lnFe","lnS"]])
    for g in ORDER[site]:
        dg=X[X.group==g]
        if not len(dg): continue
        b,th=META[site][g]; gmax=dg.pos_cm.max()
        add(dg, lambda p,b=b,th=th,gmax=gmax: toH(site,b+(p/gmax)*th)+OFF[site])
    for g,(fb,ft) in XTRA.get(site,{}).items():
        dg=X[X.group==g]
        if not len(dg): continue
        gmax=dg.pos_cm.max()
        add(dg, lambda p,fb=fb,ft=ft,gmax=gmax: (ft-(p/gmax)*(ft-fb))+OFF[site])
    return pd.concat(frames).sort_values("H",ascending=False).reset_index(drop=True)
XD={s:xrf_site(s) for s in ("TAM","SC")}
def broken(df,col):
    """arrays with NaN inserted at group changes or pos gaps>GAPCM (no bridging)."""
    H=[]; V=[]; pg=None; pp=None
    for _,r in df.iterrows():
        if pg is not None and (r.group!=pg or abs(r.pos_cm-pp)>GAPCM):
            H.append(np.nan); V.append(np.nan)
        H.append(r.H); V.append(r[col]); pg=r.group; pp=r.pos_cm
    return np.array(V),np.array(H)
allFe=np.concatenate([XD[s].lnFe.values for s in XD]); allS=np.concatenate([XD[s].lnS.values for s in XD])
feLo,feHi=np.nanpercentile(allFe,[2,98]); sLo,sHi=np.nanpercentile(allS,[2,98])

# ---------------- lane geometry ----------------
LITHO=(0.0,LWIDTH); MARKL=(0.76,1.26); INC=(1.34,2.48); POLX=(2.58,2.82); FEX=(2.96,3.58); SX=(3.72,4.34)
INC0,INC1=-45,45
def xinc(v): return INC[0]+(v-INC0)/(INC1-INC0)*(INC[1]-INC[0])
LANE=4.72
def fscale(v,lane,hi): return lane[0]+np.clip(v/hi,0,1)*(lane[1]-lane[0])
SYM_MS=9                        # accessory-symbol size (points)
FU_X=(0.775,0.865)              # grading-triangle strip inside the feature lane
SYM_XS=(0.97,1.15)              # point-symbol columns inside the feature lane (2nd = overflow)
def data_per_pt(ax):
    """vertical data units per typographic point for an axes whose ylim is already set."""
    h_in=ax.get_position().height*ax.figure.get_size_inches()[1]
    y0,y1=ax.get_ylim(); return abs(y1-y0)/(h_in*72)
def draw_features(ax, x0, items, grading=(), fu_x=FU_X, sym_xs=SYM_XS, size=SYM_MS):
    """items: [(height, symbol[, facecolour])]; grading: [(base, top, 'FU'|'CU')].
    Symbols closer than one symbol-height spill into the next column instead of overlapping."""
    for b,t,k in grading:
        ss.interval(ax,x0+fu_x[0],x0+fu_x[1],b,t,"fining_upward" if k=="FU" else "coarsening_upward")
    gap=size*1.05*data_per_pt(ax); last=[-np.inf]*len(sym_xs)
    for it in sorted(items,key=lambda r:r[0]):
        h,name=it[0],it[1]; fc=it[2] if len(it)>2 else None
        j=next((j for j,l in enumerate(last) if h-l>=gap),int(np.argmin(last)))
        last[j]=h; ss.place(ax,x0+sym_xs[j],h,name,size=size,fc=fc)
def grain_axis(ax,xs,ybase):
    y=ybase-0.12
    ax.plot([xs+LITHO[0],xs+LITHO[0]+0.66],[y,y],color="k",lw=0.6,zorder=7)
    for w,lab in GRAINAX:
        ax.plot([xs+w,xs+w],[y,y-0.04],color="k",lw=0.6,zorder=7)
        ax.text(xs+w,y-0.055,lab,ha="right",va="top",fontsize=5.0,rotation=40,rotation_mode="anchor")
    ax.text(xs+0.20,y-0.30,"grain size",ha="center",va="top",fontsize=5.4,style="italic")

def draw_site(ax,site,xs):
    off=OFF[site]; beds=BEDS[site].copy(); mk=MARK[site]
    ymin=beds.base_m.min()+off; ymax=beds.top_m.max()+off
    # lithology (pystrat image)
    img=plt.imread(LIMG[site])
    ax.imshow(img,extent=[xs+LITHO[0],xs+LITHO[0]+LWIDTH,ymin,ymax],origin="upper",
              aspect="auto",zorder=3,interpolation="antialiased")
    # feature lane: grading triangles (left strip) + point symbols (right), all from ./symbols/*.svg
    items=[(h+off,"bivalve") for h in mk.get("molluscs",[])]+[(h+off,"concretion") for h in mk.get("concretions",[])]
    items+=[(h+off,"concretion",C_IRON) for h in mk.get("iron_nodules",[])]
    for key in ("firmground","hcs","ripples","bioturbation","roots"):
        items+=[(h+off,key) for h in mk.get(key,[])]
    grading=[(b+off,t+off,"FU") for b,t in mk.get("fu",[])]+[(b+off,t+off,"CU") for b,t in mk.get("cu",[])]
    draw_features(ax,xs,items,grading)
    # parasequence brackets (left of the column) and their basal surfaces across the column
    for b,t,lab in PSEQ.get(site,[]):
        b0,t0=b+off,t+off
        ax.plot([xs-0.035,xs-0.075,xs-0.075,xs-0.035],[b0+0.015,b0+0.015,t0-0.015,t0-0.015],color="#6d4c2f",lw=0.7,zorder=6)
        if lab: ax.text(xs-0.09,(b0+t0)/2,lab,rotation=90,ha="right",va="center",fontsize=5.6,color="#6d4c2f",weight="bold")
        if b>0: ax.plot([xs+LITHO[0],xs+LITHO[0]+LWIDTH],[b0,b0],color="#6d4c2f",lw=0.9,ls=(0,(3,1.5)),zorder=7)
    # grain-size axis under the column
    grain_axis(ax,xs,ymin)
    # inclination lane
    exp=np.degrees(np.arctan(2*np.tan(np.radians(3.7))))
    for xv,c,lw,st in [(0,"#c0392b",1.0,"-"),(-exp,"0.6",.7,":"),(exp,"0.6",.7,":")]:
        ax.plot([xs+xinc(xv)]*2,[ymin,ymax],color=c,lw=lw,ls=st,zorder=2)
    d=DIR[site].copy(); d["H"]=toH(site,d.comp_cm)+off; d=d.sort_values("H")
    ax.plot(xs+xinc(d.inc.clip(INC0,INC1)),d.H,color="0.6",lw=0.6,zorder=3)
    ax.plot(xs+xinc(d.inc_smooth.clip(INC0,INC1)),d.H,color="k",lw=1.4,zorder=4)
    ax.scatter(xs+xinc(d.inc.clip(INC0,INC1)),d.H,s=5,
               c=np.where(d.inc<0,"#2166ac","#b2182b"),zorder=5,lw=0)
    # polarity
    z=POL[site]
    # polarity column spans only the measured (cored) interval
    zh=[toH(site,r.top_cm)+off for _,r in z.iterrows()]+[toH(site,r.bottom_cm)+off for _,r in z.iterrows()]
    pol_lo,pol_hi=min(zh),max(zh)
    for fb,ft in XTRA.get(site,{}).values():      # XRF-only u-channels -> hatched no-data polarity
        pol_lo,pol_hi=min(pol_lo,fb+off),max(pol_hi,ft+off)
    ax.add_patch(Rectangle((xs+POLX[0],pol_lo),POLX[1]-POLX[0],pol_hi-pol_lo,
                 facecolor="white",ec="0.6",lw=.5,hatch="////",zorder=3))
    for _,r in z.iterrows():
        y0,y1=sorted((toH(site,r.top_cm)+off,toH(site,r.bottom_cm)+off))
        ax.add_patch(Rectangle((xs+POLX[0],y0),POLX[1]-POLX[0],y1-y0,
                     facecolor=("k" if r.polarity=="N" else "w"),ec="k",lw=.4,zorder=4))
    for _,r in beds.iterrows():   # no-data over lignites, within the measured interval only
        if "lignite" in r.facies:
            g0,g1=r.base_m+off,r.top_m+off
            if g1>pol_lo and g0<pol_hi:
                ax.add_patch(Rectangle((xs+POLX[0],max(g0,pol_lo)),POLX[1]-POLX[0],
                             min(g1,pol_hi)-max(g0,pol_lo),facecolor="white",ec="0.6",lw=.4,hatch="////",zorder=6))
    for _,r in z.iterrows():
        if r.polarity=="R":
            yc=(toH(site,r.top_cm)+toH(site,r.bottom_cm))/2+off
            ax.annotate("C5AAr",(xs+POLX[1],yc),xytext=(xs+POLX[1]+0.05,yc),
                        fontsize=7,color="#c0392b",va="center",weight="bold")
    if mk.get("puru") is not None:
        h=mk["puru"]+off; ss.place(ax,xs-0.30,h,"purussaurus",size=24,zorder=8)
        ax.annotate(f"Purussaurus\ninferred {PURU_P[1]:.2f} Ma\n(95%% CI {PURU_P[0]:.2f}–{PURU_P[2]:.2f})".replace("%%","%"),
                    (xs-0.46,h),xytext=(xs-0.62,h),fontsize=6.5,color="#8e44ad",
                    va="center",ha="right",arrowprops=dict(arrowstyle="-",color="#8e44ad",lw=.8))
    # XRF lanes (both sites) — cleaned, log-ratio normalised, own axis
    def sc(v,lane,lo,hi): return lane[0]+np.clip((v-lo)/(hi-lo),0,1)*(lane[1]-lane[0])
    for lane,name,col,lo,hi,lab in [(FEX,"lnFe",C_FE,feLo,feHi,"ln(Fe/coh)"),(SX,"lnS",C_S,sLo,sHi,"ln(S/coh)")]:
        ax.add_patch(Rectangle((xs+lane[0],pol_lo),lane[1]-lane[0],pol_hi-pol_lo,fc="none",ec="0.8",lw=.4,zorder=2))
        v,h=broken(XD[site],name); xv=xs+sc(v,lane,lo,hi)
        ax.plot(xv,h,color=col,lw=.5,zorder=4)
        yax=ymin-0.10
        ax.plot([xs+lane[0],xs+lane[1]],[yax,yax],color="k",lw=0.6,zorder=7)
        for frac,val in [(0.0,lo),(1.0,hi)]:
            xt=xs+lane[0]+frac*(lane[1]-lane[0])
            ax.plot([xt,xt],[yax,yax-0.04],color="k",lw=0.6,zorder=7)
            ax.text(xt,yax-0.075,f"{val:.1f}",ha="center",va="top",fontsize=4.6)
        ax.text(xs+(lane[0]+lane[1])/2,yax-0.26,lab,ha="center",va="top",fontsize=5.2,style="italic")
    # headers
    yh=ymax+0.16
    ax.text(xs+sum(LITHO)/2+0.02,yh,"lith.",ha="center",fontsize=7)
    ax.text(xs+sum(MARKL)/2,yh,"feat.",ha="center",fontsize=6)
    ax.text(xs+sum(INC)/2,yh,"inclination (°)",ha="center",fontsize=7.5)
    ax.text(xs+sum(POLX)/2,yh,"pol.",ha="center",fontsize=7)
    ax.text(xs+sum(FEX)/2,yh,"Fe",ha="center",fontsize=7.5,color="#8a5a00")
    ax.text(xs+sum(SX)/2,yh,"S",ha="center",fontsize=7.5,color="#5b3f5b")
    ax.text(xs+(LITHO[0]+SX[1])/2,ymax+0.42,{"TAM":"Tamshiyacu (TAM)","SC":"Santa Carina (SC)"}[site],
            ha="center",fontsize=11,weight="bold")
    for v in (-40,0,40):
        ax.text(xs+xinc(v),ymin-0.18,f"{v}",ha="center",fontsize=6)

# ============ Purussaurus inferred age (Monte Carlo via La Venta age model) ============
# TAM Purussaurus (2.71 m) sits just above the measured C5AAr reversal top -> bracketed C5AAr–C5AAn.
# Propagate the U-Pb-calibrated La Venta age model (Mora-Rojas Table 4) at the C5AAr-top level (176.6 m),
# intersected with the magnetostratigraphic C5AAr–C5AAn bracket (conservative).
_rng=np.random.default_rng(42)
_LV=np.array([(63.5,16.00,15.94,16.00),(120,15.20,13.18,16.00),(154,13.26,12.72,13.77),
 (156.6,13.13,12.70,13.46),(205,13.13,12.70,13.46),(218.8,13.13,12.69,13.47),
 (242,13.13,12.68,13.45),(249,13.11,12.68,13.46),(303,12.40,11.81,12.92),(329,12.40,11.82,12.91)])
_o=np.argsort(_LV[:,0]); _dd,_mm,_llo,_hhi=_LV[_o,0],_LV[_o,1],_LV[_o,2],_LV[_o,3]; _sg=(_hhi-_llo)/3.92
# correlative La Venta level: the C5AAr top in the composite section (LV_C5AAr[1] = 176.6 m, from the
# s15 polarity column), drawn -20/+15 m to cover the no-polarity gap above the TAM record and cross-basin mapping
_dep=_rng.uniform(LV_C5AAr[1]-20,LV_C5AAr[1]+15,300000); _al=_rng.normal(np.interp(_dep,_dd,_mm),np.interp(_dep,_dd,_sg))
PURU_AGE=_al[(_al>=13.032)&(_al<=13.363)]; PURU_ALL=_al
PURU_P=np.percentile(PURU_AGE,[2.5,50,97.5])   # ~13.04, 13.18, 13.35
def draw_age_pdf(a,compact=False):
    fs=6 if compact else 9
    bins=np.linspace(12.6,13.9,90 if compact else 150)
    a.axvspan(13.183,13.363,color="#c0392b",alpha=0.12); a.axvspan(13.032,13.183,color="#2166ac",alpha=0.08)
    if not compact:
        a.hist(PURU_ALL,bins=bins,density=True,color="#8e44ad",alpha=0.28,lw=0,label="La Venta model only")
    a.hist(PURU_AGE,bins=bins,density=True,color="#6c3483",alpha=0.75,lw=0,
           label=("posterior" if compact else "La Venta model ∩ C5AAr–C5AAn"))
    for x in PURU_P: a.axvline(x,color="#4a235a",ls=":",lw=0.8)
    a.set_xlim(13.9,12.6); a.set_yticks([])
    a.text(13.27,a.get_ylim()[1]*0.86,"C5AAr",ha="center",fontsize=fs-0.5,color="#c0392b")
    a.text(13.11,a.get_ylim()[1]*0.86,"C5AAn",ha="center",fontsize=fs-0.5,color="#2166ac")
    a.tick_params(labelsize=fs-0.5)
    for sp in ("top","left","right"): a.spines[sp].set_visible(False)

# ================= figure =================
fig=plt.figure(figsize=(15,11.5))
gs=fig.add_gridspec(1,4,width_ratios=[3.15,1.05,0.42,0.62],wspace=0.20)
ax=fig.add_subplot(gs[0]); axlv=fig.add_subplot(gs[1]); axg=fig.add_subplot(gs[2]); axu=fig.add_subplot(gs[3])
for _a in (ax,axlv,axg,axu): _a.set_facecolor("none")   # lets the figure-level tie ribbons show through
ymin=min(0,SC_BEDS.base_m.min()+offSC)-0.72; ymax=max(TAM_BEDS.top_m.max(),SC_BEDS.top_m.max()+offSC)+0.78
ax.set_ylim(ymin,ymax); ax.set_xlim(-1.2,LANE+SX[1]+0.35)   # limits first: symbol spacing needs them
draw_site(ax,"TAM",0.0); draw_site(ax,"SC",LANE)
ax.set_ylabel("Height in reversal-aligned frame (m)  ·  0 = TAM base; negative = SC below the C5AAr datum",
              fontsize=9,labelpad=26)
ax.set_yticks(np.arange(-2.0,ymax,0.5)); ax.set_xticks([]); ax.tick_params(axis="y",pad=4)
for sp in ("top","right","bottom"): ax.spines[sp].set_visible(False)
ax.set_title("Peru — Pebas Fm (TAM, SC)",fontsize=11,weight="bold")

# inferred-age posterior inset in the blank space beneath TAM
axin=ax.inset_axes([0.15,-2.45,3.9,1.15],transform=ax.transData)
draw_age_pdf(axin,compact=True)
axin.set_xlabel("inferred Purussaurus age (Ma)",fontsize=6)
axin.set_title("Age via La Venta U–Pb model + C5AAr tie",fontsize=6.2,weight="bold")
axin.text(0.5,-0.30,f"{PURU_P[1]:.2f} Ma (95%% CI {PURU_P[0]:.2f}–{PURU_P[2]:.2f})".replace("%%","%"),
          transform=axin.transAxes,ha="center",fontsize=6,color="#6c3483",weight="bold")

C5AAr=(13.183,13.363)
NORM={"C5An.1n":(12.049,12.174),"C5An.2n":(12.272,12.474),"C5Ar.1n":(12.735,12.770),
 "C5Ar.2n":(12.829,12.887),"C5AAn":(13.032,13.183),"C5ABn":(13.363,13.608),"C5ACn":(13.739,14.070)}
scr_top=rev_top("SC")+offSC; scr_bot=rev_base("SC")+offSC
# Right-hand block (La Venta metres | GPTS | U-Pb, the last two in Ma) has its own vertical scales;
# columns are tied by ribbons/lines instead of being forced to share one band height.
AGE_YOUNG,AGE_OLD=11.95,13.75
C_TIE="#c0392b"; C_LVP="#8B0000"

# ---------------- La Venta composite section (CS metres, continuous 0-625 m) ----------------
LV_LITHO=(0.0,LWIDTH); LV_FEATL=(0.76,1.22); LV_POLX=(1.28,1.48)
axlv.set_ylim(-18,LV_TOP+5); axlv.set_xlim(-0.95,2.02)
axlv.imshow(plt.imread(LV_IMG),extent=[LV_LITHO[0],LV_LITHO[0]+LWIDTH,0,LV_BEDS.top_m.max()],
            origin="upper",aspect="auto",zorder=3,interpolation="antialiased")
for b,t,pol,lab in LV_POL:
    axlv.add_patch(Rectangle((LV_POLX[0],b),LV_POLX[1]-LV_POLX[0],t-b,fc=("k" if pol=="N" else "w"),ec="k",lw=.4,zorder=4))
    if t-b>=12: axlv.text(LV_POLX[1]+0.04,(b+t)/2 if lab!="C5An.1n" else b+12,lab,va="center",fontsize=5.4,
                          color=(C_TIE if lab=="C5AAr" else "0.25"),weight=("bold" if lab=="C5AAr" else "normal"))
axlv.text(LV_POLX[1]+0.04,277,"C5Ar.3r–1r",va="center",fontsize=5.0,color="0.35")
# marker horizons: Monkey Beds (StL 9) and Fish Bed (StL 11), bracketed left of the column
for (b,t),lab in [(LV_MB,"Monkey Beds (StL 9)\nPurussaurus neivensis\n13 m above C5An.1n base\n≤12.17 Ma"),(LV_FISH,"Fish Bed\n(StL 11)")]:
    axlv.plot([-0.06,-0.10,-0.10,-0.06],[b,b,t,t],color=C_LVP,lw=0.8,zorder=6)
    axlv.text(-0.13,(b+t)/2,lab,ha="right",va="center",fontsize=5.2,color=C_LVP,weight="bold",linespacing=0.95)
    axlv.add_patch(Rectangle((LV_LITHO[0],b),LWIDTH,t-b,fc=C_LVP,alpha=0.22,ec="none",zorder=5))
    axlv.add_patch(Rectangle((LV_LITHO[0],b),LWIDTH,t-b,fc="none",ec=C_LVP,lw=0.9,zorder=5))
ss.place(axlv,1.02,LV_MONKEY,"purussaurus",size=19,zorder=8,fc=C_LVP)
ss.place(axlv,1.02,sum(LV_FISH)/2,"fish",size=8.5,zorder=8)
for _x,_t,_f in [(LV_LITHO[0]+LWIDTH/2,"lith.",7),(sum(LV_FEATL)/2,"feat.",6),(sum(LV_POLX)/2,"pol.",7)]:
    axlv.text(_x,1.004,_t,ha="center",va="bottom",fontsize=_f,transform=axlv.get_xaxis_transform(),clip_on=False)
# grain-size axis under the La Venta lithology column
_gy=-6
axlv.plot([LV_LITHO[0],LV_LITHO[0]+0.76],[_gy,_gy],color="k",lw=0.6,clip_on=False,zorder=7)
for w,lab in [(0.40,"mud"),(0.54,"f sand"),(0.66,"c sand"),(0.74,"gvl")]:
    axlv.plot([w,w],[_gy,_gy-2.5],color="k",lw=0.6,clip_on=False,zorder=7)
    axlv.text(w,_gy-4,lab,ha="right",va="top",fontsize=4.6,rotation=40,rotation_mode="anchor",clip_on=False)
axlv.text(0.30,_gy-17,"grain size",ha="center",va="top",fontsize=5.0,style="italic",clip_on=False)
axlv.set_title("La Venta, Honda Gp\n(composite section)",fontsize=9.5,weight="bold",pad=12)
# ---- accessory features from the CS lithofacies codes (Mora-Rojas et al. 2023, suppl. s6) ----
_lv_items=[]
for _,r in LV_CS.iterrows():
    mid=(r.base_m+r.top_m)/2
    if abs(mid-LV_MONKEY)<9 or abs(mid-sum(LV_FISH)/2)<5: continue      # keep the marker horizons clear
    for k in ("trough_xbed","planar_xbed","epsilon_xstrat","pebbles","concretion","ripples"):
        if r[k]: _lv_items.append((mid,k))
    if r.mottling and r.top_m-r.base_m>=3: _lv_items.append((mid,"mottling"))
    elif r.planar_lamination and r.top_m-r.base_m>=4: _lv_items.append((mid,"planar_lamination"))
_lv_grad=[]   # merge touching beds with the same grading sense; keep >=3 m packages
for _,r in LV_CS[LV_CS.grading.isin(["FU","CU"])].sort_values("base_m").iterrows():
    if _lv_grad and _lv_grad[-1][2]==r.grading and r.base_m-_lv_grad[-1][1]<0.5: _lv_grad[-1][1]=r.top_m
    else: _lv_grad.append([r.base_m,r.top_m,r.grading])
_lv_grad=[tuple(g) for g in _lv_grad if g[1]-g[0]>=3]
draw_features(axlv,0.0,_lv_items,_lv_grad,fu_x=(0.775,0.865),sym_xs=(0.98,1.14),size=7)
axlv.set_yticks(np.arange(0,LV_TOP+1,100)); axlv.tick_params(labelsize=6)
axlv.set_ylabel("La Venta composite section (m)",fontsize=8,labelpad=1); axlv.yaxis.set_label_position("left"); axlv.yaxis.tick_left()
axlv.spines["left"].set_position(("data",-0.95))
axlv.set_xticks([])
for sp in ("top","right","bottom"): axlv.spines[sp].set_visible(False)
axlv.axhspan(LV_C5AAr[0],LV_C5AAr[1],xmin=0.30,xmax=0.83,color=C_TIE,alpha=0.10,zorder=0.5)

# ---------------- GPTS (Ma) ----------------
MZ8=(12.8,13.3)
gx0,gx1=0.05,0.45
axg.set_ylim(AGE_OLD,AGE_YOUNG); axg.set_xlim(0,1.25)
axg.add_patch(Rectangle((gx0,AGE_YOUNG),gx1-gx0,AGE_OLD-AGE_YOUNG,fc="white",ec="k",lw=.8))
for nm,(t,b) in NORM.items():
    if b>AGE_YOUNG and t<AGE_OLD:
        axg.add_patch(Rectangle((gx0,max(t,AGE_YOUNG)),gx1-gx0,min(b,AGE_OLD)-max(t,AGE_YOUNG),fc="k",zorder=2))
        if AGE_YOUNG<=(t+b)/2<=AGE_OLD and b-t>0.08: axg.text(gx1+0.04,(t+b)/2,nm,va="center",fontsize=5.4)
axg.text(gx1+0.04,12.80,"C5Ar",va="center",fontsize=5.0,color="0.35")
axg.axhspan(C5AAr[0],C5AAr[1],xmax=0.36,color=C_TIE,alpha=0.12,zorder=1)
axg.text(gx1+0.04,sum(C5AAr)/2,"C5AAr",va="center",fontsize=5.6,color=C_TIE,weight="bold")
axg.axhspan(MZ8[0],MZ8[1],xmin=0.90,xmax=0.97,color="#8a7d1f",alpha=0.35,zorder=1)
axg.text(1.24,sum(MZ8)/2,"Mollusc Zone 8 (~13 Ma)",rotation=90,va="center",ha="right",fontsize=5.0,color="#6b6017")
# --- Crassoretitriletes pollen zone (Rasanen et al. 1998; Gingras et al. 2002; Amazon calibration Jaramillo et al. 2025) ---
POLLEN=(12.4,14.1)  # young top = FAD Grimsdalea magnaclavata; old base = FAD Crassoretitriletes (off-axis)
_py0=max(POLLEN[0],AGE_YOUNG); _py1=min(POLLEN[1],AGE_OLD)
axg.axhspan(_py0,_py1,xmin=0.80,xmax=0.87,color="#2e8b57",alpha=0.28,zorder=1)
axg.plot([0.80*1.25,0.87*1.25],[POLLEN[0],POLLEN[0]],color="#1e6b40",lw=0.9,zorder=2)  # 12.4 Ma Grimsdalea-FAD top
axg.text(0.985,(_py0+_py1)/2,"palyn. Mid-Miocene [14.1\u201312.4 Ma pending]",rotation=90,va="center",ha="right",fontsize=4.6,color="#1e6b40")
axg.set_xticks([]); axg.set_yticks([])
for sp in ("top","left","right","bottom"): axg.spines[sp].set_visible(False)
axg.set_title("GPTS\n(GTS2020)",fontsize=9.5,weight="bold",pad=12)

# ---------------- U-Pb age-model column (La Venta detrital zircons), shares the Ma scale ----------------
inwin=sorted([u for u in LV_UPB if AGE_YOUNG<=u[2]<=AGE_OLD], key=lambda u:u[2])  # youngest->oldest, left->right
xs_u=np.linspace(0.10,0.90,len(inwin))
axu.set_ylim(AGE_OLD,AGE_YOUNG); axu.set_xlim(0,1)
axu.axhspan(C5AAr[0],C5AAr[1],color=C_TIE,alpha=0.10,zorder=1)
for (nm,sm,med,ylo,yhi),xu in zip(inwin,xs_u):
    big=(nm=="44011")
    axu.plot([xu,xu],[max(ylo,AGE_YOUNG),min(yhi,AGE_OLD)],color="#8e44ad",lw=(1.8 if big else 1.0),zorder=3)
    axu.plot(xu,med,"o",color="#6c3483",ms=(5 if big else 3.2),zorder=4)
    lab=(f"{nm} ({sm:g} m): {med:.2f} ({ylo:.2f}–{yhi:.2f}); YSP 13.9±0.3" if big
         else f"{nm} ({sm:g} m): {med:.2f} ({ylo:.2f}–{yhi:.2f})")
    _top=ylo>AGE_YOUNG+0.05                  # label at the young end unless the bar runs off the top
    axu.text(xu,(ylo-0.02) if _top else (min(yhi,AGE_OLD)+0.02),lab,rotation=90,va=("bottom" if _top else "top"),ha="center",
             fontsize=(4.8 if big else 4.3),color=("#6c3483" if big else "0.4"),
             weight=("bold" if big else "normal"),clip_on=False)
axu.set_xticks([])
axu.set_ylabel("Age (Ma)",fontsize=10,labelpad=10); axu.yaxis.set_label_position("right")
axu.yaxis.tick_right(); axu.tick_params(axis="y",pad=3)
for sp in ("top","left","bottom"): axu.spines[sp].set_visible(False)
axu.set_title("U–Pb age model\n(La Venta zircons)",fontsize=8.5,weight="bold",pad=12)

# ---------------- ties: Peru C5AAr -> La Venta C5AAr -> GPTS; La Venta reversals -> GPTS ----------------
from matplotlib.patches import ConnectionPatch
# C5AAr band in the Peru panel: each site shaded over its own full reversed interval
# (TAM: whole measured R zone up to the top of the magnetic record; SC: its R zone), stepping at the site gap
_x0,_x1=ax.get_xlim(); _fx=lambda x:(x-_x0)/(_x1-_x0); _xsplit=LANE-0.28
tam_bot=rev_base("TAM"); tam_top=rev_top("TAM")
ax.axhspan(tam_bot,tam_top,xmin=0,xmax=_fx(_xsplit),color=C_TIE,alpha=0.09,zorder=0.5)
ax.axhspan(scr_bot,scr_top,xmin=_fx(_xsplit),xmax=1,color=C_TIE,alpha=0.09,zorder=0.5)
fig.canvas.draw()                                   # freeze axes positions before mapping to figure coords
def F(a,x,y): return fig.transFigure.inverted().transform(a.transData.transform((x,y)))
_xr=ax.get_xlim()[1]
rib1=[F(ax,_xr,scr_bot),F(ax,_xr,scr_top),F(axlv,LV_LITHO[0]-0.02,LV_C5AAr[1]),F(axlv,LV_LITHO[0]-0.02,LV_C5AAr[0])]
rib2=[F(axlv,LV_POLX[1],LV_C5AAr[0]),F(axlv,LV_POLX[1],LV_C5AAr[1]),F(axg,gx0,C5AAr[0]),F(axg,gx0,C5AAr[1])]
for rb in (rib1,rib2):
    fig.add_artist(MplPolygon(rb,closed=True,transform=fig.transFigure,fc=C_TIE,alpha=0.13,ec=C_TIE,lw=0.4,zorder=-1))
# every La Venta reversal to its GPTS boundary (Mora-Rojas 2023 reinterpretation of Flynn et al. 1997)
_GB={"C5ABn":13.363,"C5AAr":13.183,"C5AAn":13.032,"C5Ar.3r":12.887,"C5Ar.2n":12.829,"C5Ar.2r":12.770,
     "C5Ar.1n":12.735,"C5Ar.1r":12.474,"C5An.2n":12.272,"C5An.1r":12.174}   # age of each zone's TOP
for b,t,pol,lab in LV_POL[:-1]:
    if lab in _GB:
        hot=(lab=="C5An.1r")
        fig.add_artist(ConnectionPatch(xyA=(LV_POLX[1],t),coordsA=axlv.transData,xyB=(gx0,_GB[lab]),coordsB=axg.transData,
                       color=(C_LVP if hot else "0.55"),lw=(0.9 if hot else 0.45),ls=("-" if hot else (0,(2,2))),zorder=-0.5))

# ---------------- legend ----------------
_used=set(TAM_BEDS.facies)|set(SC_BEDS.facies)   # only facies actually drawn (shell beds -> symbol)
leg=[Patch(fc=r.color,ec="k",label=leglabel(r)) for _,r in FAC.iterrows() if r.label in _used]
leg+=[Line2D([0],[0],color="0.6",lw=1,label="raw ChRM inc"),
      Line2D([0],[0],color="k",lw=1.6,label="7-cm median inc"),
      Line2D([0],[0],marker="o",color="w",mfc="#2166ac",ms=6,label="normal (inc<0)"),
      Line2D([0],[0],marker="o",color="w",mfc="#b2182b",ms=6,label="reversed (inc>0)"),
      Patch(fc="k",ec="k",label="N polarity"),Patch(fc="w",ec="k",label="R polarity"),
      Patch(fc="white",ec="0.6",hatch="////",label="no data"),
      ss.handle("purussaurus","Purussaurus (Pebas, TAM)",size=20),
      ss.handle("purussaurus","Purussaurus (La Venta)",size=20,fc="#8B0000"),
      ss.handle("bivalve","mollusc / shell bed (F3B)"),
      ss.handle("concretion","carbonate concretion / nodule"),
      ss.handle("concretion","ironstone-nodule horizon",fc=C_IRON),
      ss.handle("fining_upward","fining-upward (normal grading)"),
      ss.handle("coarsening_upward","coarsening-upward (inverse grading)"),
      ss.handle("firmground","Glossifungites firmground (TSE)"),
      Line2D([0],[0],color="#6d4c2f",lw=0.9,ls=(0,(3,1.5)),label="parasequence base (Gingras et al. 2002)"),
      ss.handle("bioturbation","bioturbation (Planolites, Teichichnus)"),
      ss.handle("roots","root traces"),
      ss.handle("hcs","hummocky cross-stratification"),
      ss.handle("ripples","ripple lamination"),
      ss.handle("planar_lamination","planar-parallel lamination"),
      ss.handle("trough_xbed","trough cross-bedding"),
      ss.handle("planar_xbed","planar / tabular cross-bedding"),
      ss.handle("epsilon_xstrat","epsilon cross-strat. (point bar)"),
      ss.handle("pebbles","conglomeratic lenses"),
      ss.handle("mottling","pedogenic mottling (palaeosol)"),
      Line2D([0],[0],color=C_FE,lw=1.2,label="ln(Fe/coh)"),Line2D([0],[0],color=C_S,lw=1.2,label="ln(S/coh)"),
      Line2D([0],[0],color="#8e44ad",lw=1,marker="o",ms=4,label="U–Pb age model (95% CI)"),
      *[Patch(fc=r.color,ec="k",label="La Venta "+r.label.replace("mudstone_red","red mudstone")) for _,r in LV_FAC.iterrows()
        if r.label in set(LV_CS.litho)],
      ss.handle("fish","fish remains (Fish Bed)"),
      Patch(fc="#8a7d1f",alpha=.30,label="Mollusc Zone 8 (~13 Ma)"),
      Patch(fc="#2e8b57",alpha=.28,label="palynology (Mid-Miocene; 14.1\u201312.4 Ma pending)"),
      Patch(fc="none",ec="#8B0000",lw=1,label="La Venta marker horizon"),
      Patch(fc="#c0392b",alpha=.18,label="C5AAr tie")]
leg=[h if isinstance(h,tuple) and len(h)==2 and isinstance(h[1],str) else (h,h.get_label()) for h in leg]
fig.legend(handles=[h for h,_ in leg],labels=[l for _,l in leg],handler_map=ss.HANDLER_MAP,loc="lower center",ncol=7,fontsize=7.0,frameon=False,bbox_to_anchor=(0.5,-0.15))
fig.suptitle("Figure 2. Pebas magnetostratigraphy (TAM, SC) correlated to the U–Pb-calibrated La Venta composite section and the GPTS",
             fontsize=12.5,weight="bold",y=0.99)
import textwrap
_cap=("Lithology after pystrat (grain-size widths + FGDC swatches); TAM & SC parasequences, grading, ichnology & physical structures after Gingras et al. 2002 (Tamshiyacu TmA/TmB; Santa Carina = their Santa Julia SJA-SJC; tied on the lignite-capped flooding surfaces); La Venta composite section (lithology, lithofacies codes, grading) digitised from Mora-Rojas et al. 2023 suppl. s6, polarity from their s15 reinterpretation of Flynn et al. 1997, reversals tied to Hilgen et al. 2012; Peru columns aligned on the shared C5AAr reversal; "
         "XRF cleaned (cps≥80k, MSE≤3.5, gaps unbridged, 3-pt median) and shown as log-ratios to coherent scatter (Weltje & Tjallingii 2008); "
         "Monkey Beds (StL 9) & Fish Bed (StL 11) from Mora-Rojas et al. 2023 (suppl. s6, s10) with U–Pb Bayesian age model (detrital zircons, 95% CI); "
         "La Venta Purussaurus (P. neivensis, Monkey Beds) after Carrillo et al. 2023; Mollusc Zone 8 after Wesselingh et al. 2006 & Salas-Gismondi et al. 2016.  DRAFT.")
fig.text(0.5,-0.25,"\n".join(textwrap.wrap(_cap,260)),ha="center",va="top",fontsize=6.6,style="italic",color="0.35")
fig.savefig(os.path.join(OUT,"fig2_master.png"),dpi=200,bbox_inches="tight",facecolor="w")
print("wrote fig2_master.png datum=%.2f offSC=%.2f SCrev %.2f–%.2f m feHi=%.0f sHi=%.0f nTAMxrf=%d nSCxrf=%d"
      %(datum,offSC,scr_top,scr_bot,feHi,sHi,len(XD['TAM']),len(XD['SC'])))

# ================= standalone Purussaurus inferred-age figure =================
figA,axA=plt.subplots(figsize=(7,4.2)); draw_age_pdf(axA,compact=False)
axA.set_xlabel("Age (Ma)"); axA.set_ylabel("probability density")
axA.set_title("Inferred age of the Tamshiyacu Purussaurus\nLa Venta U–Pb age model propagated through the C5AAr tie",
              fontsize=10,weight="bold")
axA.legend(fontsize=8,loc="upper left")
axA.text(0.98,0.60,f"posterior\n{PURU_P[1]:.2f} Ma\n95%% CI {PURU_P[0]:.2f}–{PURU_P[2]:.2f}".replace("%%","%"),
         transform=axA.transAxes,ha="right",va="top",fontsize=8.5,color="#4a235a",weight="bold")
figA.tight_layout(); figA.savefig(os.path.join(OUT,"purussaurus_age_pdf.png"),dpi=170,bbox_inches="tight",facecolor="w")
print("wrote purussaurus_age_pdf.png  puru age %.2f (%.2f–%.2f)"%(PURU_P[1],PURU_P[0],PURU_P[2]))
