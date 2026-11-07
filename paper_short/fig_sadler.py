import os
#!/usr/bin/env python3
"""Sadler plot: accumulation rate vs measurement timespan (log-log).
Amazonian/Andean records + condensed-setting references, placing the Iquitos Pebas
magnetostratigraphic rate in context. Each point is labelled with its depositional
SETTING (line 1) and the source study 'Author et al. date' (line 2).

Sadler effect (Sadler 1981; Schumer & Jerolmack 2009): apparent rate falls as the
averaging timespan grows, because longer records contain more hiatuses.

Importable: draw_sadler(ax) renders the whole plot onto a supplied Axes so the
same code backs both the standalone figure and panel (c) of the combined Fig 4."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=os.environ.get("PEBAS_OUT","/mnt/user-data/outputs")

C_FLU="#17a2b8"; C_AND="#7e3f9e"; C_BAS="#b8860b"; C_MAR="#2166ac"; C_CON="#555555"; C_THIS="#c0392b"; C_WET="#2e8b57"
MARK={"modern wetland analogue":"D"}   # diamond for wetland analogues

# label = "<setting>\n<Author et al. date>" ; rate_lo, rate_hi (m/Myr); tspan_lo, tspan_hi (yr); colour; group; dy
D=[
 # ---- modern analogues that fall ON the trend (validate the framework) ----
 ("Fluvial overbank (várzea)\nAalto et al. 2003",500,5000,3,15,C_FLU,"Amazonian fluvial / fan",+1),
 ("Freshwater peat mire\nGlaser et al. 2012",200,800,1.5e3,5e3,C_WET,"modern wetland analogue",+1),
 # ---- Amazon deep-sea fan (same system, two timescales) ----
 ("Deep-sea fan (lowstand)\nFlood et al. 1995",5000,15000,1e4,1e5,C_FLU,"Amazonian fluvial / fan",+1),
 ("Deep-sea fan (Neogene)\nFigueiredo et al. 2009",300,1000,7e6,1.1e7,C_FLU,"Amazonian fluvial / fan",+1),
 # ---- Andean lacustrine (long, well-dated drill cores) ----
 ("Altiplano lacustrine\nFritz et al. 2007",370,370,3e5,4.5e5,C_AND,"Andean lacustrine",-1),
 ("Andean intermontane lake\nTorres et al. 2005",200,300,1.8e6,2.5e6,C_AND,"Andean lacustrine",+1),
 # ---- Amazon foreland / basin (Miocene Pebas) ----
 ("Wetland, orbital cyclicity\nHoorn et al. 2022",100,100,3e4,5e4,C_BAS,"Amazon foreland/basin (Miocene)",-1),
 ("Foreland-basin fill\nWesselingh et al. 2006",45,130,6e6,1e7,C_BAS,"Amazon foreland/basin (Miocene)",+1),
 # ---- open-marine references (GPTS/astro-tuned) ----
 ("Pelagic, Nazca Ridge\nMix et al. 2003",7.7,7.7,1.5e6,2.5e6,C_MAR,"open-marine reference",-1),
 ("Pelagic carbonate, Ceará Rise\nBackman & Raffi 1997",16,16,1.5e6,2.5e6,C_MAR,"open-marine reference",+1),
 ("Pelagic, Carnegie Ridge\nMix et al. 2003",10,72,1.5e6,2.5e6,C_MAR,"open-marine reference",+1),
 # ---- modern sediment-starved shelves (living condensed analogues) ----
 ("Drowned, starved shelf\nHanebuth et al. 2000",10,60,6e3,1.5e4,C_CON,"condensed / starved",+1),
 ("Relict/palimpsest shelf\nEmery 1968",3,40,4e3,1.2e4,C_CON,"condensed / starved",-1),
 # ---- ancient condensed anchor ----
 ("Abyssal pelagic clay\nDunlea et al. 2015",1,4,1e6,1e7,C_CON,"condensed / starved",-1),
 # ---- this study ----
 ("Bay-margin wetland (C5AAr)\nThis study",2.4,9,1.5e5,2.1e5,C_THIS,"THIS STUDY",0),
]

# per-label offset overrides (dx,dy pts, ha, va) keyed by the exact label string
OVR={
 "Altiplano lacustrine\nFritz et al. 2007":(0,-13,"center","top"),
 "Andean intermontane lake\nTorres et al. 2005":(10,10,"left","bottom"),
 "Pelagic, Nazca Ridge\nMix et al. 2003":(-9,-9,"right","top"),
 "Pelagic carbonate, Ceará Rise\nBackman & Raffi 1997":(11,2,"left","center"),
 "Pelagic, Carnegie Ridge\nMix et al. 2003":(-9,10,"right","bottom"),
 "Abyssal pelagic clay\nDunlea et al. 2015":(0,-12,"center","top"),
 "Foreland-basin fill\nWesselingh et al. 2006":(10,11,"left","bottom"),
 "Freshwater peat mire\nGlaser et al. 2012":(0,12,"center","bottom"),
 "Drowned, starved shelf\nHanebuth et al. 2000":(-8,12,"right","bottom"),
 "Relict/palimpsest shelf\nEmery 1968":(-8,-11,"right","top"),
}

def draw_sadler(ax, title=True, fs_lab=5.6):
    """Render the Sadler plot onto Axes `ax`."""
    # global Sadler reference band (schematic): R = A * T^-0.5
    T=np.logspace(0,8.2,100); A=2e4
    ax.fill_between(T,0.32*A*T**-0.5,3.2*A*T**-0.5,color="0.82",alpha=0.35,zorder=0,lw=0)
    ax.plot(T,A*T**-0.5,color="0.55",lw=1.2,ls="--",zorder=1)
    ax.text(3.5e6,3.2,"global Sadler trend (schematic;\nSadler 1981; Schumer & Jerolmack 2009)",
            fontsize=7.2,color="0.5",rotation=-19,ha="center",va="center")
    seen=set()
    for lab,rlo,rhi,tlo,thi,c,g,dy in D:
        tc=np.sqrt(tlo*thi); rc=np.sqrt(rlo*rhi); this=g=="THIS STUDY"
        ax.plot([tlo,thi],[rc,rc],color=c,lw=1.1,alpha=0.6,zorder=3)
        ax.plot([tc,tc],[rlo,rhi],color=c,lw=1.1,alpha=0.6,zorder=3)
        mk="*" if this else MARK.get(g,"o")
        ax.scatter([tc],[rc],s=(300 if this else (85 if mk=="D" else 62)),marker=mk,
                   color=c,ec="k",lw=(1.0 if this else 0.5),zorder=5,label=(g if g not in seen else None))
        seen.add(g)
        if this:
            ax.annotate(lab,(tc,rc),fontsize=fs_lab+1.2,color=c,ha="right",va="center",
                        xytext=(-18,-1),textcoords="offset points",weight="bold")
        elif lab in OVR:
            dx,dyy,ha,va=OVR[lab]
            ax.annotate(lab,(tc,rc),fontsize=fs_lab,color=c,ha=ha,va=va,
                        xytext=(dx,dyy),textcoords="offset points")
        else:
            ax.annotate(lab,(tc,rc),fontsize=fs_lab,color=c,ha="center",
                        va=("bottom" if dy>0 else "top"),
                        xytext=(0,(10 if dy>0 else -10)),textcoords="offset points")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(2,1.5e8); ax.set_ylim(0.2,3e4)
    ax.set_xlabel("Measurement timespan (yr)",fontsize=11)
    ax.set_ylabel("Accumulation rate (m / Myr)",fontsize=11)
    if title:
        ax.set_title("Accumulation rate vs measurement timespan — the Iquitos Pebas margin "
                     "sits with the condensed end-members",fontsize=10.2,weight="bold")
    for x,lab in [(1e2,"century"),(1e4,"10 kyr"),(1e6,"1 Myr"),(1e8,"100 Myr")]:
        ax.text(x,2.2e4,lab,fontsize=6.3,color="0.55",ha="center",va="top")
    ax.legend(fontsize=7.6,loc="upper right",frameon=True,framealpha=0.92,
              title="setting",title_fontsize=8.0)
    ax.grid(True,which="major",color="0.92",lw=0.5,zorder=0)
    ax.text(0.014,0.03,
        "Modern healthy systems (Amazon floodplain, Everglades) fall on the global trend; sediment-starved shelves (Sunda, US Atlantic) fall below it.\n"
        "For its timescale the Iquitos Pebas margin plots with the starved shelves and abyssal pelagic clay — genuinely condensed, not just old.",
        transform=ax.transAxes,fontsize=6.8,color="0.35",va="bottom")

if __name__=="__main__":
    fig,ax=plt.subplots(figsize=(12.4,7.8))
    draw_sadler(ax)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_sadler.png",dpi=190,bbox_inches="tight",facecolor="w")
    fig.savefig(f"{OUT}/fig_sadler.pdf",bbox_inches="tight")
    print("wrote fig_sadler.png with %d points"%len(D))
