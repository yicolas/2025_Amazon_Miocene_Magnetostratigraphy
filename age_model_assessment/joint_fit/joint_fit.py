#!/usr/bin/env python3
"""Joint TAM+SC fit to candidate GPTS chrons.

Shared datum: the N->R boundary going up-section (TAM 50.5 cm; SC 149.5 cm)
= the BASE of the candidate reversed chron. Per-site free rates.
With the datum and r_s free, SC's two dated boundaries (R base & R top) fix
r_s = R_thickness / chron_duration exactly; the information is in the
INEQUALITIES: each flanking normal interval must fit inside its adjacent
normal chron, and TAM's truncated R (>=50 cm) must fit inside the chron.
We report, per candidate: forced rates, every overshoot (kyr), the rate
variability (fold-change between intervals) the candidate demands, and the
pollen/fossil penalties. ATNTS2012 ages."""
import numpy as np, csv, os

# chron, top, base, D_above (normal chron above), D_below (normal chron below)
CAND = [
 ("C5Ar.2r",12.770,12.829,0.035,0.058),
 ("C5Ar.3r",12.887,13.032,0.058,0.151),
 ("C5AAr", 13.183,13.363,0.151,0.245),
 ("C5ABr", 13.608,13.739,0.245,0.331),
 ("C5ACr", 14.070,14.163,0.331,0.446),
]
# observed thicknesses (cm), boundary-midpoint based (polarity_zones.csv)
SC_R, SC_NUP, SC_NLO = 45.0, 103.0, 77.5     # R 104.5-149.5; N above 0-103 (open); N below 149.5-227 (open)
TAM_R, TAM_NLO = 50.0, 48.5                   # R 0-50 truncated above; N 51-99 (open)
TAM_HORIZON_ABOVE_R_TOP = 17.0                # Purussaurus cranial horizon, cm above measured R top
PALYNO=(12.7,14.2); FOSSIL=13.8

rows=[]
print(f"{'chron':8s} {'r_SC':>6s} {'ovr_Nup':>8s} {'ovr_Nlo':>8s} {'SCvar':>6s} "
      f"{'r_TAM_min':>9s} {'TAMovr@rSC':>10s} {'span_SC':>13s} {'palyno':>7s} {'fossil':>7s} {'horizon@C':>12s}")
for nm,T,B,Dup,Dlo in CAND:
    D=B-T
    r_s = SC_R/D                                 # cm/Myr forced by the two SC boundaries
    t_up, t_lo = SC_NUP/r_s, SC_NLO/r_s          # Myr spanned by flanking normals at uniform r_s
    ovr_up = max(0.0, t_up - Dup)                # overshoot into forbidden polarity (Myr)
    ovr_lo = max(0.0, t_lo - Dlo)
    # rate fold-change needed in the upper normal for it to fit its chron
    var_up = (SC_NUP/Dup)/r_s
    r_t_min = TAM_R/D                            # TAM R must fit inside the chron
    tam_ovr_at_rs = max(0.0, TAM_R/r_s - D)      # if TAM shared SC's rate
    sec_top = T - t_up; sec_bot = B + t_lo       # SC section span at uniform r_s
    pal_ok = (sec_top>=PALYNO[0]) and (sec_bot<=PALYNO[1])
    fos_ok = sec_bot<=FOSSIL+0.05
    # Purussaurus horizon age bracket at TAM: r in [r_t_min, 800 cm/Myr]
    h_young = T - TAM_HORIZON_ABOVE_R_TOP/ r_t_min   # slowest rate -> youngest
    h_old   = B - (TAM_R+TAM_HORIZON_ABOVE_R_TOP)/800.0
    rows.append(dict(chron=nm,r_SC=round(r_s),ovr_Nup_kyr=round(ovr_up*1e3),ovr_Nlo_kyr=round(ovr_lo*1e3),
        SC_rate_var_needed=round(var_up,2),r_TAM_min=round(r_t_min),TAM_ovr_at_rSC_kyr=round(tam_ovr_at_rs*1e3),
        SC_span=f"{sec_top:.2f}-{sec_bot:.2f}",palyno_ok=pal_ok,fossil_ok=fos_ok,
        horizon_bracket=f"{h_young:.2f}-{h_old:.2f}"))
    print(f"{nm:8s} {r_s:6.0f} {ovr_up*1e3:8.0f} {ovr_lo*1e3:8.0f} {var_up:6.2f} "
          f"{r_t_min:9.0f} {tam_ovr_at_rs*1e3:10.0f} {sec_top:6.2f}-{sec_bot:5.2f} {str(pal_ok):>7s} {str(fos_ok):>7s} "
          f"{h_young:6.2f}-{h_old:5.2f}")
out=os.path.join(os.path.dirname(os.path.abspath(__file__)),"joint_fit_results.csv")
with open(out,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nColumns: r_SC = rate forced by SC's two reversal boundaries (cm/Myr, compacted);")
print("ovr_Nup/ovr_Nlo = kyr by which SC's flanking normals overrun the adjacent normal chron at uniform rate;")
print("SCvar = fold-increase in rate the upper normal requires to fit; r_TAM_min = slowest TAM rate allowed;")
print("TAMovr@rSC = kyr TAM's R overruns the chron IF TAM shared SC's rate; horizon@C = Purussaurus age bracket.")
