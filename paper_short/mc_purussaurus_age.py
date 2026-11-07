#!/usr/bin/env python3
"""Monte-Carlo inferred age of the Tamshiyacu Purussaurus, propagating the
U-Pb-calibrated La Venta Bayesian age model (Mora-Rojas et al. 2023) through the
shared C5AAr magnetostratigraphic tie.

Logic:
- The Purussaurus cranial (TAM 13) is at 2.71 m, ~0.17 m ABOVE the top of the
  measured C5AAr reversal (2.54 m = C5AAr/C5AAn boundary level) in the no-polarity
  Group-1 lignite. It is therefore bracketed by the C5AAr–C5AAn interval.
- Its chron-equivalent position in the La Venta section is the C5AAr-top region:
  176.6 m in the Mora-Rojas et al. (2023) composite section (their s15 polarity column,
  digitised to age_model_assessment/laventa/laventa_cs_polarity.csv), drawn -20/+15 m
  (156.6-191.6 m) to reflect the no-data gap + cross-basin mapping.
- The La Venta age model (Table 4) gives median age + 95% CI vs La Venta metre;
  we sample it, then (optionally) truncate to the magnetostratigraphic bracket.
"""
import numpy as np
rng=np.random.default_rng(42)

# La Venta age model, Mora-Rojas 2023 Table 4: (metre, median_Ma, CI_young, CI_old)
LV=np.array([
 (63.5,16.00,15.94,16.00),(120,15.20,13.18,16.00),(154,13.26,12.72,13.77),
 (156.6,13.13,12.70,13.46),(205,13.13,12.70,13.46),(218.8,13.13,12.69,13.47),
 (242,13.13,12.68,13.45),(249,13.11,12.68,13.46),(303,12.40,11.81,12.92),(329,12.40,11.82,12.91)])
d,med,lo,hi=LV[:,0],LV[:,1],LV[:,2],LV[:,3]
order=np.argsort(d); d,med,lo,hi=d[order],med[order],lo[order],hi[order]
sig=(hi-lo)/3.92   # 95% CI -> 1 sigma

# ATNTS2012 (Hilgen et al. 2012) chron ages
C5AAr=(13.183,13.363); C5AAn=(13.032,13.183)

N=200000
# (1) chron-equivalent La Venta depth for the Purussaurus (at/just above C5AAr top, 176.6 m)
LV_C5AAR_TOP=176.6
depth=rng.uniform(LV_C5AAR_TOP-20,LV_C5AAR_TOP+15,N)
mu=np.interp(depth,d,med); s=np.interp(depth,d,sig)
age_lv=rng.normal(mu,s)                       # pure La Venta age-model estimate

# (2) magnetostratigraphic bracket: Purussaurus in C5AAr–C5AAn (13.032–13.363)
mask=(age_lv>=C5AAn[0])&(age_lv<=C5AAr[1])
age_comb=age_lv[mask]                          # La Venta model ∩ magnetostrat bracket

def rep(name,a):
    p=np.percentile(a,[2.5,16,50,84,97.5])
    print(f"{name:32s} median {p[2]:.2f}  68% [{p[1]:.2f},{p[3]:.2f}]  95% [{p[0]:.2f},{p[4]:.2f}]  n={len(a)}")

print("Purussaurus inferred age (Ma):")
rep("La Venta age-model only",age_lv)
rep("La Venta model ∩ C5AAr–C5AAn",age_comb)
# pure magnetostratigraphy (uniform over the bracket), for comparison
rep("magnetostratigraphy only",rng.uniform(C5AAn[0],C5AAr[1],N))
import os; os.makedirs("figs",exist_ok=True); np.save("figs/puru_age_samples.npy",age_comb)
print("\nPreferred: La Venta model ∩ magnetostrat bracket.")
