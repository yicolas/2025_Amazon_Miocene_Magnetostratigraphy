# Core-azimuth reconstruction — findings
_2026-09-10. Field outward-face azimuths applied to the ChRM to recover geographic declination and test VGPs._

## Inputs (field, azimuth/dip; dip=3 flags the outward-face direction; true bedding ~0)
Bank strike ‖ river = 229. TAM 1=188, 2=241, 3=260, 4=299, 5=285, 6=276, 7=unknown.
All Santa Clara (SC) = 223. Correction convention (verified against the pipeline):
`D_geo = Decl_core + azimuth + C`, with one global reference offset **C ≈ 35°** fit here.

## Geometry note (from field/lab convention)
The measured "outward face" is the face that pointed **up (toward the ceiling)** in the
magnetometer — likely the sample **x- or y-face**. So the field azimuth is the compass direction of
the magnetometer's up-facing reference face, and the core→geographic transform is a **single fixed
rotation** — which is what the global offset **C ≈ 35°** absorbs here. If the up-face is confirmed as
specifically the x- or y-face, C can be replaced by a *physically derived* value (a multiple of 90°
plus the ~−8° IGRF declination) instead of being fit — removing one free parameter and sharpening the
per-core residuals below. (Worth checking the notes: C≈35° is not a clean 0/90/180/270+IGRF, hinting
either at an x/y-face detail or a few degrees of systematic field-compass offset.)

## What the reconstruction does — and doesn't — fix
- **It helps.** VGPs move toward the poles: fraction within 30° of a pole rises **24% → 69%**.
- **A residual per-core 180° "lab flip" remains.** About half the TAM sub-cores need an extra
  180° rotation to fit a GAD field. This is the *unrecorded tray-placement flip* we flagged earlier:
  the field azimuth records how the core sat in the outcrop, but **not** how it was loaded into the
  magnetometer. It is a binary artifact, resolvable by inclination (a 180° declination flip does not
  change inclination), so it does **not** change inclination-based polarity.
- **Net:** orientation reconstruction **confirms** the inclination polarity and gives real geographic
  declinations, but because of the flip ambiguity + shallow inclinations, declination/VGP still cannot
  *independently* set polarity — inclination stays primary. Orientation is now a validated cross-check.

## Per-core validation (after field azimuth + C + binary flip)
| Sub-core | n | meas az | Inc (pol) | lab flip | residual | VGP lat |
|---|---|---|---|---|---|---|
| TAM4  | 58 | 299 | +7.4 (R) | 180 | **9°** | −83 (clean reversed) |
| TAM5A | 2  | 285 | −2.1 (N) | 180 | 10° | +80 |
| TAM5E | 4  | 285 | −8.9 (N) | 0   | **5°** | +85 (validated, no flip) |
| TAM5B | 9  | 285 | −4.0 (N) | 180 | 25° | +65 |
| TAM5C | 7  | 285 | −13.8 (N)| 180 | 38° | +58 |
| TAM5D | 3  | 285 | −6.9 (N) | 0   | 34° | +56 |
| TAM6  | 48 | 276 | −3.1 (N) | 0   | **61°** | +30 |

- **TAM4, TAM5A, TAM5E validate to ≤10°** — the field azimuths are real and usable.
- **TAM6 is the outlier (61°)** despite excellent MAD — its field azimuth (276) or its declination is
  anomalous. Worth re-checking the TAM6 photo/measurement; its polarity (normal, Inc −3) is unaffected.
- SC used a single azimuth (223) for all cores, so its (mixed-polarity, shallow) VGPs stay smeared —
  SC would need per-sub-core azimuths to de-rotate individually.

## TAM 7 solution
TAM 7 ChRM (despiked, MAD<25): median **Inc = +4°** → inclination polarity **Reversed**.
Across the whole geologically plausible azimuth range (≈238–300°, i.e. at or near its neighbours
TAM6=276, TAM5=285), TAM 7's VGP is **reversed** (negative VGP latitude); the cleanest reversed pole
is at azimuth **≈238°**. Only an implausible azimuth near 50–90° (far from its neighbours) would make
it normal. **Conclusion: TAM 7 is Reversed for any orientation consistent with its neighbours.**

Caveat: TAM 7 is low-quality (MAD>15, shallow Inc ≈ +4°), so this is tentative; the field
photo/azimuth would confirm it, but polarity-wise the reversed reading holds across the plausible range.

## Implication for the polarity record
Orientation **confirms** TAM 4 = Reversed and TAM 5/6 = Normal, and now indicates **TAM 7 = Reversed**.
That makes the base of TAM potentially R–N–R rather than R–N — a possible second reversal low in TAM,
though it rests on low-quality TAM 7 data and should be treated as tentative pending better demag.

## Files
`fig_orientation_before_after.png`, `fig_tam7_and_azcheck.png`, `per_core_azimuth_check.csv`,
`reconstructed_directions.csv`, `orient.py`.
