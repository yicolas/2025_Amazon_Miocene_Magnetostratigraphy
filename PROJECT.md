# Amazon — Miocene Magnetostratigraphy

**Status:** Active analysis
**Stage:** u-channel data processed, age model in progress
**Target:** _(journal / deadline — fill in)_
**Manuscript:** _(not started)_

<!-- Edit this file in any text editor. The dashboard parses these sections.
     Checkboxes: "- [ ]" = open, "- [x]" = done. -->

## Session handoff (2026-07-15)

**Project:** Miocene magnetostratigraphy of Peru Amazon u-channel cores (TAM, SC) tied to ODP 1237 / the GPTS; Pebas Fm, ~12.5–14 Ma.
**Findings:** Raw sign-of-inclination "reversals" are noise (they cluster at low intensity); intensity drops are lithologic (NRM tracks susceptibility r≈0.5–0.75). The PCA ChRM fits already exist (`OSU_uchannel_magnetometer_NRM_data/dir/pca_corrected/` 300–1000; despiked under `ucheck_spike_removed/`) and give shallow, coherent, multi-polarity directions → primary Miocene DRM; the soft (100–300) component is a steep **coring overprint** (an artifact, not datable).
**In the notebook (`Uchannel_notebook.ipynb`, two new bottom sections):** "Cleaned polarity synthesis" and "Corrected polarity (S-hemisphere sign + de-noising + foam-gap depths)". KEY FIX: at ~4°S, normal polarity = NEGATIVE inclination (`POLARITY_SIGN=-1`). Corrected result — SC = dominantly Normal with one Reversed zone (NNRNNN); TAM = RNN.
**Next:** constrain `AGE_MODEL` with Pebas palynology zonation + ash dates (one reversal ties non-uniquely to the GPTS); correlate TAM↔SC internally; fold XRF Fe/Ti onto the depth axis to confirm the reversed / low-intensity interval is lithologic.
**Data:** biostrat/age supplements + models copied from iCloud into `external_data/` (taxa counts, composite section, CONOP, age-depth model, magnetostrat summary) — these are the age-anchor inputs.
**Note:** cores are unoriented → do NOT use declination for polarity.

## To-do

- [ ] Settle the polarity interpretation — sign-of-inclination is not enough, needs PCA / component fitting
- [ ] Pick the preferred age–depth model from `cross_parameter_comparison.csv` (18 candidates)
- [ ] Tie ODP 1237B section (200–319 mbsf) to the GPTS
- [ ] Consolidate the three u-channel notebooks (`Uchannel_notebook`, `_clean`, `_safety`) into one
- [ ] Commit the CT stack TIFFs (or gitignore them if they're too big)
- [ ] Integrate Pebas-XRF and susceptibility data into the interpretation

## Manuscript checklist

- [ ] Abstract
- [ ] Intro
- [ ] Methods
- [ ] Results
- [ ] Discussion
- [ ] Figures
- [ ] References

## Figures needed

- [ ] Inclination + declination vs depth with polarity column
- [ ] Correlation to GPTS
- [ ] Age–depth model with sed-rate
- [ ] Rock mag / susceptibility downhole

## Open questions

**Are the reversals real or overprint/noise? — the RAW sign-flips are noise; the FITTED ChRM is real.**
Two different things:
- The **raw single-step NRM** sign-flips ARE noise: they cluster at low intensity (median
  I(flip)/I(stable)=0.54, down to 0.07 for SC_1C/SE_3B), with implausible ±89° inclinations from a
  vertical coring overprint and 19–131 flips per core. Don't use raw sign-of-inclination.
- The **PCA ChRM fits already exist** (`dir/pca/` 100–1000, `dir/pca_corrected/` 300–1000, and the
  despiked versions under `ucheck_spike_removed/`). The cleanest set (despiked + 300–1000) gives
  **shallow inclinations** (median −18° to +4°, right for ~4°S) and **1–3 coherent polarity zones**
  per core (RNR / RN / RRN / NNR …), 1–7 flips, MAD 5–17. **This looks like real mid-Miocene
  magnetostratigraphy**, consistent with the in-situ Miocene fossils.

**Primary vs overprint.** The hard magnetization (AF>300) is shallow AND multi-polarity → a
post-depositional single overprint can't produce several reversals, so the ChRM records the field
through deposition = **primary DRM**. The **soft (100–300) component is steeper** than the ChRM
(removing it shallows incl: SC_1C 29°→20°) = a near-vertical **coring/drilling overprint** (a
mechanical artifact, not a datable remagnetization), already stripped by the 300–1000 window.
Dating the overprint by direction is not feasible: at ~4°S the field is shallow for all recent
epochs and the u-channels are **unoriented** (declination arbitrary). The datable signal is the
primary DRM via its **reversal pattern**, not the overprint.

**Does lithology explain the intensity drops? — YES, they are lithologic (measured).**
NRM intensity **tracks magnetic susceptibility** (log-intensity r = 0.49–0.75 in 3 of 4 matched
cores; TAM_3A weak at 0.23). Susceptibility ∝ magnetic-mineral concentration, so the intensity
drops are **real sediment-property changes** (magnetite-poor, likely clay-/organic-rich intervals),
not instrument artifacts. This closes the loop with the reversal question: the lithologically-
controlled low-intensity intervals are exactly where the polarity signal is noise-dominated.
XRF (Fe/Ti for magnetic minerals; clay/detrital proxies) would pin the lithology — see below.

## Focus now

- [x] **Use the existing 300–1000 despiked PCA fits** (`ucheck_spike_removed/dir/pca_corrected/`) as
  the polarity source — they're already shallow + coherent. Filter MAD<10–15, trim ~5 cm off each
  sub-core end (u-channel edge effects), and merge into polarity zones at ~5 cm response resolution.
  → Done 2026-07-22: reproduced as standalone `polarity_synthesis/polarity_synthesis.py` (MAD<12,
  5 cm edge-trim, foam removed, 7 cm median, zones ≥12 cm, POLARITY_SIGN=-1). Outputs in
  `polarity_synthesis/` (samples_corrected.csv, polarity_zones.csv, figures). Result: **TAM = RNN**,
  **SC = N R N N N** (dominantly Normal, one Reversed zone at ~106–149 cm composite).
- [~] Build a **composite stratigraphic section** (stack the sub-cores in order, common depth) then
  match the **reversal-zone pattern to the GPTS** (pattern-match the boundaries — NOT point-by-point
  waveform correlation, and NOT declination: cores are unoriented so absolute Decl is meaningless).
- [ ] Anchor with independent tie points: **ash layers** (`age_depth_model_1237`, `ash_layers_time`)
  and **palynology zonation**; require a monotonic, geologically reasonable sed rate.
- [ ] Age-model note: intensity models score poorly (Pearson 0.17–0.33); the high **declination**
  scores are likely spurious (unoriented cores) — trust **inclination polarity + reversal pattern**.
- [ ] XRF: consolidate the 27 sections onto depth and test Fe/Ti/clay vs the intensity minima
  (intensity already tracks susceptibility r≈0.5–0.75 → drops are lithologic).

## Closed / resolved

- ~~Are the low-intensity reversals real?~~ — no, they concentrate at the noise floor (see above).
- ~~Do intensity drops have a lithologic cause?~~ — yes, intensity tracks susceptibility.

## Log

- 2026-09-10 — **Discussion reframed to Hoorn mega-wetland (Everglades analog).** Recast the marine
  signal as episodic/subordinate against a dominantly freshwater-to-oligohaline background (endemic
  mollusc fauna = long-lived freshwater body); pyrite/Ca beds = brief brackish ingressions +
  condensation-favoured pyritisation. New subsection "Sedimentation rate, paleoclimate, and the road
  to the modern Amazon": ~7.8 m/Myr (both sites) is EXPECTED for a low-gradient sediment-starved
  mega-wetland (not anomalously low), consistent with tens of m/Myr regional Pebas/Solimoes; full-
  C5AAr forcing would need implausible ~2.4 m/Myr -> condensed record w/ minor diastems. Added lake-
  wetland -> transcontinental-Amazon arc (fan onset ~11-9.4 Ma). New refs Figueiredo et al. 2009,
  Hoorn et al. 2017. Rebuilt docx/pdf (15 pp). Commit follows 55fa188.

- 2026-09-10 — **Discussion + Results written; XRF paleoenvironment.** Added full Discussion
  (four-way age convergence on C5AAr ~13.2–13.4 Ma; magnetic mineralogy & reductive-diagenetic
  overprinting — XRF total Fe is authigenic, not detrital magnetite, so it cannot proxy magnetic
  concentration; redox/pyrite/Fe2+/sulfur + Ca-carbonate as marine-influence proxies in the Pebas
  mega-wetland; MMCT paleoclimate context; limitations). New "Elemental stratigraphy" results
  subsection + `paper/figs/fig_xrf_paleoenv.png` (SC composite downcore Fe,S,Ca,Fe/Ti,K/Ti,Mn +
  polarity). Script `age_model_assessment/xrf/xrf_paleoenv.py`. Rebuilt docx/pdf (14 pp).
  **Sedimentation-rate check:** barcode+GPTS age model gives ~7.8 m/Myr (compacted) at BOTH sites
  (SC 262 cm/333 kyr; TAM 99 cm/127 kyr) = ~16–23 m/Myr decompacted — low end for an active wetland;
  consistent with a condensed, low-gradient, occasionally marine-influenced lacustrine setting and/or
  uncored gaps between u-channels (barcode rate is a MINIMUM). A full-C5AAr correlation of the 43 cm
  SC-R zone would imply only ~2.4 m/Myr (implausibly condensed) — supports partial-chron capture.

- 2026-09-10 — **Core-azimuth reconstruction (orientation determined).** Applied field outward-face
  azimuths (TAM per-core 188–299; SC 223; bank 229) to raw Decl_core; convention
  D=Decl_core+az+C with global offset C≈35° (absorbs the up-face reference offset; confirm x/y-face to
  fix physically). VGP concentration 24%→69% within 30° of a pole; directions collapse onto GAD stars
  (`age_model_assessment/orientation/`). Half the sub-cores needed a binary 180° LAB tray-flip
  (resolved by inclination → doesn't change polarity). **Confirms** TAM4=R, TAM5/6=N, SC-4C→4B reversal
  (now antipodal in dec too); **TAM7 = Reversed** for any plausible azimuth (238–300°) → possible 2nd
  reversal at TAM base (tentative, low quality). Inclination still primary. Added equal-area stereonet +
  polar VGP figures; **manuscript updated** (diagnostic VGP fig replaced with reconstructed ones;
  abstract/reliability table/TAM7 revised; 12 pp). Notebook → 129 cells.
- 2026-09-10 — **La Venta regional anchor + barcode matching.** Extracted Mora-Rojas et al. 2023
  (Geodiversitas 45(6)) + Carrillo et al. 2023 from Downloads. New U-Pb + Bayesian ages: Honda Gp
  16→10.5 Ma; La Victoria <16–12.58, Villavieja 12.58–10.52; magnetostrat reinterpreted C5ABn→C5An.1n
  (ATNTS2012). **Monkey Beds / P. neivensis = chrons C5AAr–C5Ar.2n (~12.83–13.36 Ma)** — a U-Pb-anchored,
  non-circular ceiling that overlaps SC's pattern-fit reversal (C5AAr). Carrillo notes Laventan beds at
  Tamshiyacu itself (Antoine pers. comm.). Built local-vs-regional-vs-GPTS **barcode figure** + **VGP
  diagnosis figure**; U-Pb table in `age_model_assessment/laventa/`. **Manuscript now has ALL figures
  embedded** (polarity, VGP, XRF×2, fossil bracket, barcode) + new "Regional correlation to La Venta"
  section + Mora-Rojas/Carrillo refs; docx/pdf rebuilt (11 pp). Notebook → 127 cells.
- 2026-07-24 — **XRF processed** (`age_model_assessment/xrf/`; 27 clean `Results.txt`, 1906
  valid pts; note: lowercase `result.txt` are STALE/mislabeled — use `Results.txt`). **Key
  finding (corrects earlier expectation):** XRF *total Fe* does NOT track NRM intensity (pooled
  r≈+0.16; Fe/Ti +0.08; not a registration artifact — lag test peaks at 0, 5-cm-binned +0.18).
  Fe instead tracks **S (pyrite, r≤+0.69)** and **Ca (carbonate, r≤+0.88)** → most Fe is
  authigenic, not detrital magnetite. So intensity minima are lithologic/diagenetic
  (pyritisation dilutes magnetite) but diagnosed by SUSCEPTIBILITY, not XRF Fe. XRF section of
  the paper rewritten with this; XRF cell added to notebook (now 125 cells).
- 2026-07-23 — **Paper drafting started** (`paper/pebas_magnetostratigraphy.{tex,docx,pdf}`):
  Methods, Paleomagnetism, XRF, and Fossil sections drafted as living .tex (source) + .docx
  (pandoc-synced). Placeholders flagged for affiliation, XRF instrument model, SC locality.
- 2026-07-23 — **TAM Group 1 = unrecoverable.** Raw `.DAT` shows only 2 steps (NRM + one
  ~12.5 mT/125 Oe); no high-coercivity series → no ChRM fit possible (PCA degenerate, MAD=0);
  directions steep+scattered (overprint). TAM's top normal cannot be recovered without
  re-measuring with a full AF demag series. (Confirms the earlier "recover Group 1" item is
  blocked by the data, not the analysis.)
- 2026-07-23 — **VGP scatter diagnosed** (not a data fault). VGPs don't cluster at the poles
  because cores are unoriented (one fixed azimuth applied per multi-sub-core group) AND at ~4°S
  the shallow inclinations make VGP latitude hypersensitive to the unknown declination (a dec
  sweep 0–180° swings VGP_lat +89°→−83°). → Do not use VGPs for polarity/QC here; inclination
  sign only. Written into notebook + methods.
- 2026-07-23 — **New analyses integrated into `Uchannel_notebook.ipynb`** (10 cells appended in
  place on device; backup `.bak_before_append` made; 113→123 cells): polarity synthesis, TAM-1
  diagnosis, VGP diagnosis, fossil/age bracket. Helper files in `.notebook_append/`.
- 2026-07-22 — **Purussaurus fossil anchor** (found at Tamshiyacu ~0.5 m above column base).
  `age_model_assessment/fossil_anchor/`. Genus FAD ~13.8 Ma → section is YOUNGER than ~13.8 Ma:
  the independent OLDER bound the palyno floor (>12.7 Ma) lacked. Together they bracket the column
  to **~12.7–13.8 Ma** (→ ~12.7–13.0 Ma with the tight Iquitos-Pebas assemblage age, ~13 Ma,
  Salas-Gismondi 2015). Cuts candidate reversed chrons 6→3 and supports the C5Ar/C5AAn (~13 Ma)
  tie — converges with the best-fit inclination models. Caveats: needs species ID (only
  *P. neivensis* is diagnostically mid-Miocene; late-Miocene species would overturn it); watch
  circularity (lean on La Venta, not Pebas, for the ceiling); fossil sits in the invalid TAM
  Group 1 interval — re-pick its polarity to pin a zone.
- 2026-07-22 — **Internal-correlation & multi-GPTS age-model assessment** (`age_model_assessment/`).
  Each site has exactly one well-resolved reversal (SC-4C→4B is the cleanest; TAM-4 its match).
  Per-site best-fit models are mutually inconsistent → need a JOINT fit. GTS2020==GTS2012 here;
  CK95 ~50–100 kyr younger. Scoring rewards declination (spurious, unoriented) and intensity
  (lithologic) — rank by physical validity + kyr-misfit instead. Single reversal ties non-uniquely
  (~5 chrons after biostrat). Top improvements: joint two-site inversion; recover TAM Group 1
  normal; add a 2nd reversal / Sr-isotope age; drop dec+intensity from age control.
- 2026-07-22 — Polarity synthesis extracted from the notebook into a standalone, re-runnable
  script (`polarity_synthesis/`). Consumes the despiked 300–1000 PCA fits; foam-gap depths via
  `viable_positions`; S-hemisphere sign (normal = neg incl). Composite polarity **TAM = RNN**,
  **SC = NRNNN** (one Reversed zone, ~106–149 cm). Emitted tidy CSVs (per-sample + zones) and
  publication-style figures (incl+MAD+polarity per site, TAM-vs-SC, provisional GPTS tie). Fixed
  an inherited bug in the GPTS-comparison age window that clipped the SC record to the overlap.
  Next: TAM↔SC internal correlation and picking the age–depth model (`cross_parameter_comparison.csv`).
- 2026-07-15 — Reviewed the project. Measured reversal-vs-intensity (flips cluster at low
  intensity, ratio 0.54) → reversals are largely noise/overprint; and intensity-vs-susceptibility
  (r 0.5–0.75) → intensity drops are lithologic. Answered both open questions; added Focus-now.
- 2026-04-14 — notebook updates and outputs
- 2026-02-16 — CT, XRF, magnetic data added
