# Internal correlations & multi-timescale age-model assessment
### Peru Amazon Miocene magnetostratigraphy (TAM, SC — Pebas Fm, ~12.5–14.5 Ma)
_Prepared 2026-07-22. Built on `polarity_synthesis/` (corrected polarity zones), `cross_parameter_comparison.csv` (18 candidate models), the ODP 1237 age-depth model, and the Amazonian palynozonation in `external_data/`._

---

## Bottom line

Each site records **exactly one well-resolved polarity reversal** flanked by normal polarity. That is enough to say the sections are mid-Miocene and internally consistent, but **not enough to fix a unique GPTS tie**: a single reversal matches any of ~8 mid-Miocene reversed chrons, and the palynology only narrows that to ~5. The current per-site "best" age models are *mutually inconsistent* and are being ranked by a score that rewards the least reliable parameter (declination). The single biggest improvement available is to **recover a second reversal in one continuous section** and to **fit both sites jointly** rather than independently.

---

## 1. What the polarity data actually support

| Site | Polarity (top→bottom) | Reversed zone | Quality of the reversal |
|------|----------------------|---------------|-------------------------|
| **TAM** | *(N lost)* → **R** → N → N | TAM-4, +8.3° incl | High: 50+ pts, MAD 4–7, gradual N→R transition through TAM 3A→4 |
| **SC** | N → **R** → N → N → N | SC-4C→4B/4A | **Unambiguous**: within one measurement group, both incl AND decl flip, ~157° antipodal dec offset, MAD ~4–5 on both sides |

Two things to carry forward:

- **TAM's uppermost normal is missing.** Group 1 (TAM 1/2/3B) has invalid PCA (MAD = 0; the 300–1000 mT window doesn't apply because that group was demagnetised with NRM + 10 mT only, and TAM 1 was mostly lignite). So TAM's true pattern is almost certainly **N–R–N**, but we only see **R–N** — we've lost the datum that would give TAM a *second* boundary.
- **SC-4C→4B is the single best reversal in the whole dataset** (same core group, same tray, excellent MAD, full inc+dec flip, near-antipodal). It is the anchor everything else should hang from. The TAM-4 reversed interval is its TAM counterpart.

---

## 2. Internal correlation (TAM ↔ SC)

The parsimonious reading is that **both sites record the same reversal** — TAM-4 ≙ SC-4C/4B — and the flanking normals correlate. The problem is that the *automated per-site fits contradict each other*:

| Model (inclination only) | Sed rate | Reversed zone → age | Nearest GPTS boundary (GTS2020) |
|--------------------------|----------|---------------------|--------------------------------|
| TAM best (rank 1) | 778 cm/Ma | 12.94–13.00 Ma | C5Ar.2n.base / C5AAn.top (≈ **C5AAr**) |
| SC rank 1 | 414 cm/Ma | 13.53–13.64 Ma | C5ABn.base (≈ **C5ABr**) |
| SC rank 2 | 786 cm/Ma | 13.08–13.14 Ma | **falls *inside* C5AAn** (a normal chron) — internal contradiction |

So the per-site optimisers put the *same physical reversal* at three different ages, and SC's rank-2 solution places its reversal inside a normal chron. Only the **coeval pair** (TAM 778 / SC 786 cm/Ma, top ≈ 12.94 Ma) is mutually self-consistent — and even that leaves SC's reversal ~50 kyr off the nearest boundary. **The per-site ranking is not a reliable guide; a joint fit is required** (see §7, item 1).

**Within-site continuity.** The large declination jumps at core boundaries (TAM 4→5E ~182°, TAM 5A→6 ~140°, TAM 6→7 ~188°, SC 7B) are **tray-rotation artifacts of unoriented cores, not reversals** — the diagnostic is that inclination keeps the same sign across them. A *real* reversal needs inclination **and** relative declination to flip together and be near-antipodal; only TAM 3A→4 and SC 4C→4B pass that test. Important nuance for the write-up: **absolute** declination is meaningless (unoriented), but **relative** declination *within a single measurement run* is usable and is exactly what confirms SC-4C→4B.

---

## 3. The 18 candidate age models — ranked by reliability, not by score

| Site | Parameter | Score (max) | Sed-rate spread | Verdict |
|------|-----------|-------------|-----------------|---------|
| TAM | Declination | **0.960** | 470–508 cm/Ma | **Discard** — highest score but unoriented cores make absolute dec meaningless |
| SC | Declination | **0.985** | 548–552 cm/Ma | **Discard** — same reason |
| TAM | Inclination | 0.963 | 774–778 cm/Ma | **Keep** — but correlation only +0.52 (weak) |
| SC | Inclination | 0.904 | 414–788 cm/Ma | **Keep** — correlation +0.78 to +0.84; ~2× sed-rate ambiguity |
| TAM | Intensity | 0.168 | 692–696 cm/Ma | Discard for age — lithologic, not geomagnetic |
| SC | Intensity | 0.328 | 310–328 cm/Ma | Discard — implies 15.5–16.8 Ma, geologically wrong |

The headline methodological problem: **the scoring metric rewards the physically invalid parameter.** Declination scores highest precisely because, with arbitrary per-core azimuths, the optimiser has extra free wiggle room to fit anything. Intensity "models" return ages up to 16.8 Ma because intensity tracks susceptibility/lithology (r ≈ 0.5–0.75), not the field. **Score ≠ reliability.**

---

## 4. Timescale sensitivity (GTS2020 vs GTS2012 vs CK95)

- **GTS2020 = GTS2012 in this window.** Both adopt the astronomically tuned ATNTS2012 for C5n–C5ADn, so the chron ages are *identical*. Switching between them changes nothing — worth stating explicitly so no one over-interprets a "GTS2020 update."
- **Cande & Kent 1995 is ~50–100 kyr younger** for these chrons. Under CK95 the same reversed zone that sits near C5AAn.top in GTS2020 shifts to near C5AAn.base.
- **This matters because the misfits are the same size as the differences.** The reversed-zone-to-boundary misfits in the best models are 30–77 kyr — the *same order* as (a) the GTS2020↔CK95 offset and (b) the durations of the short chrons themselves (C5Ar.1n is only 35 kyr long). A single reversal therefore **cannot discriminate among adjacent short chrons on any timescale.** Full multi-scale chron table: `gpts_multiscale.csv`.

---

## 5. Independent Amazon age control (what pins it, what doesn't)

- **Palynology (the useful constraint).** The section sits in zone **T-15** (*Crassoretitriletes vanraadshooveni*, 14.2–12.7 Ma) and the base of **T-16**; the **T15/T16 boundary at 12.7 Ma** is a firm *younger* limit. That rules out C5An (12.05–12.47 Ma) and everything younger, narrowing the reversal to **C5Ar.1nr … C5ACnr** (candidates 4–8 in `fig_nonuniqueness.png`).
- **ME-4 marine incursion (14.6–11.62 Ma)** spans the entire study window → permissive, adds no constraint.
- **No datable ash.** The Pebas is notoriously ash-poor, so there is no radiometric anchor; magnetostratigraphy + palynology + Sr-isotopes must carry the age control. **Sr-isotope ages on mollusc shells (Vonhof method) are the untapped independent numeric anchor.**
- **Net:** biostratigraphy narrows the tie to ~5 reversed chrons and a ≳12.7 Ma floor, but does **not** uniquely fix it.

---

## 6. Why it's still non-unique

A single "reversed zone flanked by normals" is compatible with any of ~8 mid-Miocene reversed chrons (`fig_nonuniqueness.png`); the palyno floor removes the youngest three. To collapse the remaining ambiguity you need **≥2 reversals in one continuous section** so the *spacing* between boundaries — not just their existence — does the matching. That is the crux of every recommendation below.

---

## 7. Ways to improve the approach

### A. Correlation method
1. **Fit the two sites *jointly*, with the shared reversal age as a hard constraint.** Solve simultaneously for (TAM sed rate, SC sed rate, common boundary age) instead of optimising each site alone. This removes the contradiction in §2 and is the single most impactful change. Your table already computes a "Joint" score — make coevality a constraint, not just a label.
2. **Recover TAM's lost top normal (Group 1).** Re-pick TAM 1/2/3B over the *correct* demag range (NRM + 10 mT), even at low resolution. Restoring the N above TAM-4 gives TAM a **second boundary** and turns a floating R–N into a diagnostic N–R–N.
3. **Use VGP latitude for polarity, and relative declination as a second discriminant.** Replace raw inclination-sign with the sign of VGP latitude (you already have VGP files), and use within-run relative declination to (a) detect/repair tray rotations and (b) independently confirm each reversal.
4. **Carry the foam-gap depth model (`viable_positions`) into the age fitting.** The polarity synthesis already registers foam gaps; make sure the age-model code uses that same sediment-depth axis rather than a raw concatenated one.

### B. Age-model fitting
5. **Rank by physical validity and boundary misfit, not by a unitless score.** Report each candidate as *(implied sed rate, boundary misfit in kyr, monotonicity check)*. Retire the "Joint score" as the primary ranking key.
6. **Drop declination and intensity from age control entirely.** Keep intensity only as a rock-mag/lithology proxy (it already correlates with susceptibility). Declination's high scores are an unoriented-core artifact.
7. **Propagate uncertainty and show the model as a band.** The reversal depth carries a ± (≈5 cm response function + zone-picking) error; convert to ±kyr and require a *monotonic, shared* sedimentation rate across both sites. A line implies precision you don't have.
8. **Report the tie on both GTS2020 and CK95.** State that it is robust to the 2012↔2020 choice and shifts ~50–100 kyr under CK95 — this is a cheap, honest robustness statement for the manuscript.

### C. New constraints to break non-uniqueness (highest scientific value)
9. **Get a second reversal into one continuous section** (extend/splice cores, or exploit ODP 1237). Two boundaries with a plausible sed rate usually pin a *chron pair* uniquely — this is worth more than any reprocessing.
10. **Add Sr-isotope (⁸⁷Sr/⁸⁶Sr) ages on Pebas mollusc shells.** Given zero datable ash, even one good Sr age converts the floating polarity model into an *anchored* one.
11. **Put the palyno datums on the core depth axis.** You have taxa counts + the composite section; place the T15/T16 (12.7 Ma) and T14/T15 datums as *depth* tie-points, so biostratigraphy becomes a hard tie rather than a permissive age range.
12. **Tie to ODP 1237 as the regional reference, not directly to the abstract GPTS.** 1237 gives a local, well-dated reversal + age-depth stratigraphy (`age_depth_model_1237.xlsx`); pattern-match the Amazon reversal to 1237 first, then let 1237 carry the GPTS calibration.

### D. Data quality
13. **Mask weak-mineralogy intervals explicitly.** SC 7C (~100× susceptibility drop) and SC 2 (diamagnetic) carry unreliable directions that low MAD alone won't flag; exclude them from polarity picking and say so.
14. **Publish a per-zone confidence table** (n points, median MAD, does inc *and* dec flip, antipodality angle) beside the polarity column, so reviewers can see which boundaries are load-bearing.

---

## Files produced
- `fig_internal_corr.png` — TAM vs SC polarity on sediment depth (the shared reversal)
- `fig_multiscale_gpts.png` — observed polarity vs GTS2020 / GTS2012 / CK95 + biostrat bands
- `fig_nonuniqueness.png` — the single-reversal tie is non-unique; biostrat narrows to ~5 chrons
- `gpts_multiscale.csv` — chron boundary ages, three timescales
- `age_scenarios.csv` — every candidate model → implied reversal ages + misfit to nearest boundary, both timescales
- `age_model_analysis.py` — reproducible source

## Key sources
GTS2012/ATNTS2012: Hilgen, Lourens & Van Dam (2012), *Geologic Time Scale 2012*, ch. 29 (verified vs IODP Exp. 342 Table T7). GTS2020: Ogg (2020) GPTS — adopts ATNTS2012 unchanged for the mid-Miocene. CK95: Cande & Kent (1995), *JGR* 100(B4). Amazon palynozonation: Jaramillo, Rueda & Torres (2011), *Palynology* 35(1). Marine incursions (ME-1–ME-5): Parra et al. (2026), *Fossil Studies* 4(1). Sr-isotope method: Vonhof et al. (1998, *Palaeo3* 141; 2003, *GSA Bull.* 115). Contamana age framework: Antoine et al. (2016), *Gondwana Research* 31.
