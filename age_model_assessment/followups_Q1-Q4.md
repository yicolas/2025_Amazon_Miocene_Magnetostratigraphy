# Follow-ups: figure content, updated scores, calibrations, terrace model
_2026-07-22. Answers to four questions raised after the Purussaurus anchor._

---

## Q1 — What data would make the figures more understandable?

The polarity-column figures currently show *only* polarity + age. The biggest gains come from putting the **evidence and its uncertainty on the same axis** as the polarity:

**On the polarity/age panels**
- **Both fossil + palyno constraints as bands on every panel** (floor 12.7 Ma, ceiling ~13.8 Ma), not just the GPTS figure — so each site's admissible window is visible in place.
- **An age–depth line with sedimentation rate** annotated, and drawn as an **uncertainty band** (from the ±5 cm response function + zone-picking error → ±kyr), not a single line.
- **Reversal-boundary error bars** (the depth uncertainty of each N/R transition) — reviewers need to see how sharp each tie is.
- **The ODP 1237 reference polarity column** beside the sites, since 1237 is the local calibrated record you should tie *through* (rather than directly to the abstract GPTS).
- **Chron-duration ticks** on the GPTS column (label each chron's length in kyr) — this is what the pattern-fit below keys on.

**As adjacent tracks (share the depth axis)**
- **Inclination + MAD** (already have) — keep alongside.
- **Susceptibility / NRM intensity** downhole — flags the weak-mineralogy intervals (SC 7C ~100× drop; SC 2 diamagnetic) where directions are unreliable.
- **Lithology column** (clay / lignite / sand) and **XRF Fe/Ti** — shows which polarity wiggles are lithologic.
- **The Purussaurus horizon** drawn as a symbol at its stratigraphic height (~0.5 m above base) with its age bar.
- **A per-zone confidence strip** (n points, median MAD, does dec also flip, antipodality) so load-bearing boundaries are visually distinct from weak ones.

**Readability**
- Use a colourblind-safe fill (not pure black/white) and label N/R directly (done).
- Where records are short relative to the age axis (the multiscale figure), plot them on a **broken/zoomed age axis** so they aren't compressed to slivers.

I can rebuild the figures with any of these — the fossil+palyno bands, the age–depth band, and the ODP 1237 column are the top three.

---

## Q2 — Updated scores for TAM and SC

Two things changed how "score" should be read: (a) the fossil now supplies an admissibility filter, and (b) a **polarity-pattern-fit score** (how well the observed normal/reversed *thickness ratios* match a candidate chron) is far more diagnostic than the old wiggle-match correlation. Full table: `updated_scores.csv`.

**SC** (has a clean Normal–Reversed–Normal: 103 cm / 43 cm / 88 cm):

| Candidate reversed chron | Reversal age | Sed rate | Pattern misfit | Section age | Palyno OK | Fossil (≤13.8) |
|---|---|---|---|---|---|---|
| **C5ABr** (best pattern) | 13.61–13.74 Ma | 328 cm/Ma | **43 cm** | 13.29–14.08 Ma | ✓ | ✗ (base 14.08) |
| **C5AAr** (best fossil-compatible) | 13.18–13.36 Ma | 239 cm/Ma | 96 cm | 12.74–13.84 Ma | ✓ | ≈ (base 13.84) |
| C5Ar.2r | 12.89–13.03 Ma | 297 cm/Ma | 129 cm | 12.53–13.41 | ✗ | partial |
| C5Ar.1r | 12.77–12.83 Ma | 729 cm/Ma | 123 cm | 12.62–12.98 | ✗ | ✓ |

- **SC's polarity pattern points to C5ABr (~13.6–13.7 Ma)** — a short reversed zone between two long normals, a distinctly better fit (misfit 43 cm) than anything younger. This agrees with SC's *rank-1* wiggle-match model (13.28–14.30 Ma, r = +0.84), so two independent methods converge there.
- **The Purussaurus genus ceiling (≤13.8 Ma) demotes the pure-pattern winner** (C5ABr's section base runs to 14.08 Ma) and elevates **C5AAr (~13.2–13.4 Ma)** as the best fossil-compatible compromise.
- **Net: SC's reversal is C5AAr–C5ABr, ~13.2–13.7 Ma** — *older* than the ~13.0 Ma the earlier inclination wiggle-match alone suggested. This is the main revision from adding the fossil-aware pattern test.

**TAM** (only Reversed-top → Normal preserved; the upper normal was lost with the invalid Group 1):
- **Under-constrained** — with one boundary and no top normal, TAM can't be pattern-scored. It's compatible with C5AAr or C5ABr at plausible sed rates (~350–410 cm/Ma if its R is the same chron as SC's).
- Its value now is the **fossil**, which pins TAM's base ≤ 13.8 Ma. **Recovering TAM Group 1's polarity is the single fix that would give TAM a real score.**

**Bottom line:** the fossil didn't change the intrinsic waveform correlations, but it (1) ruled out the older half of SC's best pattern fit and (2) shifted the favored tie from ~13.0 to **~13.2–13.7 Ma (C5AAr/C5ABr)**. **Species ID is the tiebreaker:** *P. neivensis* (11.8–13.8 Ma) is compatible with both; only the circular "tight Pebas ~13.0" ceiling would force the younger, poorer-fitting option.

---

## Q3 — Other Miocene calibrations to compare against

GTS2020 ≈ GTS2012 (ATNTS2012) at 13 Ma (<10 kyr apart), so those two are one column. Worth adding:

| Calibration | Year | Basis | Δ vs ATNTS2012 at ~13 Ma |
|---|---|---|---|
| **CK92** (Cande & Kent) | 1992 | Marine magnetic anomalies | **~0.2–0.4 Myr older** (legacy) |
| **CK95** (Cande & Kent) | 1995 | Anomalies + one astro anchor | ~0.1–0.2 Myr (mostly younger) — *already in our figures* |
| **ATNTS2004 / GTS2004** (Lourens et al.) | 2004 | First astro-tuned Neogene GPTS | <0.05 Myr (parent of ATNTS2012) |
| **Abdul Aziz et al. (APTS, Orera, Spain)** | 2003 | Astro tuning of continental cycles | ≈0 — *this is the source* of the mid-Miocene ATNTS ages |
| **Mediterranean / deep-sea tunings** (Mourik 2011; Holbourn 2013–14) | 2010s | Astro-tuned δ¹⁸O across the MMCT | within a few kyr |
| **CENOGRID** (Westerhold et al.) | 2020 | Astro benthic isotope stack | consistent (climate scale, not a reversal revision) |
| GITS (Singer) | 2014 | Ar-Ar dated reversals | **Quaternary only — not applicable** |

**Practical recommendation:** add **CK92** as a third anomaly-based column (to bracket the legacy spread) and keep **CK95**; note that ATNTS2004, Abdul Aziz 2003, Holbourn, and CENOGRID are all effectively the same tie-points as ATNTS2012/GTS2020 in this window — so the *real* calibration spread your tie must survive is only **~0.1–0.4 Myr (CK92↔ATNTS)**, which is the same order as the reversal-to-boundary misfits. (Transcribe exact decimals from the primary tables before finalizing.)

---

## Q4 — Terrace model, or is it truly Miocene?

The classic fossiliferous Pebas exposures at Iquitos/Tamshiyacu are **genuinely middle Miocene in-situ Pebas** — there's no mainstream claim they're reworked or Pleistocene. But there *is* younger cover in the region, and a real stratigraphic-context question worth addressing explicitly:

**The younger units above the Pebas**
- **Nauta Fm** (Nauta 1 tide-influenced, Nauta 2 oxidized cross-bedded fluvial sands) — **Late Miocene**, gradational above Pebas (Rebata et al. 2006).
- **Iça / upper Solimões** (Brazil equivalent) — late Miocene–Pliocene fluvial.
- **Madre de Dios Fm** (Campbell et al. 2006) — **9.0–2.5 Ma** fluvial fill above the "Ucayali Unconformity."
- **Quaternary Amazon terraces / alluvium** — Pleistocene–Holocene, unconsolidated, tied to modern drainage.

**How to be sure your column is Pebas, not younger cover:** in-situ Pebas = lithified **blue-grey clays, lignite, pyrite, Pebasian mollusc coquinas, brackish/marine palynology**; Nauta/terrace = **oxidized/mottled clays, cross-bedded loose sands, no Pebasian fauna**. Your section already shows the Pebas signature (lignite, blue-grey clay, the palyno zonation, and now a Pebas-fauna vertebrate — *Purussaurus*), so it is Pebas, not terrace.

**Two caveats to name in the manuscript:**
1. **Campbell's Ucayali Peneplain model** posits paleosols capping a truncated Pebas and more late-Miocene–Pliocene fluvial cover than the Hoorn/Wesselingh continuous-succession model — so confirm each sampled interval is in-situ Pebas and not the gradational Nauta above it.
2. The **marine vs brackish-groundwater debate** affects the palynology criterion you use to confirm "Pebas," but not the Miocene age itself.

**So:** everything you're dating is Miocene Pebas; the terrace/Nauta/Madre-de-Dios units are real but sit *above* it, and the fossil + lithology + palynology together rule out a younger-terrace origin for your column.

## Sources
Timescales: Cande & Kent 1992 (*JGR* 97), 1995 (*JGR* 100); Lourens et al. 2004 (GTS2004); Abdul Aziz et al. 2003 (*JGR* 108); Holbourn et al. 2013–14; Westerhold et al. 2020 (*Science*, CENOGRID); Singer 2014 (GITS, Quaternary). Regional stratigraphy: Rebata-H. et al. 2006 (*J. South Am. Earth Sci.* 21; *Sedimentology* 53 — Nauta); Campbell et al. 2006 (*Palaeo³* 239 — Ucayali Peneplain / Madre de Dios Fm 9.0–2.5 Ma); Jaramillo et al. 2017 (*Sci. Adv.* — palyno T-15); Frontiers 2020 8:116 (saline-groundwater alternative).
