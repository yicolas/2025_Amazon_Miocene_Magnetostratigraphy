# Pebas paper figures — regeneration bundle

Run `Pebas_paper_figures.ipynb` (set BASE + OUT in the config cell, Run All) to regenerate every
manuscript and supplement figure. All scripts read the repo root from `$PEBAS_BASE` and write to
`$PEBAS_OUT` (the notebook sets both).

Manuscript figure -> script:
- Fig 1  study-area map        : EXTERNAL (QGIS) — not scripted; keep fig1_location.png in the repo
- Fig 2  master correlation    : fig2_master.py         (GPTS panel labelled GTS2020)
- Fig 3  paleomag/rock-mag     : fig3_final.py  (+ fig3_rockmag_dense.py, fig3_bottom_candidates.py)
- Fig 4  accumulation/Sadler   : fig4_combined.py       (+ fig_sadler.py)
- Fig S1 GRM/MDF greigite screen: make_supp_figs.py
- Fig S2 reversal stereonet    : make_supp_figs.py
- Fig S3 Zijderveld            : make_supp_figs.py
- Fig S4 alt-GPTS exclusion    : fig_altcorr.py

Helpers: strat_symbols.py + symbols/*.svg (used by fig2_master.py).

Notes:
- Fig 2 GPTS header is GTS2020 in the SOURCE (fig2_master.py), not just the exported PNG.
- Fig 3 panel D = SC inclination -> paleolatitude (orientation-independent); VGPs removed
  (not valid for azimuthally unoriented u-channels).
