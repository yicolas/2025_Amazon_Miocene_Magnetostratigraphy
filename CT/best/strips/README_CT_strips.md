# CT core-face strips — Pebas TAM/SC u-channels

Source: Toshiba Aquilion CT of two boards of u-channels (side-by-side channels),
DICOM -> TIFF via the lab converter. Best coronal planes chosen from CT/best/.

## Scale (in-plane pixel spacing, from DICOM)
- SE1 = "Bone 2.0 Coron TOP.Ref"  : 0.615547 mm/px  (16.24 px/cm)  = UPPER board, channels 1-4
- SE2 = "Bone 2.0 Coron BTTM.Ref" : 0.468000 mm/px  (21.37 px/cm)  = LOWER board, channels 5-8
Strip length_cm in each filename is computed = height_px * mm_per_px / 10.

## Board layout (channels = paleomag groups)  [read from board photo]
UPPER board (SE1):
  Ch1 = TAM1 - TAM2 - TAM3B      -> TAM Group 1  (NRM-only, NO polarity; hosts Purussaurus lignite)
  Ch2 = TAM3A ... TAM5D          -> TAM Group 2  (polarity: mid)
  Ch3 = TAM6 - TAM7 (+5A/5B)     -> TAM Group 3  (polarity: base)
  Ch4 = SC2 - SC3C (+SC1x)       -> SC Group 4
LOWER board (SE2):
  Ch5 = SC3B - SC3A (+SC4x)      -> SC Group 5
  Ch6 = SC6 (+SC7x)              -> SC Group 6
  Ch7 = SC5                      -> SC Group 6/7
  Ch8 = empty

## Field-metre registration (from build_nb.py / polarity_synthesis)
TAM (base->top): Group3 0.00-1.13 m ; Group2 1.13-2.54 m ; Group1 2.54-3.78 m (no polarity)
  Purussaurus cranial ~2.71 m (in Group1 lignite).  Section total ~3.78 m.
SC : paleomag cored interval 1.60-5.00 m ; groups stack G4(top)->G7 by core lengths.
Orientation: use the "TOP" arrows on each channel to set up-section direction.

## Files
<TAG>_<PLANE>_b<k>_<len>cm.png : one channel band cropped from that plane, contrast-stretched.
extract_ct_strips.py           : regenerate/adjust strips from any chosen plane.
