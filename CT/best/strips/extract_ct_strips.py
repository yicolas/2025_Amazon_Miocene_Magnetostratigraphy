#!/usr/bin/env python3
"""Extract per-channel core-face strips from a chosen CT coronal plane.
Each vertical bright band = one board channel (= one paleomag group). Adjust THR."""
import tifffile, numpy as np, sys, os
from PIL import Image
PXMM={"SE1":0.615547,"SE2":0.468000}   # mm/px from DICOM PixelSpacing
def extract(path, series, thr_frac=0.28, outdir="."):
    a=tifffile.imread(path).astype(float); mmpp=PXMM[series]
    colmax=a.max(0); thr=colmax.max()*thr_frac; on=colmax>thr
    runs=[]; i=0; N=len(on)
    while i<N:
        if on[i]:
            j=i
            while j<N and on[j]: j+=1
            runs.append((i,j)); i=j
        else: i+=1
    k=0
    for c0,c1 in runs:
        if c1-c0<25: continue
        band=a[:,c0:c1]; rm=band.max(1)>thr
        if rm.sum()<40: continue
        r0,r1=np.where(rm)[0][[0,-1]]; band=band[r0:r1+1]
        if band.mean()<colmax.max()*0.05: continue
        lo,hi=np.percentile(band[band>0],[2,99.5])
        b=np.clip((band-lo)/(hi-lo)*255,0,255).astype('uint8'); k+=1
        L=band.shape[0]*mmpp/10.0
        Image.fromarray(b).save(os.path.join(outdir,f"{series}_{os.path.basename(path).split('_')[0]}_b{k}_{L:.0f}cm.png"))
    print(f"{path}: {k} bands")
if __name__=="__main__":
    # example: python extract_ct_strips.py SE1 CT/best/SE1/0031_IM31_8bit.tif
    extract(sys.argv[2], sys.argv[1])
