'''
detrend.py
Created by Xue, J. C.
Parameters:
    fils: preliminarily aligned files by flca.py
    ref: index of aligned files with other reference files
    xm: shift-X of the aligned files with the same number of elements as ref
    ym: shift-Y of the aligned files with the same number of elements as ref
    outdir: output directory ended with a slash "/"
    interpolate: order of transform.warp.
         - 0: Nearest-neighbor
         - 1: Bi-linear (default)
         - 2: Bi-quadratic
         - 3: Bi-cubic
         - 4: Bi-quartic
         - 5: Bi-quintic
'''
from glob import glob
import os
import numpy as np
from skimage import transform
from astropy.io import fits
def readfits(f):
    hdul = fits.open(f)
    data = hdul[0].data
    hdr = hdul[0].header
    hdul.close()
    return data, hdr

# parameters:
fils = glob('/home/xuejc/data/181110/CENT/test_align/Ha*fits'); fils.sort()
ref = [0,  10,  13,  31,  37,  51,  59,  66,  72,  78,   88,  93,  101,  127,  143,  170,  184]
xm = [-0.,-0.8,-0.6,-0.4,-1.8,-2.5,-2.7,-7.7,-8.7,-10.8,-9.2,-10.4,-10.8,-15.0,-17.6,-23.3,-27.0]
ym = [-0.,0.7, 1.6, -0.3,-0.5,-4.0,-7.0,-5.2,-6.5,-7.0, -9.4,-8.40,-8.9, -14.2,-19.3,-20.7,-24.7]
outdir = '/home/xuejc/data/181110/CENT/test_align/'
interpolate = 1

j = 0
if ref[0] > 0:
    for i in range(ref[0]):
        data, hdr = readfits(fils[i])
        tform = transform.SimilarityTransform(translation=(xm[0],ym[0]))
        print(i, tform.params[0,2], tform.params[1,2], end='\r')
        data = transform.warp(data, tform, order=interpolate)
        fname = os.path.basename(fils[i])
        mysp.fitswrite(outdir+fname, data, hdr)
for i in range(ref[0], ref[-1]+1):
    data, hdr = readfits(fils[i])
    if i == ref[j]:
        tform = transform.SimilarityTransform(translation=(xm[j],ym[j]))
        j = j+1
    else:
        xms = (i-ref[j-1])/(ref[j]-ref[j-1])*(xm[j]-xm[j-1])+xm[j-1]
        yms = (i-ref[j-1])/(ref[j]-ref[j-1])*(ym[j]-ym[j-1])+ym[j-1]
        tform = transform.SimilarityTransform(translation=(xms, yms))
    print(i, tform.params[0,2], tform.params[1,2], end='\r')
    data = transform.warp(data, tform, order=interpolate)
    fname = os.path.basename(fils[i])
    fits.writeto(outdir+fname, data, hdr, output_verify='fix', overwrite=True, checksum=False)
if ref[-1]+1 < len(fils):
    for i in range(ref[-1]+1, len(fils)):
        data, hdr = readfits(fils[i])
        tform = transform.SimilarityTransform(translation=(xm[-1],ym[-1]))
        print(i, tform.params[0,2], tform.params[1,2], end='\r')
        data = transform.warp(data, tform, order=interpolate)
        fname = os.path.basename(fils[i])
        mysp.fitswrite(outdir+fname, data, hdr)
