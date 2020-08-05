'''
flca.py
created by Xue, xuejc@pmo.ac.cn
History:
    2020-08-04: v0, created.
'''
import os
import time
import numpy as np
from astropy.io import fits
from glob import glob
from skimage import transform
import ctypes as ct
flcas = ct.CDLL('./codetest/flca.so')

def run(inputf):
    starttime = time.time()
    if os.path.exists(inputf) == False:
        print('No input file!')
        return
    if os.path.exists('./source/flca.so') == False:
        print('./source/flca.so is necessary!')
        return
    f=open(inputf, 'r')
    lines = f.readlines()
    f.close()
    tmp = []
    for i in range(len(lines)):
        if lines[i][0] == '#': continue
        tmp1 = lines[i].split()
        tmp.extend(tmp1)
    
    keywords = ['*indir0', '*indir','*outdir', '*firstfits', '*x0', '*x1', '*y0', '*y1', 
            '*time0a', '*time0b', '*time1a', '*time1b', '*every', '*skip',
            '*xoffset', '*yoffset', '*sigma', '*threshold', '*kr', '*biascorrect', '*interpolate', '*twochannel']
    for item in tmp:
        if item[0] == '*':
            if keywords.count(item) == 0:
                print('The keyword is unrecognized: {}'.format(item))
                print('The allowed keywords are {}'.format(keywords))
                return
    # default values
    indir0 = ''; indir=''; outdir = ''; firstfits = ''; 
    x0=300; x1=700; y0=100; y1=400; time0a=0; time0b=0; time1a=None; time0b=None;
    every=1; skip=1; xoffset=0; yoffset=0; sigma=0
    threshold='0.'; kr=0.; biascorrect=0; interpolate=0; verbose=0; twochannel=0

    if tmp.count('*indir0'):
        ind = tmp.index('*indir0')
        if tmp[ind+1][0] != '*': indir0=tmp[ind+1]
    if tmp.count('*indir'):
        ind = tmp.index('*indir')
        if tmp[ind+1][0] != '*': indir=tmp[ind+1]
    if tmp.count('*outdir'):
        ind = tmp.index('*outdir')
        if tmp[ind+1][0] != '*': outdir=tmp[ind+1]
    if tmp.count('*firstfits'):
        ind = tmp.index('*firstfits')
        if tmp[ind+1][0] != '*': firstfits=tmp[ind+1]
    if tmp.count('*x0'):
        ind = tmp.index('*x0')
        if tmp[ind+1][0] != '*': x0=int(tmp[ind+1])
    if tmp.count('*x1'):
        ind = tmp.index('*x1')
        if tmp[ind+1][0] != '*': x1=int(tmp[ind+1])
    if tmp.count('*y0'):
        ind = tmp.index('*y0')
        if tmp[ind+1][0] != '*': y0=int(tmp[ind+1])
    if tmp.count('*y1'):
        ind = tmp.index('*y1')
        if tmp[ind+1][0] != '*': y1=int(tmp[ind+1])
    if tmp.count('*time0a'):
        ind = tmp.index('*time0a')
        if tmp[ind+1][0] != '*': time0a=int(tmp[ind+1])
    if tmp.count('*time0b'):
        ind = tmp.index('*time0b')
        if tmp[ind+1][0] != '*': time0b=int(tmp[ind+1])
    if tmp.count('*time1a'):
        ind = tmp.index('*time1a')
        if tmp[ind+1][0] != '*': time1a=int(tmp[ind+1])
    if tmp.count('*time1b'):
        ind = tmp.index('*time1b')
        if tmp[ind+1][0] != '*': time1b=int(tmp[ind+1])
    if tmp.count('*every'):
        ind = tmp.index('*every')
        if tmp[ind+1][0] != '*': every=int(tmp[ind+1])
    if tmp.count('*skip'):
        ind = tmp.index('*skip')
        if tmp[ind+1][0] != '*': skip=int(tmp[ind+1])
    if tmp.count('*xoffset'):
        ind = tmp.index('*xoffset')
        if tmp[ind+1][0] != '*': xoffset=int(tmp[ind+1])
    if tmp.count('*yoffset'):
        ind = tmp.index('*yoffset')
        if tmp[ind+1][0] != '*': yoffset=int(tmp[ind+1])
    if tmp.count('*sigma'):
        ind = tmp.index('*sigma')
        if tmp[ind+1][0] != '*': sigma=int(tmp[ind+1])
    if tmp.count('*threshold'):
        ind = tmp.index('*threshold')
        if tmp[ind+1][0] != '*': threshold=tmp[ind+1]
    if tmp.count('*kr'):
        ind = tmp.index('*kr')
        if tmp[ind+1][0] != '*':
            kr=float(tmp[ind+1])
    if tmp.count('*biascorrect'): biascorrect=1
    if tmp.count('*interpolate'): interpolate=1
    #if tmp.count('*quiet'): verbose=0
    if tmp.count('*twochannel'): twochannel=1
    
    # check and adjust parameters
    if indir == '' or outdir == '':
        print("indir and outdir must be defined.")
        return
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
        print(f"{outdir} is created")
    if twochannel == 1:
        if indir0 == '':
            print('indir0 must be defined if twochannel is turned on.')
            return
        if time0a is None or time0b is None:
            print('time0a and time0b must be set if twochannel is turned on.')
            print('For example:')
            print('*time0a -23 *time0b -8')
            return
        if time1a is None:
            time1a = time0a
            time1b = time1a+(time0b-time0a)
    if sigma < 0:
        print("Sigma cannot be negative, = {}".format(sigma))
        return
    elif sigma == 0:
        print("The sigma value is suggested to be larger than 0.")
        tmp=input("Do you want to reset sigma (y or n)? ")
        if tmp.lower() == 'y': return
        elif tmp.lower() == 'n': pass
        else:
            print('unrecognized, return')
            return
    
    if twochannel == 0:
        error = onechannel(indir, outdir, firstfits, x0, x1, y0, y1, every, skip, xoffset, 
               yoffset, sigma, threshold, kr, biascorrect, interpolate, verbose)
    else:
        error = doublechannel(indir0, indir, outdir, x0, x1, y0, y1, time0a, time0b, time1a, time1b, skip, xoffset, 
               yoffset, sigma, threshold, kr, biascorrect, interpolate, verbose)
    
    if error == 0:
        endtime = time.time()
        print('Finish')
        print(f"The time cost is {round(endtime-starttime)} s.")
    else:
        print("Error!")
        
def onechannel(indir, outdir, firstfits='', x0=0, x1=-1, y0=0, y1=-1, every=1, skip=1, xoffset=0, 
               yoffset=0, sigma=20, threshold='0.', kr=0., biascorrect=0, interpolate=0, verbose=0):
    print('Single-channel coalignment.')
    absflag=0; filterflag=1; marg = 0.2 # Margins of the calculated shifts that are excluded.
    if indir[-1] != '/': indir += '/'
    if outdir[-1] != '/': outdir += '/'
    filin = glob(indir+"*.f*ts")
    if len(filin) == 0:
        print(f"The input folder{indir} is empty!")
        return -1
    filin.sort()
    ind0 = 0
    if firstfits == '': 
        firstfits = filin[0]
        ind0 = 0
    else:
        if filin.count(firstfits):
            ind0 = filin.index(firstfits)
        else:
            fname0 = os.path.basename(firstfits)
            for i in range(len(filin)):
                if os.path.basename(filin[i]) == fname0:
                    ind0 = i
                    break
            if ind0 != i:
                print("The first fits file should have the same name as one of the files in input folder!")
                return -1
            
    xoffset = xoffset % skip; yoffset = yoffset % skip
    if threshold[-1] == 'a':
        absflag = 1;
        threshold = float(threshold[0:-1])
    else:
        threshold = float(threshold)
    if kr <= 0. or kr > 20.:
        print("Nonsense value of kr, = {}".format(kr))
        return -1
    sigma = ct.c_double(sigma);
    thresh = ct.c_double(threshold);
    kr = ct.c_double(kr);
    
    # Notice that for 2-D array in python and C languages, the 1st dimension is indicated by y and the 2nd by x.
    # It is consistent with images.
    # But in flcasubs.c, the 1st dimension is marked with x and 2nd with y, following FLCT.
    ny = y1-y0
    nx = x1-x0
    nys = int(np.ceil((ny-yoffset)/skip)) # dimensions for the results
    nxs = int(np.ceil((nx-xoffset)/skip))
    indx0 = int(nxs*marg); indx1 = int(np.ceil(nxs-indx0)); indy0 = int(nys*marg); indy1 = int(np.ceil(nys-indy0))
    ArrayType = ct.c_double*(nx*ny)
    f1ct = ArrayType();
    f2ct = ArrayType();

    print("{} files will be coaligned".format(len(filin)))
    sxarr = np.zeros(len(filin)); syarr = np.zeros(len(filin))
    sxarr[ind0] = 0; syarr[ind0] = 0;
    
    if ind0 > 0:
        kk = 0 # whether change the reference fits file, correspondes to every
        data0, hdr = readfits(firstfits)
        referf = firstfits
        shiftx0 = 0; shifty0 = 0
        for ii in range(ind0, 0, -1):
            print(f"{ii-1}", end="\r")    
            data, hdr = readfits(filin[ii-1])
            tform = transform.SimilarityTransform(translation=(shiftx0,shifty0))
            data1 = transform.warp(data, tform, order=interpolate)
            f1 = data0[y0:y1, x0:x1]
            f2 = data1[y0:y1, x0:x1]
            # tate care, all the normal data are assumed to be larger than 0
            f1[np.logical_or(np.isnan(f1), np.isnan(f2))] = threshold-10
            f2[np.logical_or(np.isnan(f1), np.isnan(f2))] = threshold-10
            #f1 = detrend(f1)
            #f2 = detrend(f2)
            for i in range(ny):
                for j in range(nx):
                    f1ct[i*nx+j] = f1[i,j]
                    f2ct[i*nx+j] = f2[i,j]
            vxct = ArrayType()
            vyct = ArrayType()
            vmct = ArrayType()
            # The 1st dimension is y and the 2nd is x.
            tmp = flcas.flca(f1ct, f2ct, ny, nx, sigma, vyct, vxct, vmct, thresh, absflag, 
                             filterflag, kr, skip, yoffset, xoffset, biascorrect, verbose)
            
            # reshape and select
            vx = np.array(vxct).reshape([ny, nx])
            vy = np.array(vyct).reshape([ny, nx])
            vm = np.array(vmct).reshape([ny, nx]);
            vx = vx[yoffset::skip, xoffset::skip]
            vy = vy[yoffset::skip, xoffset::skip]
            vm = vm[yoffset::skip, xoffset::skip];
            vx = vx[indy0:indy1, indx0:indx1]
            vy = vy[indy0:indy1, indx0:indx1]
            vm = vm[indy0:indy1, indx0:indx1];
            vm[np.logical_or(np.isnan(vx), np.isnan(vy))] = 0
            shiftx = np.median(vx[vm == 1.])
            shifty = np.median(vy[vm == 1.])
            if shiftx == np.nan or shifty == np.nan:
                print("Shift X or shift Y is Nan, return.")
                return -1
            shiftx0 = shiftx+shiftx0; shifty0 = shifty+shifty0
            if interpolate == 0:
                shiftx0 = round(shiftx0); shifty0 = round(shifty0)
            tform = transform.SimilarityTransform(translation=(shiftx0,shifty0))
            data = transform.warp(data, tform, order=interpolate)
            # order = 0: Nearest-neighbor; order = 1: Bi-linear (default). Other choices are not included.
            fname0 = os.path.basename(filin[ii-1])
            fits.writeto(outdir+fname0, data, hdr, output_verify='fix', overwrite=True, checksum=False)
            sxarr[ii-1] = shiftx0; syarr[ii-1] = shifty0
            kk = kk+1
            if kk % every == 0:
                data0 = data.copy()
                referf = filin[ii-1]
    
    data0, hdr = readfits(firstfits)
    fname0 = os.path.basename(firstfits)
    referf = firstfits
    # copy the first fits file to the output folder.
    fits.writeto(outdir+fname0, data0, hdr, output_verify='fix', overwrite=True, checksum=False)
    
    if ind0 < len(filin)-1:
        kk = 0 # whether change the reference fits file, correspondes to every
        shiftx0 = 0; shifty0 = 0
        for ii in range(ind0+1, len(filin)):
            #print("{} {} {}".format(ii, os.path.basename(filin[ii]), os.path.basename(referf)), end="\n")  
            print(f"{ii}", end="\r")  
            data, hdr = readfits(filin[ii])
            tform = transform.SimilarityTransform(translation=(shiftx0,shifty0))
            data1 = transform.warp(data, tform, order=interpolate)
            f1 = data0[y0:y1, x0:x1]
            f2 = data1[y0:y1, x0:x1]
            # tate care, all the normal data are assumed to be larger than 0
            f1[np.logical_or(np.isnan(f1), np.isnan(f2))] = threshold-10
            f2[np.logical_or(np.isnan(f1), np.isnan(f2))] = threshold-10
            #f1 = detrend(f1)
            #f2 = detrend(f2)
            for i in range(ny):
                for j in range(nx):
                    f1ct[i*nx+j] = f1[i,j]
                    f2ct[i*nx+j] = f2[i,j]
            vxct = ArrayType()
            vyct = ArrayType()
            vmct = ArrayType()
            # The 1st dimension is y and the 2nd is x.
            tmp = flcas.flca(f1ct, f2ct, ny, nx, sigma, vyct, vxct, vmct, thresh, absflag, 
                             filterflag, kr, skip, yoffset, xoffset, biascorrect, verbose)
            
            # reshape and select
            vx = np.array(vxct).reshape([ny, nx])
            vy = np.array(vyct).reshape([ny, nx])
            vm = np.array(vmct).reshape([ny, nx]);
            vx = vx[yoffset::skip, xoffset::skip]
            vy = vy[yoffset::skip, xoffset::skip]
            vm = vm[yoffset::skip, xoffset::skip];
            vx = vx[indy0:indy1, indx0:indx1]
            vy = vy[indy0:indy1, indx0:indx1]
            vm = vm[indy0:indy1, indx0:indx1];
            vm[np.logical_or(np.isnan(vx), np.isnan(vy))] = 0
            shiftx = np.median(vx[vm == 1.])
            shifty = np.median(vy[vm == 1.])
            if shiftx == np.nan or shifty == np.nan:
                print("Shift X or shift Y is Nan, return.")
                return -1
            shiftx0 = shiftx+shiftx0; shifty0 = shifty+shifty0
            if interpolate == 0:
                shiftx0 = round(shiftx0); shifty0 = round(shifty0)
            tform = transform.SimilarityTransform(translation=(shiftx0,shifty0))
            data = transform.warp(data, tform, order=interpolate)
            # order = 0: Nearest-neighbor; order = 1: Bi-linear (default).
            # Other choices are not included.
            fname0 = os.path.basename(filin[ii])
            fits.writeto(outdir+fname0, data, hdr, output_verify='fix', overwrite=True, checksum=False)
            sxarr[ii] = shiftx0; syarr[ii] = shifty0
            kk = kk+1
            if kk % every == 0:
                data0 = data.copy()
                #referf = filin[ii]
    np.savez_compressed('shiftxy.npz', sx=sxarr, sy=syarr)
    return 0
    
def doublechannel(indir0, indir, outdir, x0=0, x1=-1, y0=0, y1=-1, time0a=-23, time0b=-8, time1a=None,
               time1b=None, skip=1, xoffset=0, yoffset=0, sigma=20, threshold='0.', kr=0., 
               biascorrect=0, interpolate=0, verbose=0):
    print('Double-channel coalignment.')
    import datetime
    def is_digit(stri):
        return stri.isdigit()
    absflag=0; filterflag=1; marg = 0.2 # Margins of the calculated shifts that are excluded.
    if indir0[-1] != '/': indir0 += '/'
    if indir[-1] != '/': indir += '/'
    if outdir[-1] != '/': outdir += '/'
    filin0 = glob(indir0+"*.f*ts")
    nfilin0 = len(filin0)
    if nfilin0 == 0:
        print(f"The reference folder {indir0} is empty!")
        return -1
    filin0.sort()
    filin = glob(indir+"*.f*ts")
    nfilin = len(filin)
    if nfilin == 0:
        print(f"The input folder {indir} is empty!")
        return -1
    filin.sort()
    #ind0 = 0
    
    tims0 = "".join(filter(is_digit, filin0[0][time0a:time0b]))
    tims1 = "".join(filter(is_digit, filin[0][time1a:time1b]))
    if len(tims0) != 14 or len(tims1) != 14:
        print("time0a, time0b, time1a, time1b: time should be from year to second")
        return -1
    timarr = [];
    for ii in range(nfilin0):
        tims = "".join(filter(is_digit, filin0[ii][time0a:time0b]))
        timdate=datetime.datetime(int(tims[0:4]), int(tims[4:6]), int(tims[6:8]), int(tims[8:10]), 
                                  int(tims[10:12]), int(tims[12:14]))
        timarr.append(timdate)
    timarr = np.array(timarr, dtype=datetime.datetime)
            
    xoffset = xoffset % skip; yoffset = yoffset % skip
    if threshold[-1] == 'a':
        absflag = 1;
        threshold = float(threshold[0:-1])
    else:
        threshold = float(threshold)
    if kr <= 0. or kr > 20.:
        print("Nonsense value of kr, = {}".format(kr))
        return -1
    sigma = ct.c_double(sigma);
    thresh = ct.c_double(threshold);
    kr = ct.c_double(kr);
    
    # Notice that for 2-D array in python and C languages, the 1st dimension is indicated by y and the 2nd by x.
    # It is consistent with images.
    # But in flcasubs.c, the 1st dimension is marked with x and 2nd with y, following FLCT.
    ny = y1-y0
    nx = x1-x0
    nys = int(np.ceil((ny-yoffset)/skip)) # dimensions for the results
    nxs = int(np.ceil((nx-xoffset)/skip))
    indx0 = int(nxs*marg); indx1 = int(np.ceil(nxs-indx0)); indy0 = int(nys*marg); indy1 = int(np.ceil(nys-indy0))
    ArrayType = ct.c_double*(nx*ny)
    f1ct = ArrayType();
    f2ct = ArrayType();

    print("{} files will be coaligned".format(nfilin))
    
    shiftx0 = 0; shifty0 = 0
    sxarr = np.zeros(nfilin); syarr = np.zeros(nfilin)
    i0 = 0
    for ii in range(nfilin):
        tims1 = "".join(filter(is_digit, filin[ii][time1a:time1b]))
        timdate=datetime.datetime(int(tims1[0:4]), int(tims1[4:6]), int(tims1[6:8]), int(tims1[8:10]), 
                                  int(tims1[10:12]), int(tims1[12:14]))
        for i1 in range(i0, nfilin0):
            if timarr[i1]>timdate: break
        if i1 > i0:
            timd1 = timdate-timarr[i1-1]
            timd2 = timarr[i1]-timdate
            if timd1.seconds < timd2.seconds:
                i0 = i1-1
            else:
                i0 = i1
        print("{} {} {}".format(ii, os.path.basename(filin[ii]), os.path.basename(filin0[i0])), end="\n")    
        data, hdr = readfits(filin0[i0])
        f1 = data[y0:y1, x0:x1]
        data, hdr = readfits(filin[ii])
        tform = transform.SimilarityTransform(translation=(shiftx0,shifty0))
        data1 = transform.warp(data, tform, order=interpolate)
        f2 = data1[y0:y1, x0:x1]
        # tate care, all the normal data are assumed to be larger than 0
        f1[np.logical_or(np.isnan(f1), np.isnan(f2))] = threshold-10
        f2[np.logical_or(np.isnan(f1), np.isnan(f2))] = threshold-10
        #f1 = detrend(f1)
        #f2 = detrend(f2)
        for i in range(ny):
            for j in range(nx):
                f1ct[i*nx+j] = f1[i,j]
                f2ct[i*nx+j] = f2[i,j]
        vxct = ArrayType()
        vyct = ArrayType()
        vmct = ArrayType()
        # The 1st dimension is y and the 2nd is x.
        tmp = flcas.flca(f1ct, f2ct, ny, nx, sigma, vyct, vxct, vmct, thresh, absflag, 
                         filterflag, kr, skip, yoffset, xoffset, biascorrect, verbose)
        
        # reshape and select
        vx = np.array(vxct).reshape([ny, nx])
        vy = np.array(vyct).reshape([ny, nx])
        vm = np.array(vmct).reshape([ny, nx]);
        vx = vx[yoffset::skip, xoffset::skip]
        vy = vy[yoffset::skip, xoffset::skip]
        vm = vm[yoffset::skip, xoffset::skip];
        vx = vx[indy0:indy1, indx0:indx1]
        vy = vy[indy0:indy1, indx0:indx1]
        vm = vm[indy0:indy1, indx0:indx1];
        vm[np.logical_or(np.isnan(vx), np.isnan(vy))] = 0
        shiftx = np.median(vx[vm == 1.])
        shifty = np.median(vy[vm == 1.])
        if shiftx == np.nan or shifty == np.nan:
            print("Shift X or shift Y is Nan, return.")
            return -1
        shiftx0 = shiftx+shiftx0; shifty0 = shifty+shifty0
        if interpolate == 0:
            shiftx0 = round(shiftx0); shifty0 = round(shifty0)
        tform = transform.SimilarityTransform(translation=(shiftx0,shifty0))
        data = transform.warp(data, tform, order=interpolate)
        # order = 0: Nearest-neighbor; order = 1: Bi-linear (default). Other choices are not included.
        fname0 = os.path.basename(filin[ii])
        fits.writeto(outdir+fname0, data, hdr, output_verify='fix', overwrite=True, checksum=False)
        sxarr[ii] = shiftx0; syarr[ii] = shifty0
    np.savez_compressed('shiftxy.npz', sx=sxarr, sy=syarr)
    return 0

def readfits(f):
    hdul = fits.open(f)
    data = hdul[0].data
    hdr = hdul[0].header
    hdul.close()
    return data, hdr

# flca.twoimg
# Purpoes:
#     Comparing two images and get shifts of each pixel.
# return:
#     vx, shifts in X direction
#     vy, shifts in Y direction
#     vm, mesh of vx and vy. If the image intensity is higher than threshold, vm=1;
#         else vm=0 and corresponding vx and vy values shouldn't be used.
# Example: vx, vy, vm = flca.twoimg(f1, f2, sigma=20, threshold=0.3, biascorrect=1)
def twoimg(f1, f2, sigma=0, threshold=0, absflag=1, filterflag=1, kr=0.4, skip=0, xoff=0, yoff=0, biascorrect=0, verbose=0):
    if f1.shape != f2.shape:
        print("Two images must have the same size.")
        return
    ny, nx = f1.shape
    xoff = xoff % skip; yoff = yoff % skip
    sigmact = ct.c_double(sigma);
    thresh = ct.c_double(threshold);
    kr = ct.c_double(kr);
    
    ArrayType = ct.c_double*(nx*ny)
    f1ct = ArrayType();
    f2ct = ArrayType();
    vxct = ArrayType();
    vyct = ArrayType();
    vmct = ArrayType();

    for i in range(ny):
        for j in range(nx):
            f1ct[i*nx+j] = f1[i,j]
            f2ct[i*nx+j] = f2[i,j]
    tmp = flcas.flca(f1ct, f2ct, ny, nx, sigmact, vyct, vxct, vmct, thresh, absflag, filterflag, kr, skip, yoff, xoff, biascorrect, verbose)
    
    vx = np.array(vxct).reshape([ny, nx]); vy = np.array(vyct).reshape([ny, nx]); vm = np.array(vmct).reshape([ny, nx])
    vx = vx[yoff::skip, xoff::skip]
    vy = vy[yoff::skip, xoff::skip]
    vm = vm[yoff::skip, xoff::skip]
    
    return vx, vy, vm
