/*
       FLCA Local Correlation Tracking software
       http://cgem.ssl.berkeley.edu/~fisher/public/software/FLCA
       Copyright (C) 2007-2019, Regents of the University of California
  
       This software is based on the concepts described in Welsch & Fisher
       (2008, PASP Conf. Series 383, 373), with updates described in
       Fisher et al. 2019, "The PDFI_SS Electric Field Inversion Software",
       in prep.
       If you use the software in a scientific 
       publication, the authors would appreciate a citation to these papers 
       and any future papers describing updates to the methods.
 
       This is free software; you can redistribute it and/or
       modify it under the terms of the GNU Lesser General Public
       License version 2.1 as published by the Free Software Foundation.
 
       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
       See the GNU Lesser General Public License for more details.
 
       To view the GNU Lesser General Public License visit
       http://www.gnu.org/copyleft/lesser.html
       or write to the Free Software Foundation, Inc.,
       59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

# include <stdio.h>
# include <string.h>
# include <ctype.h>
# include <stdlib.h>
# include <math.h>

/* To include C99 complex arithmetic, uncomment the following line defining
   COMPLEXH.  To not include C99 complex arithmetic, leave this definition
   commented out. */

/* # define COMPLEXH 1 */

# ifdef COMPLEXH
   # include <complex.h>   
#endif

# include <fftw3.h> 

/* To write files deriv2.dat and deriv1.dat, containing 2nd derivatives of
the cross-correlation function, and the peak value and first derivatives,
uncomment the line below defining CCDATA: */

/* # define CCDATA 1 */

/* global declarations */

/* i4 and f4 are supposed to be definitions that give rise to 4 byte integers
 * and 4 byte floats */


typedef int i4;
typedef float f4;

/* function prototypes: */

/* flca function calling arguments not yet completely defined */

i4 flca (double * f1, double * f2, i4 nx, i4 ny, 
    double sigma, double * vx, double * vy, double * vm,
    double thresh, i4 absflag, i4 filter, double kr, i4 skip,
    i4 poffset, i4 qoffset, i4 biascor, i4 verbose);
i4 where (char *cond, i4 xsize, i4 ** index, i4 * length_index);
i4 cross_cor (i4 init, double *arr, double *barr,
	   double **absccor, i4 nx, i4 ny, double *shiftx, double *shifty, 
           i4 filterflag, double kr);
i4 writeimage (char *fname, double *arr, i4 nx, i4 ny);
i4 write2images (char *fname, double *arr, double *barr, i4 nx, i4 ny);
i4 write3images (char *fname, double *arr, double *barr, double *carr,
		 i4 nx, i4 ny);
i4 shift2d (double *arr, i4 nx, i4 ny, i4 ishift, i4 jshift);
i4 maxloc (double *arr, i4 xsize);
i4 minloc (double *arr, i4 xsize);
i4 iminloc (i4 * arr, i4 xsize);
i4 imaxloc (i4 * arr, i4 xsize);
double r (double t);
i4 interpcc2d (double *fdata, double xmiss, i4 nx, i4 ny, 
double *xwant, i4 nxinterp, double *ywant, i4 nyinterp, double **finterp);
i4 gaussfilt(double *filter, double *kx, double *ky, i4 nx, i4 ny, double kr);
i4 make_freq(double *k, i4 ndim);
i4 filter_image(double *arr, double *barr, double *outarr, double *outbarr,
        i4 nx, i4 ny, double kr);
i4 is_large_endian ();
i4 byteswapflca (unsigned char *arr, i4 arrsize, i4 nbpw);

/* end function prototypes */
