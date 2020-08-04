/*
       FLCA Local Correlation Alignment software
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

# include "flcasubs.h"

i4 flca (double * f1, double * f2, i4 nx, i4 ny, 
    double sigma, double * vx, double * vy, double * vm,
    double thresh, i4 absflag, i4 filter, double kr, i4 skip,
    i4 poffset, i4 qoffset, i4 biascor, i4 verbose) 
{

/* BEGIN FLCA FUNCTION */

/*  char *version ="1.07    "; */

  i4 sigmaeq0;
  i4 nxorig, nyorig, i, j, ii, jj;
  i4 icc, nsize=0, nt, ndmin, init;
  i4 iloc1, iloc2, imin0, jmax0, jmin0, imax0, imin, jmin, imax, jmax,
    isize, jsize;
  i4 skipon=0, noskipx,noskipy,noskipxy;
  i4 xoffset=0,yoffset=0 ;
  double *gaussdata, *f1temp, *f2temp,
    *g1, *g2, *absccor;
  double shiftx, shifty, argx, argy, f1max, f2max, fmax, fabsbar, vmask;
  double sigminv=-999.; 
  double f1bar=0., f2bar=0.;
/*      double tol=1e-4; */
  double tol = 1e-2;
  i4 hardworkneeded, belowthresh;
  i4 maxind,ixmax,iymax,absccmax,nsig;
  double fx,fy,fxx,fyy,fxy,fpeak,hessian,hm1over2,corfac,gam2oversigma2;

/* nx, ny may get set to 1 if sigma=0. so copy original values */

  nxorig=nx;
  nyorig=ny;

/* Get information out of sigma: */

  if(sigma > 0.)
  {
     sigminv = 1. / sigma;
     sigmaeq0=0;
  }
  else
  {
     sigmaeq0 = 1;
  }

  if(!sigmaeq0)
  {
     /* This stuff done only if sigma > 0 */
     nt = (i4) sigma *sqrt (log (1. / tol));
     ndmin = (nx < ny) ? (((nx / 3) / 2)) * 2 : (((ny / 3) / 2)) * 2;
     nsize = ((nt / 2) * 2 < ndmin) ? (nt / 2) * 2 : ndmin;
     if (verbose) 
     {
       printf ("flca: nominal sliding box size = %d\n", 2 * nsize);
       fflush(stdout);
     }
     if(nsize <= 0)
     {
        printf("flca: error - illegal box size, exiting\n");
        fflush(stdout);
        exit(1);
     }
  }
  if(sigmaeq0) 
  {
     /* sigma = 0 means we'll only compute one point */
     nx=1;
     ny=1;
  }

  /* figure out if threshold is in absolute or fractional units
   * and if fractional, convert to absolute units.  If thresh is between
   * zero and 1 (not inclusive) then it's assumed to be fractional,
   * (unless the threshold string ends with 'a') and must be scaled.
   * if aloc == NULL, there's no 'a' in the threshold string. */

  if(!sigmaeq0)
  {
    if ((thresh > 0.) && (thresh < 1.) && (absflag == 0))
      {
        f1temp = (double *) malloc (sizeof (double) * nx * ny);
        f2temp = (double *) malloc (sizeof (double) * nx * ny);

        for (i = 0; i < nx * ny; i++)

  	{

	  /* compute abs value of f1,f2 arrays as f1temp,
	   * f2temp arrays */

	  *(f1temp + i) = (double) fabs (*(f1 + i));
	  *(f2temp + i) = (double) fabs (*(f2 + i));
	}

      /* now find maximum absolute value of both images */

        iloc1 = maxloc (f1temp, nx * ny);
        iloc2 = maxloc (f2temp, nx * ny);
        f1max = *(f1temp + iloc1);
        f2max = *(f2temp + iloc2);
        fmax = (f1max > f2max) ? f1max : f2max;

      /* now convert relative thresh to absolute threshhold */

        thresh *= fmax;
        if (verbose) 
           printf ("flca: relative threshold in abs. units = %g\n", thresh);
           fflush(stdout);

        free (f1temp);
        free (f2temp);
    }
  }

  /* debug: output the two input images to file "f1f2.dat" */
/*
      write2images("f1f2.dat",f1,f2,nx,ny,transp);
*/

  /* the vm array (velocity mask array, not to be confused with the
   * gaussian mask array that is used to modulate the images) 
   * will later be set to 1.0 where the velocity is computed,
   * and to 0.0 where the velocity is not computed -- because the image
   * pixels are below the noise threshold value "thresh" input from
   * the command line. */

  /* Now create master gaussian image mask: */

  gaussdata = (double *) malloc (sizeof (double) * (2 * nxorig) * (2 * nyorig));

  if(!sigmaeq0) /* this case for sigma > 0 */
  {
    for (i = 0; i < 2 * nxorig; i++)
    {
        argx = sigminv * (double) (i - nxorig);
        for (j = 0; j < 2 * nyorig; j++)
        {
	  argy = sigminv * (double) (j - nyorig);
	  *(gaussdata + i * (2 * ny) + j) = exp (-argx * argx - argy * argy);
        }
    }
  }
  else /* this case for sigma = 0. ie set gaussian to 1.0 */
  {
    for (i = 0; i < 2 * nxorig; i++)
    {
        for (j = 0; j < 2 * nyorig; j++)
        {
	  *(gaussdata + i * (2 * nyorig) + j) = (double) 1.;
        }
    }
  }

  gam2oversigma2=0.;
  nsig=0;

  /* Now do the master loop over i,j for computing velocity field: */

  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
	{
          if((nx ==1) && (ny == 1))
            {
              init=2; /* special case: must initialize AND destroy plans */
            }
	  else if ((i == 0) && (j == 0) && ((nx+ny) > 2))
	    {
	      /* 1st time through, set init to 1 so that
	       * fftw FFT plans are initialized */

	      init = 1;
	    }
	  else if ((i == (nx - 1)) && (j == (ny - 1)) && ((nx+ny) > 2))
	    {
	      /* last time through, set init to -1 so that
	       * fftw static variables are freed and plans
	       * destroyed */

	      init = -1;
	    }
	  else
	    {
	      /* the rest of the time just chunk along */

	      init = 0;
	    }

	  /* Now, figure out if image value is below
	   * threshold: */

	  /* the data is considered below theshhold if the
	   * absolute value of average of the pixel value from the 2 images 
	   * is below threshold */

	  fabsbar = 0.5 * (fabs (*(f1 + i * ny + j) + *(f2 + i * ny + j)));
	  belowthresh = (fabsbar < thresh);

	  /* Or alternatively:
	     belowthresh = ((fabs1 < thresh) || (fabs2 < thresh));
	   */

	  /* all the hard work of doing the cross-correlation
	   * needs to be done if the avg data is above the
	   * threshold OR if init != 0 */

          /* added skip logic here */

          skipon=skip+abs(qoffset)+abs(poffset);

          if(skipon)
          {
             xoffset=poffset;
             yoffset=qoffset;
             noskipx = !((i-xoffset) % skip);
             noskipy = !((j-yoffset) % skip);
             noskipxy=noskipx*noskipy;
          }
          else
          {
             noskipxy=1;
          }

	  hardworkneeded = (((!belowthresh) && (noskipxy)) || (init != 0));

	  if (hardworkneeded)
	    {

	      /* the hard work for this particular pixel starts
	       * now */


	      /* Now find where the gaussian modulated image 
	       * is
	       * chopped off by the sliding box.  The sliding
	       * box is centered at i,j, unless the edges of
	       * the box would go outside the array -- 
	       * then the
	       * sliding box just sits at edges and/or corners
	       * of the array   */

              if(!sigmaeq0) /* for sigma > 0 */
              {
	        imin0 = (0 > (i - (nsize - 1))) ? 0 : i - (nsize - 1);
                imax0 = ((nx - 1) < (i + nsize)) ? nx - 1 : i + nsize;
                imin = (imax0 == nx - 1) ? nx - 1 - (2 * nsize - 1) : imin0;
                imax = (imin0 == 0) ? 0 + (2 * nsize - 1) : imax0;

                jmin0 = (0 > (j - (nsize - 1))) ? 0 : j - (nsize - 1);
                jmax0 = ((ny - 1) < (j + nsize)) ? ny - 1 : j + nsize;
                jmin = (jmax0 == ny - 1) ? ny - 1 - (2 * nsize - 1) : jmin0;
                jmax = (jmin0 == 0) ? 0 + (2 * nsize - 1) : jmax0;

                isize = imax - imin + 1;
                jsize = jmax - jmin + 1;

                 /* If the following tests for isize, jsize fail,
	         this is very bad:  exit */

	        if (isize != 2 * nsize)
                {
		  printf ("flca: exiting, bad isize = %d\n", isize);
		  exit (1);
                }
	        if (jsize != 2 * nsize)
                {
		  printf ("flca: exiting, bad jsize = %d\n", jsize);
		  exit (1);
                }
              }
              else /* if sigma = 0. just set isize=nxorig, jsize=nyorig */
              {
                 isize=nxorig;
                 jsize=nyorig;
                 imin=0;
                 jmin=0;
              }
              /* debug:
              printf("isize = %d, jsize = %d,\n",isize,jsize);
              */

              /* Compute sub-image means of f1 and f2: */

              f1bar=0.;
              f2bar=0.;
              for (ii = 0; ii < isize; ii++)
                { 
                   for (jj = 0; jj < jsize; jj++)
                      {
                         f1bar=f1bar+ *(f1 + (ii+imin)*nyorig + (jj+jmin));
                         f2bar=f2bar+ *(f2 + (ii+imin)*nyorig + (jj+jmin));
                      }
                }

              f1bar=f1bar/((double)isize*jsize);
              f2bar=f2bar/((double)isize*jsize);

	      g1 = (double *) malloc (sizeof (double) * isize * jsize);
	      g2 = (double *) malloc (sizeof (double) * isize * jsize);

	      /* Now fill the reduced size arrays (sub-images) with the 
	       * appropriate values from the 
	       * full-sized arrays: */

	      for (ii = 0; ii < isize; ii++)
		{
		  for (jj = 0; jj < jsize; jj++)
		    {
		      *(g1 + ii * jsize + jj) = 
                          *(gaussdata + (nxorig-i+(ii+imin))*2*nyorig
                           +nyorig-j+(jj+jmin)) *
                          (*(f1 + (ii + imin) * nyorig + (jj + jmin))-f1bar) ;

		      *(g2 + ii * jsize + jj) = 
                          *(gaussdata + (nxorig-i+(ii+imin))*2*nyorig
                           +nyorig-j+(jj+jmin)) *
                          (*(f2 + (ii + imin) * nyorig + (jj + jmin))-f2bar) ;
		    }
		}

	      /* Call to cross_cor is used to find the 
	       * relative
	       * shift of image g2 relative to g1: */

	      icc = cross_cor (init, g1, g2, &absccor,
			       isize, jsize, &shiftx, &shifty, filter, 
                               kr);

/*                              debug:  output of absccor */

/*
                              writeimage("absccor.dat",absccor,isize,jsize,
                                 transp);
*/


              absccmax=1;
              maxind = maxloc (absccor, isize * jsize);
              if( *(absccor+maxind) == (double)0.)
              {
                 absccmax=0;
              }
              if(absccmax == 1)
              {
                 ixmax = maxind / jsize;
                 iymax = maxind % jsize;
              }
              else
              {
                 ixmax = isize/2;
                 iymax = jsize/2;
              }
              /* printf("flca: ixmax = %d, iymax = %d\n",ixmax,iymax); */
              if( (ixmax > 0) && (ixmax < (isize-1)) && (iymax > 0) && 
                  (iymax < (jsize-1)) && (absccmax == 1))
              {

                fx=0.5* ( *(absccor+(ixmax+1)*jsize+iymax) - 
                     *(absccor+(ixmax-1)*jsize+iymax) );
                fy=0.5* ( *(absccor+ixmax*jsize+iymax+1) - 
                   *(absccor+ixmax*jsize+iymax-1) );
                fxx = ( *(absccor+(ixmax+1)*jsize+iymax)+ 
                   *(absccor+(ixmax-1)*jsize+iymax) 
                  -2.*( *(absccor+ixmax*jsize+iymax))  );
                fyy = ( *(absccor+ixmax*jsize+iymax+1) + 
                   *(absccor+ixmax*jsize+iymax-1)
                  -2.*( *(absccor+ixmax*jsize+iymax)) );
                fxy = 0.25*( *(absccor+(ixmax+1)*jsize+iymax+1) + 
                   *(absccor+(ixmax-1)*jsize+iymax-1) -
                   *(absccor+(ixmax+1)*jsize+iymax-1) - 
                   *(absccor+(ixmax-1)*jsize+iymax+1) );
                fpeak=*(absccor+ixmax*jsize+iymax);

                hessian=fxx*fyy/(fpeak*fpeak)-fxy*fxy/(fpeak*fpeak);
                hm1over2=1./sqrt(hessian);
                /* Don't let ratio of hm1over2 to sigma^2 approach 1 
                 or get singularity in corfac */
                if((hm1over2 > 0.95*sigma*sigma) || (hessian == 0))
                {
                  hm1over2=0.95*sigma*sigma;
                }
                if(!sigmaeq0 && (hessian > 0))
                {
                  /* develop statistics for mean value of gamma^2/sigma^2 */
                  gam2oversigma2+=(hm1over2/(sigma*sigma));
                  nsig++;
                }
                if( (sigmaeq0) || (!biascor)) 
                {
                  corfac=1.;
                }
                else 
                {
                  corfac=1./(1.-0.8*hm1over2/(sigma*sigma));
                }
                /* Now add corrections to shiftx, shifty 
                based on bias correction factor */
                shiftx*=corfac;
                shifty*=corfac;
              }

/* free temporary arrays created during loop */

	      free (g1);
	      free (g2);
	      free (absccor);

	      /* all the hard work for this pixel is now done */

	    }

	  /* default value for vmask is 1. */

	  vmask = 1.;

	  if ((belowthresh || !noskipxy) && !sigmaeq0)

	    /* If data below threshold, set shiftx, shifty to 0. 
	     * and vmask to 0, meaning vel not computed. */
            /* If sigma=0 just ignore the threshold and compute anyway */

	    {
	      shiftx = 0.;
	      shifty = 0.;
	      vmask = 0.;
	    }


	  *(vx + i * ny + j) = shiftx;
	  *(vy + i * ny + j) = shifty;
	  *(vm + i * ny + j) = vmask;
	}
    }
   if(!sigmaeq0)
   {
/*   printf("flca debug: gam2oversigma2 = %g, nsig=%d\n",gam2oversigma2,nsig);*/
     gam2oversigma2/=nsig;
     if(verbose)
     {
       printf ("flca: mean value of gamma^2/sigma^2 = %g\n",gam2oversigma2);
       fflush(stdout);
     }
   }

  /* Outer loops over i,j finally done! */
  /* Output the vx, vy arrays to the output file 'outfile': */

  /* Output step now must be done in main program */
  /* write3images (outfile, vx, vy, vm, nx, ny, transp); */

  /* free the gaussian mask array, the original images, and the
   * velocity arrays */

  free (gaussdata);

/* Don't want to free these things anymore, they get returned */
/*
  free (f1);
  free (f2);
  free (vx);
  free (vy);
  free (vm);
*/

  if (verbose)
  {
    printf ("\nflca: finished\n");
    fflush(stdout);
  }

  /* we're done! */
  return 0;
  /*  END FLCA FUNCTION */
}

i4 where (char *cond, i4 xsize, i4 ** index, i4 * length_index)
/* This function serves a similar purpose to the "where" function in IDL.
 * Given the array *cond (char) of length xsize 
 * containing a pre-computed condition, the function finds **index, a double
 * pointer to an array of indices which reference those values
 * where *cond is
 * non-zero.  The integer *length_index returns the length of *index. 
 */
/* Note - this function no longer used in vel_ccor */
{
  i4 ier;	/* function return value - not thought through yet */
  i4 i, ii;	/* counter variables */
  i4 *indtmp;	/* temporary local array of indices of *x */
  ier = 0;	/*return value of function */
  *length_index = 0;	/* initialize length of *index array to 0 */

/*	printf("\nxsize = %d",xsize); */

  indtmp = (i4 *) malloc (sizeof (i4) * xsize);	/* create temp. ind. array */

  /* Ready to start */

  ii = 0;		/* set initial counter of temp index array to 0 */
  for (i = 0; i < xsize; i++)	/* start incrementing the *cond array: */
    {
      if (*(cond + i))
	{
	  /* if the condition is true, record i into temp index */
	  *(indtmp + ii) = (i4) i;
	  ii++;		/* and then increment ii */
	}
      /* otherwise just keep incrementing i and doing nothing */
    }

/*	printf("\nii= %d\n", ii) ;
	fflush (stdout) ; */

  /* Now create right amount of space for *index: */

  *index = (i4 *) malloc (sizeof (i4) * ii);

  /* Now copy index values from temp array into *index array */

  memcpy ((void *) *index, (void *) indtmp, ii * sizeof (i4));

  /* Now set the length of the *index array */

  *length_index = (i4) ii;

  /* Now free memory from temp. index array */

  free (indtmp);
  return ier;			/* always 0 at the moment */
}

i4 cross_cor (i4 init, double *arr, double *barr,
	   double **absccor, i4 nx, i4 ny, double *shiftx, double *shifty, 
           i4 filterflag, double kr)
{
  i4 i, j, maxind, ixmax, iymax, ishft, absccmax;
  double normfac, shiftx0, shifty0, shiftxx, shiftyy;
  double shiftsubx, shiftsuby, fx, fy, fxx, fyy, fxy;

  /* following variables must be saved between calls; declared static */

  static double *ina, *inb, *ccor;
  static double *filter, *kx, *ky;
  static fftw_complex *outa, *outb, *ccorconj;
  static fftw_plan pa, pb, pback;

  /* absccor is a double pointer containing abs. value of cc function */

  /* debug:
  printf("cross_cor: nx = %d, ny = %d\n",nx,ny);
  */

  *absccor = malloc (sizeof (double) * nx * ny);

/*	printf("initialization stuff done in cross_cor\n"); */
  if ((init == 1) || (init == 2))
    {
      /* First time through: */
      /* Initialization of FFT variables and FFTW plans. */
      /* NOTE -- empirically had to add 1 to "y" dimensions of outa,
       * outb, and ccorconj to
       * avoid a memory leak and a seg fault at fftw_free */

      /* should check to see if still a problem */

      outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
					   nx * ((ny / 2) + 2));
      outb = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * 
              nx * ((ny / 2) + 2));	/* valgrind sometimes complains */
      ccorconj = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
					       nx * ((ny / 2) + 2));

      ina = (double *) fftw_malloc (sizeof (double) * nx * ny);
      inb = (double *) fftw_malloc (sizeof (double) * nx * ny);
      ccor = (double *) fftw_malloc (sizeof (double) * nx * ny);
      filter = (double *) fftw_malloc (sizeof (double) * nx * ny);
      kx=(double *) fftw_malloc (sizeof(double)*nx);
      ky=(double *) fftw_malloc (sizeof(double)*ny);
      if(filterflag)
      {
         make_freq(kx,nx);
         make_freq(ky,ny);
         gaussfilt(filter,kx,ky,nx,ny,kr);
      }

      for (i = 0; i < nx * ny; i++)
	{
	  *(ina + i) = (double) 0.;
	  *(inb + i) = (double) 0.;
	}
      for (i = 0; i < nx * ((ny / 2 + 1)); i++)
      {
#ifdef COMPLEXH
        /* If complex.h included, do this: */
        *(ccorconj+i)=0.+I*0.;
#else
        /* If complex.h not included, do this: */
	ccorconj[i][0] = 0.;
	ccorconj[i][1] = 0.;
#endif
      }
      pa = fftw_plan_dft_r2c_2d (nx, ny, ina, outa, FFTW_MEASURE);
      pb = fftw_plan_dft_r2c_2d (nx, ny, inb, outb, FFTW_MEASURE);
      pback = fftw_plan_dft_c2r_2d (nx, ny, ccorconj, ccor, FFTW_MEASURE);
    }

    for (i = 0; i < nx * ny; i++)

    {
/*		printf("1st loop: i = %d, *(arr+i)= %g, *(barr+i) = %g\n",
				i,*(arr+i),*(barr+i)); */

      /* copy from input doubles to fftw variables */

      *(ina + i) = (double) (*(arr + i));
      *(inb + i) = (double) (*(barr + i));
    }

  /* actually do the forward FFTs: */

  fftw_execute (pa);
  fftw_execute (pb);

  /* calculate normalization factor */

  normfac = (1. / ((double) nx * ny));
  normfac *= normfac; /* square of above line */
 
/* Now apply the gaussian filter to the FFT'd data in frequency space */

    if(filterflag)
    {
      for (i=0;i<nx;i++)
      {
         for (j=0;j<(ny/2)+1;j++)
         {
#ifdef COMPLEXH
/* If complex.h is invoked, just multiply outa and outb by filter function */
            outa[i*((ny/2)+1)+j]=outa[i*((ny/2)+1)+j]*filter[i*ny+j];
            outb[i*((ny/2)+1)+j]=outb[i*((ny/2)+1)+j]*filter[i*ny+j];
#else
/* If complex.h not invoked, multiply real and im parts by filter function: */
            outa[i*((ny/2)+1)+j][0]=outa[i*((ny/2)+1)+j][0]*filter[i*ny+j];
            outa[i*((ny/2)+1)+j][1]=outa[i*((ny/2)+1)+j][1]*filter[i*ny+j];
            outb[i*((ny/2)+1)+j][0]=outb[i*((ny/2)+1)+j][0]*filter[i*ny+j];
            outb[i*((ny/2)+1)+j][1]=outb[i*((ny/2)+1)+j][1]*filter[i*ny+j];
#endif            
         }
      }
    }

  /* Now calculate product of conj(outa) * outb */

  for (i = 0; i < nx * ((ny/2) + 1); i++)
    {

#ifdef COMPLEXH
        /* if complex.h included, do this */
        *(ccorconj+i)=(conj(*(outa+i))*(*(outb+i)))*normfac;
#else
        /* if complex.h not invoked, do this: */
        ccorconj[i][0] = (outa[i][0] * outb[i][0] + outa[i][1] * outb[i][1])
	* normfac;
        ccorconj[i][1] = (outa[i][0] * outb[i][1] - outa[i][1] * outb[i][0])
	* normfac;
#endif
    }

  /* now do the inverse transform to get cc function */

  fftw_execute (pback);

  /* now calculate the absolute value of cc function */

  for (i = 0; i < nx * ny; i++)
    {
      *(*absccor + i) = (double) fabs(*(ccor+i));
    }

  if ((init == -1) || (init == 2))
    {
      /* Last time through: free all the plans and static variables */

      fftw_free (outa);
      fftw_free (outb);
      fftw_free (ccorconj);
      fftw_free (ccor);
      fftw_free (filter);
      fftw_free (kx);
      fftw_free (ky);
      fftw_free (ina);
      fftw_free (inb);
      fftw_destroy_plan (pa);
      fftw_destroy_plan (pback);
      fftw_destroy_plan (pb);
    }

/* Now shift the absccor array by nx/2, ny/2 to avoid aliasing problems */
/* ishft is set to 0, unused, and absccor is shifted. Xue */

  ishft = shift2d (*absccor, nx, ny, nx / 2, ny / 2);

  /* Now find maximum of the shifted cross-correlation function to 1 pixel
     accuracy:  */

  absccmax=1;
  maxind = maxloc (*absccor, nx * ny);
  if( *(*absccor+maxind) == (double)0.) 
  {
     absccmax=0;
  }
  if(absccmax == 1)
  {
     ixmax = maxind / ny;
     iymax = maxind % ny;
  }
  else
  {
     ixmax = nx/2;
     iymax = ny/2;
  }
  shiftx0 = ixmax;
  shifty0 = iymax;
  shiftsubx=0.;
  shiftsuby=0.;

  if((ixmax > 0) && (ixmax < (nx-1))
     && (iymax > 0) && (iymax < (ny-1)) && (absccmax == 1))
  {
     fx=0.5* ( *(*absccor+(ixmax+1)*ny+iymax) - 
         *(*absccor+(ixmax-1)*ny+iymax) );
     fy=0.5* ( *(*absccor+ixmax*ny+iymax+1) - *(*absccor+ixmax*ny+iymax-1) );
     fxx = ( *(*absccor+(ixmax+1)*ny+iymax)+ *(*absccor+(ixmax-1)*ny+iymax)
        -2.*( *(*absccor+ixmax*ny+iymax))  );
     fyy = ( *(*absccor+ixmax*ny+iymax+1) + *(*absccor+ixmax*ny+iymax-1)
        -2.*( *(*absccor+ixmax*ny+iymax)) );
     fxy = 0.25*( *(*absccor+(ixmax+1)*ny+iymax+1) + 
            *(*absccor+(ixmax-1)*ny+iymax-1) -
            *(*absccor+(ixmax+1)*ny+iymax-1) - 
            *(*absccor+(ixmax-1)*ny+iymax+1) );
/* In following expressions for subshifts, shift is in units of pixel length */
     shiftsubx=(fyy*fx-fy*fxy)/(fxy*fxy-fxx*fyy);
     shiftsuby=(fxx*fy-fx*fxy)/(fxy*fxy-fxx*fyy);
  }

  shiftxx=shiftx0 + shiftsubx;
  shiftyy=shifty0 + shiftsuby;
/*
       printf("shiftx0-nx/2 = %g\n",(shiftx0-(double)(nx/2)));
       printf("shifty0-ny/2 = %g\n",(shifty0-(double)(ny/2)));
 
*/

/* Now, assign values to shiftx, shifty to return to calling program */

  *shiftx = shiftxx - (double) (nx / 2);
  *shifty = shiftyy - (double) (ny / 2);

/* Following expressions used if only 1 pixel accuracy needed from absccor 
 *
	*shiftx=((double)ixmax)-(double)(nx/2);
	*shifty=((double)iymax)-(double)(ny/2);
*/

  return 0;
}

i4 make_freq(double *k, i4 ndim)
{
/* k is assumed already allocated in main program, with dimension ndim */
i4 n21,i,inext;
n21=(ndim/2)-1;
for (i=0;i<n21+1;i++)
{
	k[i]=(double)i;
}

inext=n21+1;
if((ndim/2)*2 != ndim)
{
	k[inext]=(double)(ndim/2);
	inext++;
        k[inext]=-(double)(ndim/2);
        inext++;
}

else
{
	k[inext]=(double)(ndim/2);
        inext++;
}

for (i=inext;i<ndim;i++)
{
	k[i]=-(double)(n21-(i-inext));
}
/* debug */

/*
for (i=0;i<ndim;i++)
{
	printf("i = %d, k = %g\n",i,k[i]);
}
*/

/* end debug */

return 0;
}

i4 gaussfilt(double *filter, double *kx, double *ky, i4 nx, i4 ny, double kr)
{
/* Assumes kx of size nx, ky of size ny, and filter of size (nx,ny) */
double kxmax,kymax,kxroll,kyroll,smxinv,smyinv,argx,argy;
i4 i,j;
kxmax=(double)kx[nx/2];
kymax=(double)ky[ny/2];
kxroll=kr*kxmax;
kyroll=kr*kymax;
smxinv=(double)1./kxroll;
smyinv=(double)1./kyroll;
for (i=0;i<nx;i++)
{
	argx=kx[i]*smxinv;
	for(j=0;j<ny;j++)
	{
                argy=ky[j]*smyinv;
		filter[i*ny+j]=exp( -(argx*argx + argy*argy) );
	}
}
return 0;
}

i4 filter_image(double *arr, double *barr, double *outarr, double *outbarr,
        i4 nx, i4 ny, double kr)

/* Takes images arr, barr and filters them by gaussians in k-space of width
kr*kmax, where 0 < kr < 1, and kmax is the maximum wavenumber in x and y,
considered separately.  The input arrays are arr, barr, and the output
arrays are outarr, and outbarr.  They are assumed already allocated in
caller.  

This function is not used in this particular version of flca (filtering
is done within cross_cor), but is
included as it may be useful in the future.

*/

{

  i4 i,j;
  double *ina, *inb;
  double *filter, *kx, *ky;
  double normfac;
  fftw_complex *outa, *outb;
  fftw_plan pa, pb, pbacka, pbackb;
  outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
     nx * ((ny / 2) + 2));
  outb = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * 
     nx * ((ny / 2) + 2));	/* valgrind sometimes complains */
  ina = (double *) fftw_malloc (sizeof (double) * nx * ny);
  inb = (double *) fftw_malloc (sizeof (double) * nx * ny);
  filter = (double *) fftw_malloc (sizeof (double) * nx * ny);
  kx=(double *) fftw_malloc (sizeof(double)*nx);
  ky=(double *) fftw_malloc (sizeof(double)*ny);
  make_freq(kx,nx);
  make_freq(ky,ny);
  gaussfilt(filter,kx,ky,nx,ny,kr);
      for (i = 0; i < nx * ny; i++)
	{
	  *(ina + i) = (double) 0.;
	  *(inb + i) = (double) 0.;
	}
      for (i = 0; i < nx * ((ny / 2 + 1)); i++)
      {
#ifdef COMPLEXH
        /* If complex.h included, do this: */
        *(outa+i)=0.+I*0.;
        *(outb+i)=0.+I*0.;
#else
        /* if complex.h not included, do this: */
	outa[i][0] = 0.;
	outa[i][1] = 0.;
	outb[i][0] = 0.;
	outb[i][1] = 0.;
#endif
      }
      /* set up plans for FFTs: */
      pa = fftw_plan_dft_r2c_2d (nx, ny, ina, outa, FFTW_MEASURE);
      pb = fftw_plan_dft_r2c_2d (nx, ny, inb, outb, FFTW_MEASURE);
      pbacka = fftw_plan_dft_c2r_2d (nx, ny, outa, ina, FFTW_MEASURE);
      pbackb = fftw_plan_dft_c2r_2d (nx, ny, outb, inb, FFTW_MEASURE);

    for (i = 0; i < nx * ny; i++)

    {
/*		printf("1st loop: i = %d, *(arr+i)= %g, *(barr+i) = %g\n",
				i,*(arr+i),*(barr+i)); */

      /* copy from input doubles to fftw variables */

      *(ina + i) = (double) (*(arr + i));
      *(inb + i) = (double) (*(barr + i));
    }

  /* actually do the forward FFTs: */

  fftw_execute (pa);
  fftw_execute (pb);
 /* calculate normalization factor */

  normfac = (1. / ((double) nx * ny));

/* Now apply the gaussian filter to the FFT'd data in frequency space */

    for (i=0;i<nx;i++)
    {
      for (j=0;j<(ny/2)+1;j++)
      {
#ifdef COMPLEXH
        /* if complex.h included: */
        outa[i*((ny/2)+1)+j]=outa[i*((ny/2)+1)+j]*filter[i*ny+j]*normfac;
        outb[i*((ny/2)+1)+j]=outb[i*((ny/2)+1)+j]*filter[i*ny+j]*normfac;
#else
        /* If complex.h not included: */
        outa[i*((ny/2)+1)+j][0]=outa[i*((ny/2)+1)+j][0]*filter[i*ny+j]*normfac;
        outa[i*((ny/2)+1)+j][1]=outa[i*((ny/2)+1)+j][1]*filter[i*ny+j]*normfac;
        outb[i*((ny/2)+1)+j][0]=outb[i*((ny/2)+1)+j][0]*filter[i*ny+j]*normfac;
        outb[i*((ny/2)+1)+j][1]=outb[i*((ny/2)+1)+j][1]*filter[i*ny+j]*normfac;
#endif

      }
    }

/* now do the inverse transform to get filtered images: */

    fftw_execute (pbacka);
    fftw_execute (pbackb);

  for (i = 0; i < nx * ny; i++)
    {
      *(outarr+i)=(double) (*(ina+i));
      *(outbarr+i)=(double) (*(inb+i));
    }

/* Free the plans and locally created arrays */

      fftw_free (outa);
      fftw_free (outb);
      fftw_free (filter);
      fftw_free (kx);
      fftw_free (ky);
      fftw_free (ina);
      fftw_free (inb);
      fftw_destroy_plan (pa);
      fftw_destroy_plan (pbacka);
      fftw_destroy_plan (pb);
      fftw_destroy_plan (pbackb);

return 0;

}

i4 writeimage (char *fname, double *arr, i4 nx, i4 ny)
{

/* Function to write array dimensions, and then write out array */

  FILE *f1;
  i4 i, ier, ibe, ise, vcid, vcidtmp;
  f4 *farr;
  i4 nxtmp, nytmp;
  vcid = 2136967593;
  vcidtmp = vcid;
  nxtmp = nx;
  nytmp = ny;
  ibe = is_large_endian ();
  ise = 0;
  if (ibe == 0)
    ise = 1;	/* set flag for byteswapping if small endian */

  /* get the file fname open for binary write */

  f1 = fopen (fname, "wb");
  ier = 0;
  if (f1 == NULL)
    {
      printf ("writeimage: cannot open file %s\n", fname);
      exit (1);
    }

  if (ise)
    {
      byteswapflca ((void *) &vcidtmp, 1, sizeof (i4));
      byteswapflca ((void *) &nxtmp, 1, sizeof (i4));
      byteswapflca ((void *) &nytmp, 1, sizeof (i4));
    }

  fwrite (&vcidtmp, sizeof (i4), 1, f1);  /* write vel_ccor id integer 1st */

  fwrite (&nxtmp, sizeof (i4), 1, f1);
  fwrite (&nytmp, sizeof (i4), 1, f1);


/*
      printf("\n\nnx,ny wrote out to file arr = %d,%d\n",nx,ny);
*/

  /* create temporary f4 array to write out */

  farr = (f4 *) malloc (sizeof (f4) * nx * ny);

  /* now fill the array */

  for (i = 0; i < nx * ny; i++)
    {
      *(farr + i) = (f4) * (arr + i);
    }

  if (ise)
    byteswapflca ((void *) farr, nx * ny, sizeof (f4));
  /* byteswap if small endian */

  /* write out the array */
  fwrite (farr, sizeof (f4), nx * ny, f1);

  /* free temp. array and close file */
  free (farr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 write2images (char *fname, double *arr, double *barr, 
         i4 nx, i4 ny)
{

/* Function to write array dimensions, and write out 2 arrays arr and barr
 * while converting them to single precision f4 arrays  */

  FILE *f1;
  i4 ier, i, ise, ibe;
  i4 nxtmp, nytmp, vcid, vcidtmp;
  f4 *farr, *fbarr;
  nxtmp = nx;
  nytmp = ny;
  vcid = 2136967593;
  vcidtmp = vcid;
  ibe = is_large_endian ();
  ise = 0;
  if (ibe == 0) ise = 1;
  if (ise)		/* if small endian, byteswap nxtmp and nytmp */
    {
      byteswapflca ((void *) &vcidtmp, 1, sizeof (i4));
      byteswapflca ((void *) &nxtmp, 1, sizeof (i4));
      byteswapflca ((void *) &nytmp, 1, sizeof (i4));
    }

  /* open the file fname for a binary write */

  f1 = fopen (fname, "wb");
  ier = 0;
  if (f1 == NULL)
    {
      printf ("write2images: cannot open file %s\n", fname);
      exit (1);
    }

  fwrite (&vcidtmp, sizeof (i4), 1, f1);	/* write vel_ccor id flag */

  fwrite (&nxtmp, sizeof (i4), 1, f1);
  fwrite (&nytmp, sizeof (i4), 1, f1);

/*
      printf("\n\nnx,ny wrote out to file arr = %d,%d\n",nx,ny);
*/

  /* create space for the temporary f4 arrays farr and fbarr */

  farr = (f4 *) malloc (sizeof (f4) * nx * ny);
  fbarr = (f4 *) malloc (sizeof (f4) * nx * ny);

  /* fill the temporary arrays */

  for (i = 0; i < nx * ny; i++)
    {
      *(farr + i) = (f4) * (arr + i);
      *(fbarr + i) = (f4) * (barr + i);
    }

  /* now write out the 2 arrays */

  if (ise)	/* byteswap if small endian */
    {
      byteswapflca ((void *) farr, nx * ny, sizeof (f4));
      byteswapflca ((void *) fbarr, nx * ny, sizeof (f4));
    }

  fwrite (farr, sizeof (f4), nx * ny, f1);
  fwrite (fbarr, sizeof (f4), nx * ny, f1);

  /* free temp arrays and close file */

  free (farr);
  free (fbarr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 write3images (char *fname, double *arr, double *barr, double *carr,
	      i4 nx, i4 ny)
{

/* Function to write array dimensions, and write out 3 arrays arr,barr, and carr
 * while converting them to single precision f4 arrays  */

  FILE *f1;
  i4 ier, i, ibe, ise;
  i4 nxtmp, nytmp, vcid, vcidtmp;
  f4 *farr, *fbarr, *fcarr;
  nxtmp = nx;
  nytmp = ny;
  vcid = 2136967593;
  vcidtmp = vcid;
  ibe = is_large_endian ();
  ise = 0;
  if (ibe == 0) ise = 1;	/* test for small endian for doing byteswaps */
  if (ise)			/* byteswap nxtmp, nytmp if small endian */
    {
      byteswapflca ((void *) &vcidtmp, 1, sizeof (i4));
      byteswapflca ((void *) &nxtmp, 1, sizeof (i4));
      byteswapflca ((void *) &nytmp, 1, sizeof (i4));
    }

  /* open the file fname for a binary write */

  f1 = fopen (fname, "wb");
  ier = 0;
  if (f1 == NULL)
    {
      printf ("write3images: cannot open file %s\n", fname);
      exit (1);
    }

  fwrite (&vcidtmp, sizeof (i4), 1, f1);

  fwrite (&nxtmp, sizeof (i4), 1, f1);
  fwrite (&nytmp, sizeof (i4), 1, f1);

/*
      printf("\n\nnx,ny wrote out to file arr = %d,%d\n",nx,ny);
*/

  /* create space for the temporary f4 arrays farr, fbarr, and fcarr */

  farr = (f4 *) malloc (sizeof (f4) * nx * ny);
  fbarr = (f4 *) malloc (sizeof (f4) * nx * ny);
  fcarr = (f4 *) malloc (sizeof (f4) * nx * ny);

  /* fill the temporary arrays */

  for (i = 0; i < nx * ny; i++)
    {
      *(farr + i) = (f4) * (arr + i);
      *(fbarr + i) = (f4) * (barr + i);
      *(fcarr + i) = (f4) * (carr + i);
    }

  if (ise)			/* if small endian, byteswap the arrays */
    {
      byteswapflca ((void *) farr, nx * ny, sizeof (f4));
      byteswapflca ((void *) fbarr, nx * ny, sizeof (f4));
      byteswapflca ((void *) fcarr, nx * ny, sizeof (f4));
    }

  /* now write out the 3 arrays */

  fwrite (farr, sizeof (f4), nx * ny, f1);
  fwrite (fbarr, sizeof (f4), nx * ny, f1);
  fwrite (fcarr, sizeof (f4), nx * ny, f1);

  /* free temp arrays and close file */

  free (farr);
  free (fbarr);
  free (fcarr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 shift2d (double *arr, i4 nx, i4 ny, i4 ishift, i4 jshift)
{

/* Circular shift of the x,y indices of array *arr by ishift,jshift */
/* This function is similar to the shift function in IDL.  nx, ny
 * are the assumed dimensions of the array */

  double *temp;
  i4 i, j, ii, jj;
  temp = (double *) malloc (sizeof (double) * nx * ny);
  for (i = 0; i < nx; i++)
    {
      ii = (i + ishift) % nx;	/* ii = (i + ishift) modulo nx */

      for (j = 0; j < ny; j++)
	{
	  jj = (j + jshift) % ny;	/* jj = (j+jshift) modulo ny */

	  /* Now members of temp array get shifted: */

	  *(temp + ii * ny + jj) = *(arr + i * ny + j);
	}
    }

  /* Now copy temp array back into arr, then destroy temp and return */

  memcpy ((void *) arr, (void *) temp, nx * ny * sizeof (double));
  free (temp);
  return 0;
}

i4 maxloc (double *arr, i4 xsize)
{

/* finds the location of the maximum of the double array *arr and returns it. */

  i4 i, location;
  double amax;
  /* initialize amax and location to 0th element */
  amax = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) > amax)
	{
	  amax = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 imaxloc (i4 * arr, i4 xsize)
{

/* finds the location of the maximum of the i4 array *arr and returns it. */

  i4 i, location;
  i4 amax;
  /* initialize amax and location to 0th element */
  amax = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) > amax)
	{
	  amax = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 minloc (double *arr, i4 xsize)
{

/* finds the location of the minimum of the double array *arr and returns it. */

  i4 i, location;
  double amin;
  /* initialize amin and location to 0th element */
  amin = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) < amin)
	{
	  amin = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 iminloc (i4 * arr, i4 xsize)
{

/* finds the location of the minimum of the i4 array *arr and returns it. */

  i4 i, location;
  i4 amin;
  /* initialize amin and location to 0th element */
  amin = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) < amin)
	{
	  amin = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 interpcc2d (double *fdata, double xmiss, i4 nx, i4 ny, 
    double *xwant, i4 nxinterp, double *ywant, i4 nyinterp, double **finterp)
{
  /*
   * This function does cubic convolution interpolation onto an array 
   * finterp from data defined in array fdata.  nx, ny are the
   * assumed dimensions of fdata, and nxinterp, nyinterp are the
   * assumed dimensions of finterp.  The values of x,y at which
   * the interpolation is desired are passed in through the arrays
   * xwant and ywant, which are dimensioned nxinterp and nyinterp,
   * respectively.  It is assumed that xwant, ywant are in units of
   * the indices of the original data array (fdata), 
   * treated as floating point (double precision, actually) 
   * numbers. Arrays fdata, xwant, and ywant are passed in
   * as pointers; The array finterp is defined in this function
   * as a "double" pointer and the array is created and passed back to
   * the calling function.  In the calling function, finterp is declared
   * as a pointer, but when it is passed into this function as
   * an argument, the address of the pointer is used.
   * 
   * if any of the datapoints within a kernel weighting distance of
   * xwant and ywant are equal to xmiss,
   * the returned value of finterp is also set to xmiss.  xmiss is a user-
   * defineable calling argument.
   */

  double *cdata;
/*  double txt, tyt, xint, yint, ftmp, xmiss = 0.; */
  double txt, tyt, xint, yint, ftmp;

  /* Logic for a user-defined value of xmiss has been added.  Previously
   * was just set to 0 as a local variable */

  double tx, ty, rx, ry;
  i4 i, ii, j, jj, itemp, jtemp, izero, jzero, databad;
/*  i4 transp; */

  /* Now, create the cdata array, bigger by 1 gp than fdata
   * all around the borders: */

  cdata = (double *) malloc (sizeof (double) * (nx + 2) * (ny + 2));

  /* Now fill the interior of cdata with fdata */

  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
	{
	  *(cdata + (i + 1)*(ny + 2) + (j + 1)) = *(fdata + i*ny + j);
	}
    }

  /*
   * The basic concept for filling in edges and corners of cdata is this:
   * The edge point is equal to 3*(value of adjacent point)
   * -3*value(next to adjacent point) + 1*value(3rd point).  This
   * prescription yields an extrapolation which is consistent with
   * a 3rd (or is it 4th?) order Taylor expansion of the function
   * evaluated at the last real gridpoint, and extrapolated to the
   * edge point.  This procedure is followed
   * thoughout here, though I think it isn't really correct for the
   * corner points because there I think an expansion from both
   * both directions should be done.  But no harm seems to be done
   * to the results.
   */

  /* Fill in the edges of cdata: */

  for (j = 0; j < ny; j++)
    {

      /* left and right edges: */

      *(cdata + 0*(ny + 2) + (j+1)) = *(fdata + 2*ny + j)
	- 3. * (*(fdata + 1*ny + j)) + 3. * (*(fdata + 0*ny + j));

      *(cdata + (nx + 1)*(ny + 2) + (j + 1)) = *(fdata + (nx - 3)*ny + j)
	- 3. * (*(fdata + (nx - 2)*ny + j)) + 3. * (*(fdata + (nx - 1)*ny + j));
    }
  for (i = 0; i < nx; i++)
    {

      /* bottom and top edges: */

      *(cdata + (i + 1)*(ny + 2) + 0) = *(fdata + i*ny + 2)
	- 3. * (*(fdata + i*ny + 1)) + 3. * (*(fdata + i*ny + 0));

      *(cdata + (i + 1)*(ny + 2) + ny + 1) = *(fdata + i*ny + ny - 3)
	- 3. * (*(fdata + i*ny + ny - 2)) + 3. * (*(fdata + i*ny + ny - 1));
    }

  /* Now fill in the 4 corners: */

  *(cdata + 0*(nx + 2) + 0) = 
    3. * (*(cdata + 1*(ny + 2) + 0)) -
    3. * (*(cdata + 2*(ny + 2) + 0)) + *(cdata + 3*(ny + 2) + 0);

  *(cdata + (nx + 1)*(ny + 2) + 0) = 
    3. * (*(cdata + nx*(ny + 2) + 0)) -
    3. * (*(cdata + (nx - 1)*(ny + 2) + 0)) + *(cdata +
						    (nx - 2)*(ny + 2) + 0);

  *(cdata + 0*(ny + 2) + ny + 1) = 
    3. * (*(cdata + 0*(ny + 2) + ny)) -
    3. * (*(cdata + 0*(ny + 2) + ny - 1)) + *(cdata + 0*(ny + 2) + ny - 2);

  *(cdata + (nx + 1)*(ny + 2) + ny + 1) =
    3. * (*(cdata + nx*(ny + 2) + ny + 1)) -
    3. * (*(cdata + (nx - 1)*(ny + 2) + ny + 1)) + *(cdata +
						       (nx - 2)*(ny + 2) +
						       ny + 1);

  /* Now create the space for finterp */

  *finterp = (double *) malloc (sizeof (double) * nxinterp * nyinterp);

  /* Now interpolate onto the desired grid */

  for (i = 0; i < nxinterp; i++)
    {
      /* starting the outer loop over x */

      xint = *(xwant + i);

      /* make sure izero is in bounds */

      itemp = ((i4) xint > 0) ? (i4) xint : 0;
      izero = (itemp < (nx - 2)) ? itemp : nx - 2;
      for (j = 0; j < nyinterp; j++)
	{
	  /* starting the outer loop over y */

	  yint = *(ywant + j);
	  if ((yint < 0.) || (yint > (double) (ny - 1))
	      || ((xint < 0) || (xint > (double) (nx - 1))))
	    {
	      /* if data missing, set interp to xmiss */

/* Debug
              printf("interpccd2: i=%d,j=%d gets finterp[i,j] set to xmiss\n",
                    i,j);
*/
	      *(*finterp + i * nyinterp + j) = xmiss;
	    }
	  else
	    {
	      /* make sure jzero is in bounds */

	      jtemp = ((i4) yint > 0) ? (i4) yint : 0;
	      jzero = (jtemp < (ny - 2)) ? jtemp : ny - 2;

	      /* initialize the temporary finterp value */

	      ftmp = (double) 0.;

	      /* start the innermost loops over neighboring
	       * data points*/

              databad=0;
	      for (ii = -1; ii < 3; ii++)
		{
		  txt = xint - (double) (izero + ii);
		  tx = (double) fabs (txt);

		  /* evaluate kernel wt function r(tx): */

		  /* Note no testing for out of range 
		   * values of |tx| or |ty| > 2 --
		   * we assume the tx, ty are properly
		   * computed such that their absolute
		   * value never exceeds 2. */

		  rx = (tx >= (double) 1.0) ?
		    (((((double) (-0.5)) * tx +
		       ((double) 2.5)) * tx) -
		     (double) 4.) * tx + (double) 2. :
		    (((double) (1.5)) * tx -
		     ((double) (2.5))) * tx * tx + (double) 1.;

		  for (jj = -1; jj < 3; jj++)
		    {

		      tyt = yint - (double) (jzero + jj);
		      ty = (double) fabs (tyt);

		      /* evaluate kernel weighting
		       * function r(ty): */

		      ry = (ty >= (double) 1.0) ?
			(((((double) (-0.5)) * ty +
			   ((double) 2.5)) * ty) -
			 (double) 4.) * ty + (double) 2. :
			(((double) (1.5)) * ty -
			 ((double) (2.5))) * ty * ty + (double) 1.;

		      /* do the cubic convolution
		       * over the neighboring data
		       * points, using the x and
		       * y evaluated kernel weighting
		       * functions rx and ry: */

		      ftmp = ftmp +
			*(cdata + (izero + 1 + ii)*(ny + 2)
			  + jzero + 1 + jj) * rx*ry;
                      if( *(cdata+(izero+1+ii)*(ny+2)+jzero+1+jj) == xmiss)
                          databad=1;
		    }
		}
	      /* now assign this value to interpolated
	         array, unless one of the values was xmiss: */
              if(databad)
              {
/* Debug
                 printf("interpcc2d: i=%d,j=%d gives databad\n",i,j);
*/
                 *(*finterp + i*nyinterp + j) = xmiss;
              }
              else
              {
	         *(*finterp + i*nyinterp + j) = ftmp;
              }
	    }
	}
    }


/* DEBUG
  transp=1;
  writeimage("cdata.dat",cdata,nx+2,ny+2,transp);
*/

  /* free the cdata space */
  free (cdata);

  /* we're done */

  return 0;
}

i4 byteswapflca (unsigned char *arr, i4 arrsize, i4 nbpw)
/* Pretty simple:  arr is input array, which is byte-swapped in place,
   nbpw is the number of bytes per word, and arrsize is the size of the array
   (in units of nbpw bytes).  It is assumed that arr has
   already have been correctly defined and allocated in the calling program. */
{
  i4 i, j;
  unsigned char temp;
  for (i = 0; i < arrsize; i++)	/* the loop over the array elements */
    {
      for (j = 0; j < nbpw/2; j++)/* the loop over bytes in a single element */
	{
	  temp = *(arr + i*nbpw + (nbpw - j - 1));
	  *(arr + i*nbpw + (nbpw - j - 1)) = *(arr + i*nbpw + j);
	  *(arr + i*nbpw + j) = temp;
	}
    }
  return 0;
}

i4 is_large_endian ()
/* This function returns 1 if it is a large endian machine, 0 otherwise */
{
  const unsigned char fakeword[4] = { 0xFF, 0x00, 0xFF, 0x00 };
  i4 realword = *((i4 *) fakeword);
  if (realword == 0xFF00FF00)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}
