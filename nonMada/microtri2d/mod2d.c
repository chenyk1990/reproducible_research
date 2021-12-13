/* 2D acoustic wavefield modeling using the pseudo-spectral method  */
/*
  Copyright (C) 2020 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  
  This is the simplified 2D version of the 3D modeling code
  https://github.com/chenyk1990/passive_imaging/mod3d.c
*/

#include <rsf.h>

float Ricker(float t, float f0, float t0, float A) 
/*< ricker wavelet:
 * f0: peak frequency
 * t0: time lag
 * A: amplitude
 * ************************>*/
{
        float x=pow(SF_PI*f0*(t-t0),2);
        return -A*exp(-x)*(1-2*x);
}


static kiss_fft_cpx *shape;

void ricker_init(int nfft   /* time samples */, 
		 float freq /* frequency */,
		 int order  /* derivative order */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w;
    kiss_fft_cpx cw;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 1./(nfft*freq);
 
    shape = (kiss_fft_cpx*) sf_complexalloc(nw);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w *= w;

	switch (order) {
	    case 2: /* half-order derivative */
		cw.r = 2*SF_PI/nfft;
		cw.i = iw*2*SF_PI/nfft;
		cw = sf_csqrtf(cw);
		shape[iw].r = cw.r*w*expf(1-w)/nfft;
		shape[iw].i = cw.i*w*expf(1-w)/nfft;
		break;
	    case 0:
	    default:
		shape[iw].r = w*expf(1-w)/nfft;
		shape[iw].i = 0.;
		break;
	}
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_cset(shape);
}

void ricker_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

static int nx, nz, nx2, nz2, nbt, nbb, nbl, nbr;
static float ct, cb, cl, cr;
static float *wt, *wb, *wl, *wr;

void abc_cal(int abc /* decaying type*/,
             int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< find absorbing coefficients >*/
{
    int ib;
    /*const float pi=SF_PI;*/
    if(!nb) return;
    switch(abc) {
    default:
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ib)
#endif
        for(ib=0; ib<nb; ib++){
	    w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
	}
    }
}

void abc_init(int n1,  int n2    /*model size*/,
	      int n12, int n22   /*padded model size*/,
	      int nb1, int nb2   /*top, bottom*/,
	      int nb3, int nb4   /*left, right*/,
	      float c1, float c2 /*top, bottom*/,
	      float c3, float c4 /*left, right*/)
/*< initialization >*/
{
    int c;
    nz = n1;
    nx = n2;
    nz2= n12;
    nx2= n22;
    nbt = nb1;
    nbb = nb2;
    nbl = nb3;
    nbr = nb4;
    ct = c1;
    cb = c2;
    cl = c3;
    cr = c4;
    if(nbt) wt =  sf_floatalloc(nbt);
    if(nbb) wb =  sf_floatalloc(nbb);
    if(nbl) wl =  sf_floatalloc(nbl);
    if(nbr) wr =  sf_floatalloc(nbr);
    c=0;
    abc_cal(c,nbt,ct,wt);
    abc_cal(c,nbb,cb,wb);
    abc_cal(c,nbl,cl,wl);
    abc_cal(c,nbr,cr,wr);
}
   

void abc_close(void)
/*< free memory allocation>*/
{
    if(nbt) free(wt);
    if(nbb) free(wb);
    if(nbl) free(wl);
    if(nbr) free(wr);
}

void abc_apply(float *a /*2-D matrix*/) 
/*< boundary decay>*/
{
    int i;
    int iz, ix;

    /* top */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix,i)
{
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=nbl; ix < nx-nbr; ix++) {
	  i = nz2*ix + iz;
	  a[i] *= wt[iz];
        }
    }
    /* bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=nbl; ix < nx-nbr; ix++) {
	  i = nz2*ix + nz-1-iz;
	  a[i] *= wb[iz];
        }
    }
    /* left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=nbt; iz < nz-nbb; iz++) {  
        for (ix=0; ix < nbl; ix++) {
	  i = nz2*ix + iz;
	  a[i] *= wl[ix];
        }
    }
    /* right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=nbt; iz < nz-nbb; iz++) {  
        for (ix=0; ix < nbr; ix++) {
	  i = nz2*(nx-1-ix) + iz;
          a[i] *= wr[ix];
        }
    }
    /* top left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nbl; ix++) {
	  i = nz2*ix + iz;
	  a[i] *= iz>ix? wl[ix]:wt[iz];
        }
    }
    /* top right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nbr; ix++) {
	  i = nz2*(nx-1-ix) + iz;
	  a[i] *= iz>ix? wr[ix]:wt[iz];
        }
    }
    /* bottom left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nbl; ix++) {
	  i = nz2*ix + nz-1-iz;
          a[i] *= iz>ix? wl[ix]:wb[iz];
        }
    }
    /* bottom right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nbr; ix++) {
	  i = nz2*(nx-1-ix) + nz-1-iz;
          a[i] *= iz>ix? wr[ix]:wb[iz];
        }
    }

#ifdef _OPENMP
}
#endif
}






static bool cmplx;
static int n1, n2, nk;
static float wwt;

static float **ff=NULL;
static sf_complex **cc=NULL,*dd=NULL;

#ifdef SF_HAS_FFTW
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fftr_cfg cfg=NULL, icfg=NULL;
static kiss_fft_cfg cfg1=NULL, icfg1=NULL, cfg2=NULL, icfg2=NULL;
static kiss_fft_cpx **tmp=NULL, *ctrace2=NULL;
static sf_complex *trace2=NULL;
#endif

int fft2_init(bool cmplx1        /* if complex transform */,
	      int pad1           /* padding on the first axis */,
	      int nx,   int ny   /* input data size */, 
	      int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_init_threads();
    sf_warning("Using threaded FFTW3!\n");
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
#else
   int i2;
#endif
	
    cmplx = cmplx1;
	
    if (cmplx) {
	nk = n1 = kiss_fft_next_fast_size(nx*pad1);
		
#ifndef SF_HAS_FFTW
	cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
	icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);
#endif
    } else {
	nk = kiss_fft_next_fast_size(pad1*(nx+1)/2)+1;
	n1 = 2*(nk-1);
		
#ifndef SF_HAS_FFTW
	cfg  = kiss_fftr_alloc(n1,0,NULL,NULL);
	icfg = kiss_fftr_alloc(n1,1,NULL,NULL);
#endif
    }
		
    n2 = kiss_fft_next_fast_size(ny);

    if (cmplx) {
	cc = sf_complexalloc2(n1,n2);
    } else {
	ff = sf_floatalloc2(n1,n2);
    }
    dd = sf_complexalloc(nk*n2);
	
#ifndef SF_HAS_FFTW
    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);
 	
    tmp =    (kiss_fft_cpx **) sf_alloc(n2,sizeof(*tmp));
    tmp[0] = (kiss_fft_cpx *)  sf_alloc(nk*n2,sizeof(kiss_fft_cpx));
    for (i2=0; i2 < n2; i2++) {
	tmp[i2] = tmp[0]+i2*nk;
    }
	
    trace2 = sf_complexalloc(n2);
    ctrace2 = (kiss_fft_cpx *) trace2;
#endif

    *nx2 = n1;
    *ny2 = n2;
	
    wwt =  1.0/(n1*n2);
	
    return (nk*n2);
}

void fft2(float *inp      /* [n1*n2] */, 
	  sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==cfg) {
	cfg = cmplx? 
	    fftwf_plan_dft_2d(n2,n1,
			      (fftwf_complex *) cc[0],
			      (fftwf_complex *) dd,
			      FFTW_FORWARD, FFTW_MEASURE):
	    fftwf_plan_dft_r2c_2d(n2,n1,
				  ff[0], (fftwf_complex *) dd,
				  FFTW_MEASURE);
	if (NULL == cfg) sf_error("FFTW failure.");
    }
#endif

    /* FFT centering */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    if (cmplx) {
		cc[i2][i1] = sf_cmplx(((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1],0.);
	    } else {
		ff[i2][i1] = (i2%2)? -inp[i2*n1+i1]:inp[i2*n1+i1];
	    }
	}
    }
    
#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);

    for (i1=0; i1 < nk*n2; i1++)
      out[i1] = dd[i1];
#else	
    for (i2=0; i2 < n2; i2++) {
	if (cmplx) {
	    kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc[i2],tmp[i2],1);
	} else {
	    kiss_fftr (cfg,ff[i2],tmp[i2]);
	}
    }
	
    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(cfg2,tmp[0]+i1,ctrace2,nk);
	for (i2=0; i2<n2; i2++) {
	    out[i2*nk+i1] = trace2[i2];
	}
    }
#endif
}

void ifft2_allocate(sf_complex *inp /* [nk*n2] */)
/*< allocate inverse transform >*/
{
  /* for backward compatibility */
}

void ifft2(float *out      /* [n1*n2] */, 
	   sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==icfg) {
      icfg = cmplx? 
	fftwf_plan_dft_2d(n2,n1,
			  (fftwf_complex *) dd, 
			  (fftwf_complex *) cc[0],
			  FFTW_BACKWARD, FFTW_MEASURE):
	fftwf_plan_dft_c2r_2d(n2,n1,
			      (fftwf_complex *) dd, ff[0],
			      FFTW_MEASURE);
      if (NULL == icfg) sf_error("FFTW failure.");
    }
#endif

#ifdef SF_HAS_FFTW
    for (i1=0; i1 < nk*n2; i1++)
      dd[i1] = inp[i1];

    fftwf_execute(icfg);
#else
    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(icfg2,(kiss_fft_cpx *) (inp+i1),ctrace2,nk);
		
	for (i2=0; i2<n2; i2++) {
	    tmp[i2][i1] = ctrace2[i2];
	}
    }
    for (i2=0; i2 < n2; i2++) {
	if (cmplx) {
	    kiss_fft_stride(icfg1,tmp[i2],(kiss_fft_cpx *) cc[i2],1);
	} else {
	    kiss_fftri(icfg,tmp[i2],ff[i2]);
	}
    }
#endif
    
    /* FFT centering and normalization */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    if (cmplx) {
		out[i2*n1+i1] = (((i2%2==0)==(i1%2==0))? wwt:-wwt) * crealf(cc[i2][i1]);
	    } else {
		out[i2*n1+i1] = (i2%2? -wwt: wwt)*ff[i2][i1];
	    }
	}
    }
}

void fft2_finalize()
/*< clean up fftw >*/
{
/* make sure everything is back to its pristine state */
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_cleanup_threads();
#endif
    fftwf_destroy_plan(cfg);
    fftwf_destroy_plan(icfg);
    fftwf_cleanup();
    cfg=NULL;
    icfg=NULL;
#else
    if (NULL != cfg) { free(cfg); cfg=NULL; }
    if (NULL != icfg) { free(icfg); icfg=NULL; }
    if (NULL != cfg1) { free(cfg1); cfg1=NULL; }
    if (NULL != icfg1) { free(icfg1); icfg1=NULL; }
    if (NULL != cfg2) { free(cfg2); cfg2=NULL; }
    if (NULL != icfg2) { free(icfg2); icfg2=NULL; }
    if (NULL != tmp) { free(*tmp); free(tmp); tmp=NULL; }
    if (NULL != trace2) { free(trace2); trace2=NULL; }
#endif
    if (cmplx) {
      if (NULL != cc) { free(*cc); free(cc); cc=NULL; }
    } else {
      if (NULL != ff) { free(*ff); free(ff); ff=NULL; }
    }
    if (NULL != dd) { free(dd); dd=NULL; }
}



typedef struct Pspar {
  /*survey parameters*/
  int   nx, nz;
  float dx, dz;
  int   n_srcs;
  int   *spx, *spz;
  int   gpz, gpx, gpl;
  int   gpz_v, gpx_v, gpl_v;
  int   snap;
  /*fft related*/
  bool  cmplx;
  int   pad1;
  /*absorbing boundary*/
  bool abc;
  int nbt, nbb, nbl, nbr;
  float ct,cb,cl,cr;
  /*source parameters*/
  int src; /*source type*/
  int nt;
  float dt,*f0,*t0,*A;
  /*misc*/
  bool verb, ps;
  float vref;
} * pspar; /*psp parameters*/
/*^*/



int psp(float **wvfld, float **dat, float **dat_v, float *img, float *vel, pspar par, bool mig)
/*< pseudo-spectral wave extrapolation >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k,t;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*prev,*lapl;

    /*source*/
    float **rick;
    float freq;
    int fft_size;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
    src   = par->src;
    nt    = par->nt;
    dt    = par->dt;
    f0    = par->f0;
    t0    = par->t0;
    A     = par->A;
    verb  = par->verb;
    ps    = par->ps;
    vref  = par->vref;
    

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
#else
    nth = 1;
#endif
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    if (!mig && src==0) {
      rick = sf_floatalloc2(nt,n_srcs);
      for (i=0; i<n_srcs; i++) {
	for (it=0; it<nt; it++) {
	  rick[i][it] = 0.f;
	}
	rick[i][(int)(t0[i]/dt)] = A[i]; /*time delay*/
	freq = f0[i]*dt;           /*peak frequency*/
	fft_size = 2*kiss_fft_next_fast_size((nt+1)/2);
	ricker_init(fft_size, freq, 0);
	sf_freqfilt(nt,rick[i]);
	ricker_close();
      }
    } else rick = NULL;

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    if (mig) { /* time-reversal propagation */
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* modeling */
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);

	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

	if (mig) {
	  /* inject data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      curr[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat[ix][it];
	    }
	  }
	  if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
	    for (iz = 0; iz < gpl_v; iz++) {
	      curr[gpz_v+iz+(gpx_v)*nz2] += vv[gpz_v+iz+(gpx_v)*nz]*dat_v[iz][it];
	    }
	  }
	} else {
	  t = it*dt;
	  for (i=0; i<n_srcs; i++) {
	    for(ix=-1;ix<=1;ix++) {
	      for(iz=-1;iz<=1;iz++) {
		ik = spz[i]+iz+nz*(spx[i]+ix);
		j = spz[i]+iz+nz2*(spx[i]+ix);
		if (src==0) {
		  curr[j] += vv[ik]*rick[i][it]/(abs(ix)+abs(iz)+1);
		} else {
		  curr[j] += vv[ik]*Ricker(t, f0[i], t0[i], A[i])/(abs(ix)+abs(iz)+1);
		}
	      }
	    }
	  }
	}
	
	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}

	if (!mig) {
	  /* record data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      dat[ix][it] = curr[gpz+(ix+gpx)*nz2];
	    }
	  }
	  if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
	    for (iz = 0; iz < gpl_v; iz++) {
	      dat_v[iz][it] = curr[gpz_v+iz+(gpx_v)*nz2];
	    }
	  }
	}

	/* save wavefield */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    if (mig) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz)
#endif
      for (ix = 0; ix < nx1; ix++) {
	for (iz = 0; iz < nz1; iz++) {
	  img[iz + nz1*ix] = curr[iz+nbt + (ix+nbl)*nz2];
	}
      }
    }

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(prev);
    free(cwave);
    free(cwavem);
    
    return 0;
}

int psp2(float **wvfld1, float **wvfld, float **dat, float **dat_v, float *img, float *vel, pspar par, bool adj)
/*< pseudo-spectral de-migration/migration >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*prev,*lapl;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
    src   = par->src;
    nt    = par->nt;
    dt    = par->dt;
    f0    = par->f0;
    t0    = par->t0;
    A     = par->A;
    verb  = par->verb;
    ps    = par->ps;
    vref  = par->vref;
    

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
#else
    nth = 1;
#endif
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }
    if (adj) {
      for (iz=0; iz < nz1*nx1; iz++) {
        img[iz] = 0.;
      }
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    if (adj) { /* RTM */
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* de-migration */
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);

	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

	if (adj) {
	  /* inject data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      curr[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat[ix][it];
	    }
	  }
	  if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
	    for (iz = 0; iz < gpl_v; iz++) {
	      curr[gpz_v+iz+(gpx_v)*nz2] += vv[gpz_v+iz+(gpx_v)*nz]*dat_v[iz][it];
	    }
	  }
	} else {
	  /*adj of cross-correlation*/
	  if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,ik)
#endif
	    for (ix = 0; ix < nx1; ix++) {
	      for (iz = 0; iz < nz1; iz++) {
		i = iz + nz1*ix;
		j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
		ik = iz+nbt + (ix+nbl)*nz; /* padded grid */
		curr[j] += wvfld1[it/snap][i]*img[i]; // !!!IMPORTANT!!! vv[ik]*
	      }
	    }
	  }
	}

	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}
	
	if (adj) {
	  /* cross-correlation imaging condition */
	  if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,ik)
#endif
	    for (ix = 0; ix < nx1; ix++) {
	      for (iz = 0; iz < nz1; iz++) {
		i = iz + nz1*ix;
		j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
		ik = iz+nbt + (ix+nbl)*nz; /* padded grid */
		img[i] += wvfld1[it/snap][i]*curr[j]; // !!!IMPORTANT!!! vv[ik]*
	      }
	    }
	  }
	} else {
	  /* record data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      dat[ix][it] = curr[gpz+(ix+gpx)*nz2];
	    }
	  }
	  if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
	    for (iz = 0; iz < gpl_v; iz++) {
	      dat_v[iz][it] = curr[gpz_v+iz+(gpx_v)*nz2];
	    }
	  }
	}

	/* save wavefield */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(prev);
    free(cwave);
    free(cwavem);
    
    return 0;
}

int psp3(float **wvfld, float **wvfld1, float **dat, float **dat1, float *img, float *vel, pspar par)
/*< pseudo-spectral back propagation of two receiver wavefields >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*curr1,*prev,*prev1,*lapl;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
    src   = par->src;
    nt    = par->nt;
    dt    = par->dt;
    f0    = par->f0;
    t0    = par->t0;
    A     = par->A;
    verb  = par->verb;
    ps    = par->ps;
    vref  = par->vref;
    

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
#else
    nth = 1;
#endif
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    curr1  = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    prev1  = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	curr1[iz] = 0.;
	prev[iz] = 0.;
	prev1[iz] = 0.;
    }
    for (iz=0; iz < nz1*nx1; iz++) {
        img[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    /* step backward in time */
    it1 = nt-1;
    it2 = -1;
    its = -1;	

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);
	
        /* first receiver wavefield */
	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

	/* inject data */
	if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	        curr[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat[ix][it];
	    }
	}

	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}
	
	/* second receiver wavefield */
	/* matrix multiplication */
	fft2(curr1,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr1[j];
		c += c - prev1[j];
		prev1[j] = old;
		c += wave[j]*vv[i];
		curr1[j] = c;
	    }
	}

	/* inject data */
	if (NULL!=dat1) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	        curr1[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat1[ix][it];
	    }
	}

	/*apply abc*/
	if (abc) {
	  abc_apply(curr1);
	  abc_apply(prev1);
	}

	/* cross-correlation imaging condition */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	    for (ix = 0; ix < nx1; ix++) {
	        for (iz = 0; iz < nz1; iz++) {
		    i = iz + nz1*ix;
		    j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
		    img[i] += curr[j]*curr1[j]; //!!!IMPORTANT!!!
		}
	    }
	}
	
	/* save wavefield */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
              wvfld1[it/snap][i] = curr1[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(curr1);
    free(prev);
    free(prev1);
    free(cwave);
    free(cwavem);
    
    return 0;
}

int dt2v2(float **wvfld, float *vel, pspar par)
/*< inplace adjusting the background wavefield w to d^2w/dt^2/v^2 >*/
{
  int nx1, nz1, wfnt;
  int ix, iz, it, i, j;
  float wfdt;
  float *vel2;
  
  nz1 = par->nz - par->nbt - par->nbb;
  nx1 = par->nx - par->nbl - par->nbr;
  wfnt = par->nt/par->snap;
  wfdt = par->snap*par->dt;

  vel2 = sf_floatalloc(nz1*nx1);

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i,j)
#endif
  for (ix=0; ix<nx1; ix++) {
    for (iz=0; iz<nz1; iz++) {
      i = ix*nz1+iz;
      j = (ix+par->nbl)*par->nz + iz + par->nbt;
      vel2[i] = vel[j]*vel[j]*wfdt*wfdt;
    }
  }
  
  /*2nd order forward FD*/
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i)
#endif
  for (ix=0; ix<nx1; ix++) {
    for (iz=0; iz<nz1; iz++) {
      i = iz + nz1*ix;
      wvfld[0][i] = (wvfld[2][i] - 2.f*wvfld[1][i] + wvfld[0][i])/vel2[i];
    }
  }

  for (it=1; it<wfnt-1; it++) { /*2nd order central FD*/
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i)
#endif
    for (ix=0; ix<nx1; ix++) {
      for (iz=0; iz<nz1; iz++) {
	i = iz + nz1*ix;
	wvfld[it][i] = (wvfld[it+1][i] - 2.f*wvfld[it][i] + wvfld[it-1][i])/vel2[i];
      }
    }
  } /*it iteration*/

  /*2nd order backward FD*/
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i)
#endif
  for (ix=0; ix<nx1; ix++) {
    for (iz=0; iz<nz1; iz++) {
      i = iz + nz1*ix;
      wvfld[wfnt-1][i] = (wvfld[wfnt-1][i] - 2.f*wvfld[wfnt-2][i] + wvfld[wfnt-3][i])/vel2[i];
    }
  }

  return 0;
}

int psp4(float **wvfld, float **wvfld0, float **dat, float **dat_v, float *vel, pspar par)
/*< pseudo-spectral wave extrapolation using wavefield injection >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k,t;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*prev,*lapl;

    /*source*/
    float **rick;
    float freq;
    int fft_size;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
    src   = par->src;
    nt    = par->nt;
    dt    = par->dt;
    f0    = par->f0;
    t0    = par->t0;
    A     = par->A;
    verb  = par->verb;
    ps    = par->ps;
    vref  = par->vref;
    

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
#else
    nth = 1;
#endif
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    if (src==0) {
      rick = sf_floatalloc2(nt,n_srcs);
      for (i=0; i<n_srcs; i++) {
	for (it=0; it<nt; it++) {
	  rick[i][it] = 0.f;
	}
	rick[i][(int)(t0[i]/dt)] = A[i]; /*time delay*/
	freq = f0[i]*dt;           /*peak frequency*/
	fft_size = 2*kiss_fft_next_fast_size((nt+1)/2);
	ricker_init(fft_size, freq, 0);
	sf_freqfilt(nt,rick[i]);
	ricker_close();
      }
    } else rick = NULL;

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    /* modeling */
    /* step forward in time */
    it1 = 0;
    it2 = nt;
    its = +1;

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);
	
	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

        if (NULL!=wvfld0) { /* wavefield injection */
          if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,ik)
#endif
            for (ix=0; ix<nx1; ix++) {
              for (iz=0; iz<nz1; iz++) {
                i = iz + nz1*ix;
                j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
                ik = iz+nbt + (ix+nbl)*nz; /* padded grid */
                curr[j] += vv[ik]*wvfld0[it/snap][i];
              }
            }
          }
        } else { /* source injection */
          t = it*dt;
          for (i=0; i<n_srcs; i++) {
            for(ix=-1;ix<=1;ix++) {
              for(iz=-1;iz<=1;iz++) {
                ik = spz[i]+iz+nz*(spx[i]+ix);
                j = spz[i]+iz+nz2*(spx[i]+ix);
                if (src==0) {
                  curr[j] += vv[ik]*rick[i][it]/(abs(ix)+abs(iz)+1);
                } else {
                  curr[j] += vv[ik]*Ricker(t, f0[i], t0[i], A[i])/(abs(ix)+abs(iz)+1);
                }
              }
            }
          }
        }

	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}

        /* record data */
        if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix = 0; ix < gpl; ix++) {
                dat[ix][it] = curr[gpz+(ix+gpx)*nz2];
            }
        }
        if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
            for (iz = 0; iz < gpl_v; iz++) {
                dat_v[iz][it] = curr[gpz_v+iz+(gpx_v)*nz2];
            }
        }

	/* save wavefield */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(prev);
    free(cwave);
    free(cwavem);
    
    return 0;
}

int psp5(int split, float **wvfld1, float **wvfld, float **dat, float *img, float *vel, pspar par)
/*< pseudo-spectral back propagation of multiple receiver wavefields and cross-correlation >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j,ir;  /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k;
    float c, old;
    int d_gpl, o_gpl, n_gpl;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,**currs,*prev,*lapl;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
    src   = par->src;
    nt    = par->nt;
    dt    = par->dt;
    f0    = par->f0;
    t0    = par->t0;
    A     = par->A;
    verb  = par->verb;
    ps    = par->ps;
    vref  = par->vref;
    

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
#else
    nth = 1;
#endif
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    currs  = sf_floatalloc2(nzx2,nt/snap);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (it=0; it<nt/snap; it++) {
        for (iz=0; iz < nzx2; iz++) {
            currs[it][iz] = 1.;
        }
    }
    for (iz=0; iz < nz1*nx1; iz++) {
        img[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    /* step backward in time */
    it1 = nt-1;
    it2 = -1;
    its = -1;	

    d_gpl = (int) ceil(gpl*1.0/split);
    /* Loop over groups of receivers */
    for (ir=0; ir<split; ir++) {
        o_gpl = ir*d_gpl;
        n_gpl = (ir==split-1) ? gpl-ir*d_gpl : d_gpl;

        sf_warning("ir=%d/%d",ir,split);
        sf_warning("gpl=%d,split=%d,o_gpl=%d,n_gpl=%d",gpl,split,o_gpl,n_gpl);

        for (iz=0; iz < nzx2; iz++) {
            curr[iz] = 0.;
            prev[iz] = 0.;
        }
        /* MAIN LOOP */
        for (it=it1; it!=it2; it+=its) {

            if(verb) sf_warning("it=%d/%d;",it,nt);

            /* first receiver wavefield */
            /* matrix multiplication */
            fft2(curr,cwave);

            for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                cwavem[ik] = cwave[ik]*lapl[ik];
#else
                cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
            }

            ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
            for (ix = 0; ix < nx; ix++) {
                for (iz=0; iz < nz; iz++) {
                    i = iz+ix*nz;  /* original grid */
                    j = iz+ix*nz2; /* padded grid */

                    old = c = curr[j];
                    c += c - prev[j];
                    prev[j] = old;
                    c += wave[j]*vv[i];
                    curr[j] = c;
                    if (snap > 0 && it%snap==0) {
                        currs[it/snap][j] *= c;
                    }
                }
            }

            /* inject data */
            if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
                for (ix = 0; ix < n_gpl; ix++) {
                    curr[gpz+(ix+o_gpl+gpx)*nz2] += vv[gpz+(ix+o_gpl+gpx)*nz]*dat[ix+o_gpl][it];
                }
            }

            /*apply abc*/
            if (abc) {
                abc_apply(curr);
                abc_apply(prev);
            }

            if (ir == split-1) {
                /* cross-correlation imaging condition */
                if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
                    for (ix = 0; ix < nx1; ix++) {
                        for (iz = 0; iz < nz1; iz++) {
                            i = iz + nz1*ix;
                            j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
                            img[i] += wvfld1[it/snap][i]*currs[it/snap][j];
                        }
                    }
                }

                /* save wavefield */
                if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
                    for (ix=0; ix<nx1; ix++) {
                        for (iz=0; iz<nz1; iz++) {
                            i = iz + nz1*ix;
                            j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
                            wvfld[it/snap][i] = currs[it/snap][j];
                        }
                    }
                }
            }
        }
        if(verb) sf_warning(".");
    }

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(*currs); free(currs);
    free(prev);
    free(cwave);
    free(cwavem);
    
    return 0;
}



int main(int argc, char* argv[])
{

    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt,ntsnap;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps, mig;
    float vref;

    pspar par;
    int nx1, nz1; /*domain of interest*/
    int it;
    float *vel,**dat,**dat_v,**wvfld,*img; /*velocity profile*/
    sf_file Fi,Fo,Fd,Fd_v,snaps; /* I/O files */
    sf_axis az,ax; /* cube axes */

    sf_init(argc,argv);

    if (!sf_getint("snap",&snap)) snap=0; /* interval for snapshots */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex fft */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if(!sf_getbool("abc",&abc)) abc=false; /* absorbing flag */
    if (abc) {
      if(!sf_getint("nbt",&nbt)) sf_error("Need nbt!");
      if(!sf_getint("nbb",&nbb)) nbb = nbt;
      if(!sf_getint("nbl",&nbl)) nbl = nbt;
      if(!sf_getint("nbr",&nbr)) nbr = nbt;
      if(!sf_getfloat("ct",&ct)) sf_error("Need ct!");
      if(!sf_getfloat("cb",&cb)) cb = ct;
      if(!sf_getfloat("cl",&cl)) cl = ct;
      if(!sf_getfloat("cr",&cr)) cr = ct;
    } else {
      nbt = 0; nbb = 0; nbl = 0; nbr = 0;
      ct = 0; cb = 0; cl = 0; cr = 0;
    }
    if (!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if (!sf_getbool("ps",&ps)) ps=false; /* use pseudo-spectral */
    if (ps) sf_warning("Using pseudo-spectral...");
    else sf_warning("Using pseudo-analytical...");
    if (!sf_getbool("mig",&mig)) mig=false; /* use pseudo-spectral */
    if (mig) sf_warning("Time-reversal propagation");
    else sf_warning("Forward modeling");
    if (!sf_getfloat("vref",&vref)) vref=1500; /* reference velocity (default using water) */

    /* setup I/O files */
    Fi = sf_input ("in");
    Fo = sf_output("out");
    if (mig) {
      gpl = -1;
      gpl_v = -1;
      if (NULL==sf_getstring("dat") && NULL==sf_getstring("dat_v"))
	sf_error("Need Data!");
      if (NULL!=sf_getstring("dat")) {
	Fd = sf_input("dat");
	sf_histint(Fd,"n1",&nt);
	sf_histfloat(Fd,"d1",&dt);
	sf_histint(Fd,"n2",&gpl);
      } else Fd = NULL;
      if (NULL!=sf_getstring("dat_v")) {
	Fd_v = sf_input("dat_v");
	sf_histint(Fd_v,"n1",&nt);
	sf_histfloat(Fd_v,"d1",&dt);
	sf_histint(Fd_v,"n2",&gpl_v);
      } else Fd_v = NULL;
      src = -1; n_srcs = -1;
      spx = NULL; spz = NULL;
      f0 = NULL; t0 = NULL; A = NULL;
    } else {
      Fd = NULL;
      if (!sf_getint("nt",&nt)) sf_error("Need nt!");
      if (!sf_getfloat("dt",&dt)) sf_error("Need dt!");
      if (!sf_getint("gpl",&gpl)) gpl = -1; /* geophone length */
      if (!sf_getint("gpl_v",&gpl_v)) gpl_v = -1; /* geophone height */
      if (!sf_getint("src",&src)) src=0; /* source type */
      if (!sf_getint("n_srcs",&n_srcs)) n_srcs=1; /* source type */
      spx = sf_intalloc(n_srcs);
      spz = sf_intalloc(n_srcs);
      f0  = sf_floatalloc(n_srcs);
      t0  = sf_floatalloc(n_srcs);
      A   = sf_floatalloc(n_srcs);
      if (!sf_getints("spx",spx,n_srcs)) sf_error("Need spx!"); /* shot position x */
      if (!sf_getints("spz",spz,n_srcs)) sf_error("Need spz!"); /* shot position z */
      if (!sf_getfloats("f0",f0,n_srcs)) sf_error("Need f0! (e.g. 30Hz)");   /*  wavelet peak freq */
      if (!sf_getfloats("t0",t0,n_srcs)) sf_error("Need t0! (e.g. 0.04s)");  /*  wavelet time lag */
      if (!sf_getfloats("A",A,n_srcs)) sf_error("Need A! (e.g. 1)");     /*  wavelet amplitude */
    }
    if (!sf_getint("gpx",&gpx)) gpx = -1; /* geophone position x */
    if (!sf_getint("gpz",&gpz)) gpz = -1; /* geophone position z */
    if (!sf_getint("gpx_v",&gpx_v)) gpx_v = -1; /* geophone position x */
    if (!sf_getint("gpz_v",&gpz_v)) gpz_v = -1; /* geophone position z */

    if (SF_FLOAT != sf_gettype(Fi)) sf_error("Need float input");

    /* Read/Write axes */
    az = sf_iaxa(Fi,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fi,2); nx = sf_n(ax); dx = sf_d(ax);
    nz1 = nz-nbt-nbb;
    nx1 = nx-nbl-nbr;
    if (gpx==-1) gpx = nbl;
    if (gpz==-1) gpz = nbt;
    if (gpl==-1) gpl = nx1;
    if (gpx_v==-1) gpx_v = nbl;
    if (gpz_v==-1) gpz_v = nbt;
    if (gpl_v==-1) gpl_v = nz1;
    ntsnap=0;
    if (snap)
        for (it=0;it<nt;it++)
            if (it%snap==0) ntsnap++;
    if (mig) { /*output final wavefield*/
      sf_setn(az,nz1);
      sf_setn(ax,nx1);
      sf_oaxa(Fo,az,1);
      sf_oaxa(Fo,ax,2);
      sf_settype(Fo,SF_FLOAT);
    } else { /*output data*/
      sf_setn(ax,gpl);
      /*output horizontal data is mandatory*/
      sf_putint(Fo,"n1",nt);
      sf_putfloat(Fo,"d1",dt);
      sf_putfloat(Fo,"o1",0.);
      sf_putstring(Fo,"label1","Time");
      sf_putstring(Fo,"unit1","s");
      sf_oaxa(Fo,ax,2);
      sf_settype(Fo,SF_FLOAT);
      /*output vertical data is optional*/
      if (NULL!=sf_getstring("dat_v")) {
	Fd_v = sf_output("dat_v");
	sf_setn(az,gpl_v);
	sf_putint(Fd_v,"n1",nt);
	sf_putfloat(Fd_v,"d1",dt);
	sf_putfloat(Fd_v,"o1",0.);
	sf_putstring(Fd_v,"label1","Time");
	sf_putstring(Fd_v,"unit1","s");
	sf_oaxa(Fd_v,az,2);
	sf_settype(Fd_v,SF_FLOAT);	
      } else Fd_v = NULL;
    }

    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	sf_setn(az,nz1);
	sf_setn(ax,nx1);
	sf_oaxa(snaps,az,1);
	sf_oaxa(snaps,ax,2);
	sf_putint(snaps,"n3",ntsnap);
	sf_putfloat(snaps,"d3",dt*snap);
	sf_putfloat(snaps,"o3",0.);
	sf_putstring(snaps,"label3","Time");
	sf_putstring(snaps,"unit3","s");
    } else snaps = NULL;

    par = (pspar) sf_alloc(1,sizeof(*par));
    vel = sf_floatalloc(nz*nx);


    if (mig && NULL==Fd) {dat = NULL;  }
    else { dat = sf_floatalloc2(nt,gpl);}


    if (NULL!=Fd_v) dat_v = sf_floatalloc2(nt,gpl_v);
    else dat_v = NULL;

    if (mig) img = sf_floatalloc(nz1*nx1);
    else img = NULL;

    if (snap>0) wvfld = sf_floatalloc2(nx1*nz1,ntsnap);
    else wvfld = NULL;
    

    sf_floatread(vel,nz*nx,Fi);

    if (mig) {
      if (NULL!=Fd)   sf_floatread(dat[0],gpl*nt,Fd);
      if (NULL!=Fd_v) sf_floatread(dat_v[0],gpl_v*nt,Fd_v);
    }

    /*passing the parameters*/
    par->nx    = nx;  
    par->nz    = nz;
    par->dx    = dx;
    par->dz    = dz;
    par->n_srcs= n_srcs;
    par->spx   = spx;
    par->spz   = spz;
    par->gpz   = gpz;
    par->gpx   = gpx;
    par->gpl   = gpl;
    par->gpz_v = gpz_v;
    par->gpx_v = gpx_v;
    par->gpl_v = gpl_v;
    par->snap  = snap;
    par->cmplx = cmplx;
    par->pad1  = pad1;
    par->abc   = abc;
    par->nbt   = nbt;
    par->nbb   = nbb;
    par->nbl   = nbl;
    par->nbr   = nbr;
    par->ct    = ct;
    par->cb    = cb;
    par->cl    = cl;
    par->cr    = cr;
    par->src   = src;
    par->nt    = nt;
    par->dt    = dt;
    par->f0    = f0;
    par->t0    = t0;
    par->A     = A;
    par->verb  = verb;
    par->ps    = ps;
    par->vref  = vref;


    /*do the work*/
    psp(wvfld, dat, dat_v, img, vel, par, mig);

    if (mig) {
      sf_floatwrite(img,nz1*nx1,Fo);
    } else {
      sf_floatwrite(dat[0],gpl*nt,Fo);
      if (NULL!=Fd_v)
	sf_floatwrite(dat_v[0],gpl_v*nt,Fd_v);
    }

    if (snap>0)
      sf_floatwrite(wvfld[0],nz1*nx1*ntsnap,snaps);
    
    exit (0);
}
