#include "seismic.h"
#include "ewem.h"
void ewem(float **ux, float **uy, float **uz,
		float **mpp, float **mps1, float **mps2,
		float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		float sx,float sy,
		int nz, float oz, float dz, float gz, float sz,
		float **vp, float **vs, 
		float fmin, float fmax,
		int padt, int padx,
		float damping,
		bool adj, bool pade_flag, bool verbose)
/*< elastic wave equation depth migration operator. >*/
{
	int iz,ix,imx,imy,igx,igy,ik,iw,it,nw,nkx,nky,ntfft;
	float dw,dkx,dky;
	int ifmin,ifmax;
	float *d_t;
	complex *d_w;
	complex **ux_g_wx,**uy_g_wx,**uz_g_wx,**u_s_wx;
	fftwf_complex *a,*b;
	int *n;
	fftwf_plan p1,p2;
	float *po_p,**pd_p;
	float *po_s,**pd_s;
	float progress;
	int ithread,nthread;
	float sigma;
	float **mpp_threads,**mps1_threads,**mps2_threads;
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) mpp[ix][iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) mps1[ix][iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) mps2[ix][iz] = 0.;
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) ux[ix][it] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) uy[ix][it] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) uz[ix][it] = 0.;
	}
	ntfft = (int) 2*truncf(padt*((float) nt)/2);
	nw = (int) truncf(ntfft/2)+1;
	nkx = nmx > 1 ? padx*nmx : 1;
	nky = nmy > 1 ? padx*nmy : 1;
	dkx = 2*PI/((float) nkx)/dmx;
	dky = 2*PI/((float) nky)/dmy;
	dw = 2*PI/((float) ntfft)/dt;
	if(fmax*dt*ntfft+1<nw) ifmax = trunc(fmax*dt*ntfft)+1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = trunc(fmin*dt*ntfft);
	else ifmin = 0;
	ux_g_wx = alloc2complex(nw,nmx*nmy);
	uy_g_wx = alloc2complex(nw,nmx*nmy);
	uz_g_wx = alloc2complex(nw,nmx*nmy);
	u_s_wx = alloc2complex(nw,nmx*nmy);
	d_t = alloc1float(nt);
	d_w = alloc1complex(nw);
	for (it=0;it<nt;it++)  d_t[it] = 0.;  
	for (iw=0;iw<nw;iw++)  d_w[iw] = 0.; 

	/* decompose slowness into layer average, and layer purturbation */
	po_p = alloc1float(nz); 
	pd_p = alloc2float(nz,nmx*nmy); 
	for (iz=0;iz<nz;iz++){
		po_p[iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) po_p[iz] += vp[ix][iz];
		po_p[iz] /= (float) nmx*nmy;
		po_p[iz]  = 1./po_p[iz];
		for (ix=0;ix<nmx*nmy;ix++) pd_p[ix][iz] = 1.0/vp[ix][iz] - po_p[iz];
	}
	po_s = alloc1float(nz); 
	pd_s = alloc2float(nz,nmx*nmy); 
	for (iz=0;iz<nz;iz++){
		po_s[iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) po_s[iz] += vs[ix][iz];
		po_s[iz] /= (float) nmx*nmy;
		po_s[iz]  = 1./po_s[iz];
		for (ix=0;ix<nmx*nmy;ix++) pd_s[ix][iz] = 1.0/vs[ix][iz] - po_s[iz];
	}

	/* set up fftw plans and pass them to the OMP region of the code */
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	n = alloc1int(2); 
	n[0] = nkx;
	n[1] = nky;
	p1 = fftwf_plan_dft(2, n, a, a, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftwf_plan_dft(2, n, b, b, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (ik=0;ik<nkx*nky;ik++){
		a[ik] = 0.;
		b[ik] = 0.;
	} 
	fftwf_execute_dft(p1,a,a);
	fftwf_execute_dft(p2,b,b);

	/**********************************************************************/
	igx = (int) truncf((sx - omx)/dmx); /*position to inject source in x-dir*/
	igy = (int) truncf((sy - omy)/dmy); /*position to inject source in y-dir*/

	//fprintf(stderr,"adj=%d igx=%d\n",adj,igx);
	//fprintf(stderr,"adj=%d igy=%d\n",adj,igy);

	/* source wavefield*/
	for (ix=0;ix<nmx*nmy;ix++) for (iw=0;iw<nw;iw++) u_s_wx[ix][iw] = 0.;
	for (it=0;it<nt;it++) d_t[it] = wav[0][it];
	f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
	for (iw=0;iw<nw;iw++) u_s_wx[igx*nmy + igy][iw] = d_w[iw];

	/* receiver wavefield*/
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++){
			// x component
			for (it=0;it<nt;it++) d_t[it] = ux[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) ux_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) ux_g_wx[ix][iw] = d_w[iw];
			for (iw=ifmax;iw<nw;iw++) ux_g_wx[ix][iw] = 0.;
			// y component
			for (it=0;it<nt;it++) d_t[it] = uy[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) uy_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) uy_g_wx[ix][iw] = d_w[iw];
			for (iw=ifmax;iw<nw;iw++) uy_g_wx[ix][iw] = 0.;
			// z component
			for (it=0;it<nt;it++) d_t[it] = uz[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) uz_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) uz_g_wx[ix][iw] = d_w[iw];
			for (iw=ifmax;iw<nw;iw++) uz_g_wx[ix][iw] = 0.;
		}
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++){
			for (iw=0;iw<nw;iw++){
				ux_g_wx[ix][iw] = 0.;
				uy_g_wx[ix][iw] = 0.;
				uz_g_wx[ix][iw] = 0.;
			}
		}
	}

	sigma = 0.;
	for (it=0;it<nt;it++) if (sigma < fabsf(wav[0][it])) sigma = fabsf(wav[0][it]);
	sigma *= sqrtf(damping)/sqrtf((float) ntfft);
	sigma *= sigma;

	nthread = omp_thread_count();

	//fprintf(stderr,"nthread=%d\n",nthread);

	if (adj){
		mpp_threads = alloc2float(nz,nmx*nmy*nthread);
		mps1_threads = alloc2float(nz,nmx*nmy*nthread);
		mps2_threads = alloc2float(nz,nmx*nmy*nthread);
		for (imx=0;imx<nmx;imx++){
			for (imy=0;imy<nmy;imy++){
				for (ithread=0;ithread<nthread;ithread++){
					for (iz=0;iz<nz;iz++){
						mpp_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
						mps1_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
						mps2_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
					}
				}
			}
		}
	}
	else{
		mpp_threads = alloc2float(nz,nmx*nmy);
		mps1_threads = alloc2float(nz,nmx*nmy);
		mps2_threads = alloc2float(nz,nmx*nmy);
		for (imx=0;imx<nmx;imx++){
			for (imy=0;imy<nmy;imy++){
				for (iz=0;iz<nz;iz++){
					mpp_threads[imx*nmy + imy][iz] = mpp[imx*nmy + imy][iz];
					mps1_threads[imx*nmy + imy][iz] = mps1[imx*nmy + imy][iz];
					mps2_threads[imx*nmy + imy][iz] = mps2[imx*nmy + imy][iz];
				}
			}
		}
	}

	progress = 0.;
#pragma omp parallel for private(iw) shared(mpp_threads,mps1_threads,mps2_threads,ux_g_wx,uy_g_wx,uz_g_wx,u_s_wx)
	for (iw=ifmin;iw<ifmax;iw++){ 
		progress += 1./((float) ifmax - ifmin);
		if (verbose) progress_msg(progress);
		elastic_extrap1f(mpp_threads,mps1_threads,mps2_threads,
				ux_g_wx,uy_g_wx,uz_g_wx,
				u_s_wx,
				sigma,iw,nw,ifmax,ntfft,
				dw,dkx,dky,nkx,nky,
				nz,oz,dz,gz,sz,nmx,omx,dmx,nmy,omy,dmy,
				nthread,
				vp,po_p,pd_p,vs,po_s,pd_s,
				p1,p2,adj,pade_flag,verbose);
	}
	if (adj){
		// reduction over parallel axis 
		for (imx=0;imx<nmx;imx++){ 
			for (imy=0;imy<nmy;imy++){
				for (ithread=0;ithread<nthread;ithread++){
					for (iz=0;iz<nz;iz++){
						mpp[imx*nmy + imy][iz] += mpp_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
						mps1[imx*nmy + imy][iz] += mps1_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
						mps2[imx*nmy + imy][iz] += mps2_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
					}
				}
			}
		}
	}
	if (verbose) fprintf(stderr,"\n");

	else{
		for (ix=0;ix<nmx*nmy;ix++){
			// x component
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) d_w[iw] = ux_g_wx[ix][iw];
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) ux[ix][it] = d_t[it];
			// y component
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) d_w[iw] = uy_g_wx[ix][iw];
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) uy[ix][it] = d_t[it];
			// x component
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) d_w[iw] = uz_g_wx[ix][iw];
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) uz[ix][it] = d_t[it];
		}
	}

	free1int(n); 
	fftwf_free(a);
	fftwf_free(b);
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);
	free1float(d_t);
	free1complex(d_w);
	free2complex(ux_g_wx);
	free2complex(uy_g_wx);
	free2complex(uz_g_wx);
	free2complex(u_s_wx);
	free1float(po_p);
	free2float(pd_p);
	free1float(po_s);
	free2float(pd_s);
	free2float(mpp_threads);
	free2float(mps1_threads);
	free2float(mps2_threads);

	return;
} 

void elastic_extrap1f(float **mpp, float **mps1, float **mps2,
		complex **ux_g_wx, complex **uy_g_wx, complex **uz_g_wx, 
		complex **u_s_wx,
		float sigma, int iw, int nw,int ifmax,int ntfft,float dw,float dkx,float dky,int nkx,int nky,
		int nz, float oz, float dz, float gz, float sz,
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		int nthread,
		float **vp,float *po_p,float **pd_p,
		float **vs,float *po_s,float **pd_s,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, bool pade_flag, bool verbose)
/*< extrapolate 1 frequency >*/
{
	float w,factor,z;
	int iz,ix,imx,imy,ithread;
	complex *ux_xg,*uy_xg,*uz_xg;
	complex *up_xg,*us1_xg,*us2_xg;
	complex *up_xs;
	complex **smig;

	ithread = omp_get_thread_num(); 
	//fprintf(stderr,"ithread=%d\n",ithread);

	ux_xg = alloc1complex(nmx*nmy);
	uy_xg = alloc1complex(nmx*nmy);
	uz_xg = alloc1complex(nmx*nmy);
	up_xg = alloc1complex(nmx*nmy);
	us1_xg = alloc1complex(nmx*nmy);
	us2_xg = alloc1complex(nmx*nmy);
	up_xs = alloc1complex(nmx*nmy);

	for (ix=0;ix<nmx*nmy;ix++) ux_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) uy_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) uz_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) up_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) us1_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) us2_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) up_xs[ix] = 0.;
	if (iw==0) factor = 1.;
	else factor = 2.;

	w = iw*dw;
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++){ 
			up_xs[ix] = u_s_wx[ix][iw]/sqrtf((float) ntfft);
			ux_xg[ix] = ux_g_wx[ix][iw]/sqrtf((float) ntfft);
			uy_xg[ix] = uy_g_wx[ix][iw]/sqrtf((float) ntfft);
			uz_xg[ix] = uz_g_wx[ix][iw]/sqrtf((float) ntfft);
		}
		for (iz=0;iz<nz;iz++){ // extrapolate source and receiver wavefields
			z = oz + dz*iz;
			if (z >= sz){
				ssop(up_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,vp,po_p,pd_p,p1,p2,true,pade_flag,true,verbose);
			} 
			if (z >= gz){
				// elastic_separate(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,3); //3
				elastic_separate_2d(ux_xg,uz_xg,up_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,1./po_p[iz],1./po_s[iz],p1,p2,3); //3
				ssop(up_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,p1,p2,true,pade_flag,false,verbose);
				ssop(us1_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,p1,p2,true,pade_flag,false,verbose);
				ssop(us2_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,p1,p2,true,pade_flag,false,verbose);
				// elastic_separate(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,4); //4
				elastic_separate_2d(ux_xg,uz_xg,up_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,1./po_p[iz],1./po_s[iz],p1,p2,4); //4
				for (imx=0;imx<nmx;imx++){ 
					for (imy=0;imy<nmy;imy++){
						mpp[imx*nmy*nthread + imy*nthread + ithread][iz]  += factor*crealf(conjf(up_xs[imx*nmy + imy])*up_xg[imx*nmy + imy])/cabsf((up_xs[imx*nmy + imy]*conjf(up_xs[imx*nmy + imy])) + sigma);
						mps1[imx*nmy*nthread + imy*nthread + ithread][iz] += factor*crealf(conjf(up_xs[imx*nmy + imy])*us1_xg[imx*nmy + imy])/cabsf((up_xs[imx*nmy + imy]*conjf(up_xs[imx*nmy + imy])) + sigma);
						mps2[imx*nmy*nthread + imy*nthread + ithread][iz] += factor*crealf(conjf(up_xs[imx*nmy + imy])*us2_xg[imx*nmy + imy])/cabsf((up_xs[imx*nmy + imy]*conjf(up_xs[imx*nmy + imy])) + sigma);
					}
				}
			}
		}
	}

	else{
		smig = alloc2complex(nz,nmx*nmy);
		for (ix=0;ix<nmx*nmy;ix++) up_xs[ix] = u_s_wx[ix][iw]/sqrtf((float) ntfft);
		for (iz=0;iz<nz;iz++){ // extrapolate source wavefield 
			z = oz + dz*iz;
			if (z >= sz){
				ssop(up_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,vp,po_p,pd_p,p1,p2,true,pade_flag,true,verbose);
				for (ix=0;ix<nmx*nmy;ix++) smig[ix][iz]  = up_xs[ix]/cabsf((up_xs[ix]*conjf(up_xs[ix])) + sigma);
			}
			else{
				for (ix=0;ix<nmx*nmy;ix++) smig[ix][iz] = 0.;
			}
		}
		for (iz=nz-1;iz>=0;iz--){ // extrapolate receiver wavefield 
			z = oz + dz*iz;
			if (z >= gz){
				// elastic_separate(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,2); //2 
				elastic_separate_2d(ux_xg,uz_xg,up_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,1./po_p[iz],1./po_s[iz],p1,p2,2); //2 
				for (ix=0;ix<nmx*nmy;ix++) up_xg[ix]  += smig[ix][iz]*mpp[ix][iz];
				for (ix=0;ix<nmx*nmy;ix++) us1_xg[ix] += smig[ix][iz]*mps1[ix][iz];
				for (ix=0;ix<nmx*nmy;ix++) us2_xg[ix] += smig[ix][iz]*mps2[ix][iz];
				ssop(up_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,vp,po_p,pd_p,p1,p2,false,pade_flag,false,verbose);
				ssop(us1_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,vs,po_s,pd_s,p1,p2,false,pade_flag,false,verbose);
				ssop(us2_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,vs,po_s,pd_s,p1,p2,false,pade_flag,false,verbose);
				// elastic_separate(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,1); //1
				elastic_separate_2d(ux_xg,uz_xg,up_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,1./po_p[iz],1./po_s[iz],p1,p2,1); //1
			}
		}
		for (ix=0;ix<nmx*nmy;ix++){
			ux_g_wx[ix][iw] = ux_xg[ix]/sqrtf((float) ntfft);
			uy_g_wx[ix][iw] = uy_xg[ix]/sqrtf((float) ntfft);
			uz_g_wx[ix][iw] = uz_xg[ix]/sqrtf((float) ntfft);
		}
		free2complex(smig);
	}

	free1complex(ux_xg);
	free1complex(uy_xg);
	free1complex(uz_xg);
	free1complex(up_xg);
	free1complex(us1_xg);
	free1complex(us2_xg);
	free1complex(up_xs);

	return;
}

void ssop(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
		float **v,float *po,float **pd,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool pade_flag, 
		bool src,
		bool verbose)
{

	float kx,ky,s;
	complex L;
	int ik,ikx,iky,imx,imy; 
	complex *d_k;
	fftwf_complex *a,*b;
	int lmx,lmy;

	if (nmx>100) lmx=30;
	else lmx=0;
	if (nmy>100) lmy=30;
	else lmy=0;

	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	d_k = alloc1complex(nkx*nky);

	if (adj){
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) a[imx*nky + imy] = d_x[imx*nmy + imy];
				else a[imx*nky + imy] = 0.;
			}
		}
	}
	else{
		boundary_condition(d_x,nmx,lmx,nmy,lmy);    
		if (pade_flag) pade(d_x,nmx,omx,dmx,nmy,omy,dmy,dz,w,iz,v,po,pd,adj,src,verbose);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){
					L = cexpf(I*w*pd[imx*nmy + imy][iz]*dz);  
					a[imx*nky + imy] = d_x[imx*nmy + imy]*L; // SS operator
				}
				else a[imx*nky + imy] = 0.;
			}
		}
	}

	fftwf_execute_dft(p1,a,a); 
	for (ikx=0;ikx<nkx;ikx++){
		if (ikx<= (int) truncf(nkx/2)) kx = (float) dkx*ikx;
		else                          kx = -((float) dkx*nkx - dkx*ikx);
		for (iky=0;iky<nky;iky++){
			if (iky<= (int) truncf(nky/2)) ky = (float) dky*iky;
			else                          ky = -((float) dky*nky - dky*iky);
			s = (w*w)*(po[iz]*po[iz]) - (kx*kx) - (ky*ky);
			if (s>=0) L = cexpf(I*sqrtf(s)*dz); 
			else L = cexpf(-0.2*sqrtf(fabsf(s))*fabsf(dz));
			d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);        
		}
	}
	for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
	fftwf_execute_dft(p2,b,b);
	if (adj){
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){
					L = cexpf(I*w*pd[imx*nmy + imy][iz]*dz);
					d_x[imx*nmy + imy] = ((complex) b[imx*nky + imy])*L/sqrtf((float) nkx*nky); // SS operator
				}
			}
		}

		if (pade_flag) pade(d_x,nmx,omx,dmx,nmy,omy,dmy,dz,w,iz,v,po,pd,adj,src,verbose);
		boundary_condition(d_x,nmx,lmx,nmy,lmy);    

	}
	else{
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){ 
					d_x[imx*nmy + imy] = ((complex) b[imx*nky + imy])/sqrtf((float) nkx*nky);
				}
			}
		}
	}

	free1complex(d_k);
	fftwf_free(a);
	fftwf_free(b);


	return;
}

void pade(complex *d,
		int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,
		float w, int iz,
		float **v,float *po,float **pd,
		bool adj, bool src,
		bool verbose)
/*  Split-Step Padé Fourier migration correction: Huang and Fehler (2002) */
{

	int imx,imy;
	float m,ko,w2;
	complex *aX,*X,*aY,*Y;

	if (w>60.) w2 = w;
	else w2 = 60.;
	X = alloc1complex(nmx); 
	aX = alloc1complex(nmx); 
	for (imy=0;imy<nmy;imy++){
		for (imx=0;imx<nmx;imx++){ 
			X[imx] = d[imx*nmy + imy];
			m = v[imx*nmy + imy][iz]*po[iz];
			ko = w*po[iz];
			if (src) aX[imx] =-(v[imx*nmy + imy][iz]*v[imx*nmy + imy][iz]/(w2*w2))*I*(m-1.)*ko*fabsf(dz)/(2.*m*m);
			else     aX[imx] = (v[imx*nmy + imy][iz]*v[imx*nmy + imy][iz]/(w2*w2))*I*(m-1.)*ko*fabsf(dz)/(2.*m*m);
		}
		fdop(X,nmx,dmx,aX,adj);
		for (imx=0;imx<nmx;imx++){
			//fprintf(stderr,"a=%f + %fi\n",crealf(a),cimagf(a));
			if (adj) d[imx*nmy + imy] += X[imx];
			else     d[imx*nmy + imy] -= X[imx];
		}
	}
	free1complex(aX);
	free1complex(X);

	Y = alloc1complex(nmy); 
	aY = alloc1complex(nmy); 
	for (imx=0;imx<nmx;imx++){
		for (imy=0;imy<nmy;imy++){ 
			Y[imy] = d[imx*nmy + imy];
			m = v[imx*nmy + imy][iz]*po[iz];
			ko = w*po[iz];
			if (src) aY[imy] =-(v[imx*nmy + imy][iz]*v[imx*nmy + imy][iz]/(w2*w2))*I*(m-1.)*ko*fabsf(dz)/(2.*m*m);
			else     aY[imy] = (v[imx*nmy + imy][iz]*v[imx*nmy + imy][iz]/(w2*w2))*I*(m-1.)*ko*fabsf(dz)/(2.*m*m);
		}
		fdop(Y,nmy,dmy,aY,adj);
		for (imy=0;imy<nmy;imy++){
			//fprintf(stderr,"a=%f + %fi\n",crealf(a),cimagf(a));
			if (adj) d[imx*nmy + imy] += Y[imy];
			else     d[imx*nmy + imy] -= Y[imy];
		}
	}
	free1complex(aY);
	free1complex(Y);

	return;
}

void fdop(complex *X,int nx,float dx, complex *a, bool adj)
	/* centered 2nd order finite difference approximation to the 2nd derivative */
{
	int ix;
	complex *X2;

	if (nx > 3){
		X2 = alloc1complex(nx);
		if (adj){
			X2[0] = -2.*a[0]*X[0] + a[0]*X[1];
			X2[nx-1] = -2.*a[nx-1]*X[nx-1] + a[nx-1]*X[nx-2];
			for (ix=1;ix<nx-1;ix++) X2[ix] = (a[ix]*X[ix-1] - 2.*a[ix]*X[ix] + a[ix]*X[ix+1]);
		}
		else{
			X2[0] = -2.*a[0]*X[0] + a[1]*X[1];
			X2[nx-1] = -2.*a[nx-1]*X[nx-1] + a[nx-2]*X[nx-2];
			for (ix=1;ix<nx-1;ix++) X2[ix] = (a[ix-1]*X[ix-1] - 2.*a[ix]*X[ix] + a[ix+1]*X[ix+1]);
		}
		for (ix=0;ix<nx;ix++) X[ix] = X2[ix]/dx/dx;
		free1complex(X2);
	}
	else{
		for (ix=0;ix<nx;ix++) X[ix] = 0.;
	}

	return;
}

void f_op(complex *m,float *d,int nw,int nt,bool adj)
{
	fftwf_complex *out1a,*in1b;
	float *in1a,*out1b;
	int ntfft,it,iw;
	fftwf_plan p1a,p1b;

	ntfft = (nw-1)*2;

	if (adj){ /* data --> model */
		out1a = fftwf_malloc(sizeof(fftwf_complex) * nw);
		in1a = alloc1float(ntfft);
		p1a = fftwf_plan_dft_r2c_1d(ntfft, in1a, (fftwf_complex*)out1a, FFTW_ESTIMATE);
		for(it=0;it<nt;it++) in1a[it] = d[it];
		for(it=nt;it<ntfft;it++) in1a[it] = 0.;
		fftwf_execute(p1a); 
		for(iw=0;iw<nw;iw++) m[iw] = out1a[iw];
		fftwf_destroy_plan(p1a);
		fftwf_free(in1a); fftwf_free(out1a);
	}

	else{ /* model --> data */
		out1b = alloc1float(ntfft);
		in1b = fftwf_malloc(sizeof(fftwf_complex) * ntfft);
		p1b = fftwf_plan_dft_c2r_1d(ntfft, (fftwf_complex*)in1b, out1b, FFTW_ESTIMATE);
		for(iw=0;iw<nw;iw++) in1b[iw] = m[iw];
		for(iw=nw;iw<ntfft;iw++) in1b[iw] = 0.;
		fftwf_execute(p1b); 
		for(it=0;it<nt;it++) d[it] = out1b[it];
		fftwf_destroy_plan(p1b);
		fftwf_free(in1b); fftwf_free(out1b);
	}

	return;
}

void progress_msg(float progress)
{ 
	fprintf(stderr,"\r[%6.2f%% complete]      ",progress*100);
	return;
}

float signf(float a)
	/*< sign of a float >*/
{
	float b;
	if (a>0.)      b = 1.;
	else if (a<0.) b =-1.;
	else          b = 0.;
	return b;
}

int compare (const void * a, const void * b)
{
	float fa = *(const float*) a;
	float fb = *(const float*) b;
	return (fa > fb) - (fa < fb);
}

int omp_thread_count() {
	int n = 0;
#pragma omp parallel reduction(+:n)
	n += 1;
	return n;
}

void boundary_condition(complex *d_x,int nmx,int lmx,int nmy,int lmy)
{
	int imx,imy;
	float tmx,tmy; 
	tmx = 1.;tmy = 1.;
	for (imx=0;imx<nmx;imx++){
		if (imx>=0   && imx<lmx) tmx = expf(-powf(0.015*((float) lmx - imx),2.));
		if (imx>=lmx && imx<=nmx-lmx) tmx = 1.;
		if (imx>nmx-lmx && imx<nmx) tmx = expf(-powf(0.015*((float) imx - nmx + lmx),2.));
		for (imy=0;imy<nmy;imy++){
			if (imy>=0   && imy<lmy) tmy = expf(-powf(0.015*((float) lmy - imy),2.));
			if (imy>=lmy && imy<=nmy-lmy) tmy = 1.;
			if (imy>nmy-lmy && imy<nmy) tmy = expf(-powf(0.015*((float) imy - nmy + lmy),2.));
			d_x[imx*nmy + imy] *= tmx*tmy;
		}
	}  

	return;
}

void elastic_separate_2d(complex *ux, complex *uz,
		complex *up, complex *us2,
		float w, 
		float dkx, int nkx, int nmx, float omx, float dmx,
		float vp,float vs,
		fftwf_plan p1,fftwf_plan p2,
		int option)
{
	int imx,ikx,ik;
	fftwf_complex *a,*b;  
	complex *ux_k,*uz_k,*up_k,*us2_k;
	float kx,s1,s2,norm;
	complex kzp,kzs;
	ux_k = alloc1complex(nkx);
	uz_k = alloc1complex(nkx);
	up_k = alloc1complex(nkx);
	us2_k = alloc1complex(nkx);
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	if (option==2 || option==3){ /* inverse: wavefield components to wavefield potentials */

		// x-component
		for(imx=0;imx<nkx;imx++){ 
			a[imx] = imx < nmx ? ux[imx] : 0.;
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx;ik++) ux_k[ik] = a[ik]/sqrtf(nkx);

		// z-component
		for(imx=0;imx<nkx;imx++){ 
			a[imx] = imx < nmx ? uz[imx] : 0.;
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx;ik++) uz_k[ik] = a[ik]/sqrtf(nkx); 

		for (ikx=0;ikx<nkx;ikx++){
				if (ikx<nkx/2.) kx = dkx*ikx;
				else            kx = -(dkx*nkx - dkx*ikx);
				s1 = w*w/(vp*vp) - kx*kx;
				s2 = w*w/(vs*vs) - kx*kx;
				if (s1>0.) kzp = sqrtf(s1);
				else kzp = I*sqrtf(-s1);
				if (s2>0.) kzs = sqrtf(s2);
				else kzs = I*sqrtf(-s2);
				norm = kx*kx + kzp*kzs;
				if (option == 2){   // Q-inv
					up_k[ikx]  =  I *kx*ux_k[ikx] + I*kzs*uz_k[ikx];  
					us2_k[ikx] = -I*kzp*ux_k[ikx] + I*kx*uz_k[ikx];
				}  
				else{ // option 3  // (Q)*
					up_k[ikx] = (I*kx*ux_k[ikx] + I*conjf(kzp)*uz_k[ikx])/(norm + 1e-10);  
					us2_k[ikx] = (-I*conjf(kzs)*ux_k[ikx] + I*kx*uz_k[ikx])/(norm + 1e-10);
				}  
		}
		for (ik=0;ik<nkx;ik++) b[ik] = up_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx) up[imx] = b[imx]/sqrtf(nkx);
		}      
		for (ik=0;ik<nkx;ik++) b[ik] = us2_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx) us2[imx] = b[imx]/sqrtf(nkx);
		}      
	}
	else if (option==1 || option==4){ /* forward: wavefield potentials to wavefield components */
		// p-component
		for(imx=0;imx<nkx;imx++){ 
			a[imx] = (imx < nmx) ? up[imx] : 0.;
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx;ik++) up_k[ik] = a[ik]/sqrtf(nkx);
		// s2-component
		for(imx=0;imx<nkx;imx++){ 
			a[imx] = (imx < nmx) ? us2[imx] : 0.;
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx;ik++) us2_k[ik] = a[ik]/sqrtf(nkx); 
		for (ikx=0;ikx<nkx;ikx++){
				if (ikx<nkx/2.) kx = dkx*ikx;
				else         kx = -(dkx*nkx - dkx*ikx);
				s1 = w*w/(vp*vp) - kx*kx;
				s2 = w*w/(vs*vs) - kx*kx;
				if (s1>0.) kzp = sqrtf(s1);
				else kzp = I*sqrtf(-s1);
				if (s2>0.) kzs = sqrtf(s2);
				else kzs = I*sqrtf(-s2);
				norm = kx*kx + kzp*kzs;
				if (option == 1){ // Q
					ux_k[ikx] = (-I*kx*up_k[ikx] + I*kzs*us2_k[ikx])/(norm + 1e-10);  
					uz_k[ikx] = (-I*kzp*up_k[ikx] - I*kx*us2_k[ikx])/(norm + 1e-10);
				}  
				else{ // option 4 // (Q-inv)*
					ux_k[ikx] =  -I*kx*up_k[ikx] + I*conjf(kzp)*us2_k[ikx];  
					uz_k[ikx] = -I*conjf(kzs)*up_k[ikx] - I*kx*us2_k[ikx];	    
				}  
		}
		
		for (ik=0;ik<nkx;ik++) b[ik] = ux_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx) ux[imx] = b[imx]/sqrtf(nkx);
		}      

		for (ik=0;ik<nkx;ik++) b[ik] = uz_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx) uz[imx] = b[imx]/sqrtf(nkx);
		}     
		 
	}

	fftwf_free(a);
	fftwf_free(b);
	free1complex(ux_k);
	free1complex(uz_k);
	free1complex(up_k);
	free1complex(us2_k);

	return;
}

void elastic_separate_3d(complex *ux, complex *uy, complex *uz,
		complex *up, complex *us1, complex *us2,
		float w, 
		float dkx, int nkx, int nmx, float omx, float dmx,
		float dky, int nky, int nmy, float omy, float dmy,
		float vp,float vs,
		fftwf_plan p1,fftwf_plan p2,
		int option)
{
	int imx,imy,ikx,iky,ik;
	fftwf_complex *a,*b;  
	complex *ux_k,*uy_k,*uz_k,*up_k,*us1_k,*us2_k;
	float kx,ky,k,s1,s2,norm;
	complex kzp,kzs;
	ux_k = alloc1complex(nkx*nky);
	uy_k = alloc1complex(nkx*nky);
	uz_k = alloc1complex(nkx*nky);
	up_k = alloc1complex(nkx*nky);
	us1_k = alloc1complex(nkx*nky);
	us2_k = alloc1complex(nkx*nky);
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	if (option==2 || option==3){ /* inverse: wavefield components to wavefield potentials */
		// x-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? ux[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) ux_k[ik] = a[ik]/sqrtf(nkx*nky);
		// y-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? uy[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) uy_k[ik] = a[ik]/sqrtf(nkx*nky);
		// z-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? uz[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) uz_k[ik] = a[ik]/sqrtf(nkx*nky); 
		for (ikx=0;ikx<nkx;ikx++){
			for (iky=0;iky<nky;iky++){
				if (ikx<nkx/2.) kx = dkx*ikx;
				else         kx = -(dkx*nkx - dkx*ikx);
				if (iky<nky/2.) ky = dky*iky;
				else         ky = -(dky*nky - dky*iky);
				s1 = w*w/(vp*vp) - kx*kx - ky*ky;
				s2 = w*w/(vs*vs) - kx*kx - ky*ky;
				if (s1>0.) kzp = sqrtf(s1);
				else kzp = I*sqrtf(-s1);
				if (s2>0.) kzs = sqrtf(s2);
				else kzs = I*sqrtf(-s2);
				// normalize the wavenumbers
				norm = sqrtf(kx*kx + ky*ky + cabsf(kzp*kzs));
				//fprintf(stderr,"norm=%f\n",norm);
				kx = kx/(norm + 1e-10);
				ky = ky/(norm + 1e-10);
				kzp = kzp/(norm + 1e-10);
				kzs = kzs/(norm + 1e-10);
				k = sqrt(kx*kx + ky*ky);
				if (option == 2){   // Q-inv
					up_k[ikx*nky + iky]  = kx*ux_k[ikx*nky + iky] + ky*uy_k[ikx*nky + iky] + kzs*uz_k[ikx*nky + iky];  
					us1_k[ikx*nky + iky] = (k*ky - ky*kzp*kzs/(k+1e-10))*ux_k[ikx*nky + iky] + (k*kx - kx*kzp*kzs/(k+1e-10))*uy_k[ikx*nky + iky];   
					us2_k[ikx*nky + iky] = (kx*kzp/(k+1e-10))*ux_k[ikx*nky + iky] + (ky*kzp/(k+1e-10))*uy_k[ikx*nky + iky] - k*uz_k[ikx*nky + iky];
				}  
				else{ // option 3  // (Q)*
					up_k[ikx*nky + iky] = kx*ux_k[ikx*nky + iky] + ky*uy_k[ikx*nky + iky] + conjf(kzp)*uz_k[ikx*nky + iky];  
					us1_k[ikx*nky + iky] = -(ky/(k+1e-10))*ux_k[ikx*nky + iky] + (kx/(k+1e-10))*uy_k[ikx*nky + iky];  
					us2_k[ikx*nky + iky] = (kx*conjf(kzs)/(k+1e-10))*ux_k[ikx*nky + iky] + (ky*conjf(kzs)/(k+1e-10))*uy_k[ikx*nky + iky] - k*uz_k[ikx*nky + iky];  
				}  
			}
		}
		for (ik=0;ik<nkx*nky;ik++) b[ik] = up_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) up[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
		for (ik=0;ik<nkx*nky;ik++) b[ik] = us1_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) us1[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
		for (ik=0;ik<nkx*nky;ik++) b[ik] = us2_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) us2[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
	}
	else if (option==1 || option==4){ /* forward: wavefield potentials to wavefield components */
		// p-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? up[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) up_k[ik] = a[ik]/sqrtf(nkx*nky);
		// s1-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? us1[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) us1_k[ik] = a[ik]/sqrtf(nkx*nky);
		// s2-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? us2[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) us2_k[ik] = a[ik]/sqrtf(nkx*nky); 
		for (ikx=0;ikx<nkx;ikx++){
			for (iky=0;iky<nky;iky++){
				if (ikx<nkx/2.) kx = dkx*ikx;
				else         kx = -(dkx*nkx - dkx*ikx);
				if (iky<nky/2.) ky = dky*iky;
				else         ky = -(dky*nky - dky*iky);
				s1 = w*w/(vp*vp) - kx*kx - ky*ky;
				s2 = w*w/(vs*vs) - kx*kx - ky*ky;
				if (s1>0.) kzp = sqrtf(s1);
				else kzp = I*sqrtf(-s1);
				if (s2>0.) kzs = sqrtf(s2);
				else kzs = I*sqrtf(-s2);
				// normalize the wavenumbers
				norm = sqrtf(kx*kx + ky*ky + cabsf(kzp*kzs));
				//fprintf(stderr,"norm=%f\n",norm);
				kx = kx/(norm + 1e-10);
				ky = ky/(norm + 1e-10);
				kzp = kzp/(norm + 1e-10);
				kzs = kzs/(norm + 1e-10);
				k = sqrt(kx*kx + ky*ky);
				if (option == 1){ // Q
					ux_k[ikx*nky + iky] =  kx*up_k[ikx*nky + iky] - (ky/(k+1e-10))*us1_k[ikx*nky + iky] + (kx*kzs/(k+1e-10))*us2_k[ikx*nky + iky];  
					uy_k[ikx*nky + iky] =  ky*up_k[ikx*nky + iky] + (kx/(k+1e-10))*us1_k[ikx*nky + iky] + (ky*kzs/(k+1e-10))*us2_k[ikx*nky + iky];  
					uz_k[ikx*nky + iky] = kzp*up_k[ikx*nky + iky] - k*us2_k[ikx*nky + iky];
				}  
				else{ // option 4 // (Q-inv)*
					ux_k[ikx*nky + iky] =  kx*up_k[ikx*nky + iky] + (k*ky - ky*conjf(kzp)*conjf(kzs)/(k+1e-10))*us1_k[ikx*nky + iky] + (kx*conjf(kzp)/(k+1e-10))*us2_k[ikx*nky + iky];  
					uy_k[ikx*nky + iky] =  ky*up_k[ikx*nky + iky] + (k*kx - kx*conjf(kzp)*conjf(kzs)/(k+1e-10))*us1_k[ikx*nky + iky] + (ky*conjf(kzp)/(k+1e-10))*us2_k[ikx*nky + iky] ;   
					uz_k[ikx*nky + iky] = conjf(kzs)*up_k[ikx*nky + iky] - k*us2_k[ikx*nky + iky];	    
				}  
			}
		}
		for (ik=0;ik<nkx*nky;ik++) b[ik] = ux_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) ux[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
		for (ik=0;ik<nkx*nky;ik++) b[ik] = uy_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) uy[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
		for (ik=0;ik<nkx*nky;ik++) b[ik] = uz_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) uz[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
	}

	fftwf_free(a);
	fftwf_free(b);
	free1complex(ux_k);
	free1complex(uy_k);
	free1complex(uz_k);
	free1complex(up_k);
	free1complex(us1_k);
	free1complex(us2_k);

	return;
}


