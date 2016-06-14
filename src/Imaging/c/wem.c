#include "seismic.h"
#include "wem.h"
void wem(float **d, float **m, float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		float sx,float sy,
		int nz, float oz, float dz, float gz, float sz,
		float **vel, int nref, float fmin, float fmax,
		int padt, int padx,
		bool adj, bool pspi, bool verbose)
/*< wave equation depth migration operator. Can specify different velocities for src and rec side wavefields. >*/
{
	int iz,ix,imx,imy,igx,igy,ik,iw,it,nw,nkx,nky,ntfft;
	float dw,dkx,dky,w;
	int ifmin,ifmax;
	float *d_t;
	complex *d_w,**d_g_wx,**d_s_wx;
	fftwf_complex *a,*b;
	int *n;
	fftwf_plan p1,p2;
	float *po,**pd;
	float progress;
	int ithread,nthread;
	float **m_threads;

	/* decompose slowness into layer average, and layer purturbation */
	po = alloc1float(nz); 
	pd = alloc2float(nz,nmx*nmy); 
	for (iz=0;iz<nz;iz++){
		po[iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) po[iz] += vel[ix][iz];
		po[iz] /= (float) nmx*nmy;
		po[iz]  = 1./po[iz];
		for (ix=0;ix<nmx*nmy;ix++) pd[ix][iz] = 1.0/vel[ix][iz] - po[iz];
	}


	/****************************************************************************************/
	float **vref,vmin,vmax,v;
	int **iref1,**iref2;
	int iref;
	/* generate reference velocities for each depth step */
	vref = alloc2float(nz,nref); /* reference velocities for each layer */
	iref1 = alloc2int(nz,nmx*nmy); /* index of nearest lower reference velocity for each subsurface point */
	iref2 = alloc2int(nz,nmx*nmy); /* index of nearest upper reference velocity for each subsurface point */
	for (iz=0;iz<nz;iz++){
		vmin=vel[0][iz];
		for (ix=0;ix<nmx;ix++) if (vel[ix][iz] < vmin) vmin = vel[ix][iz];
		vmax=vel[nmx-1][iz];
		for (ix=0;ix<nmx*nmy;ix++) if (vel[ix][iz] > vmax) vmax = vel[ix][iz];
		for (iref=0;iref<nref;iref++) vref[iref][iz] = vmin + (float) iref*(vmax-vmin)/((float) nref-1);
		for (ix=0;ix<nmx*nmy;ix++){
			v = vel[ix][iz];
			if (vmax>vmin+10){
				iref = (int) (nref-1)*(v-vmin)/(vmax-vmin);
				iref1[ix][iz] = iref;
				iref2[ix][iz] = iref+1;
				if (iref>nref-2){
					iref1[ix][iz] = nref-1;
					iref2[ix][iz] = nref-1;
				}
			}
			else{
				iref1[ix][iz] = 0;
				iref2[ix][iz] = 0;
			}
		}

	}
	/****************************************************************************************/

	if (adj){
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) m[ix][iz] = 0.;
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) d[ix][it] = 0.;
	}
	ntfft = (int) 2*padt*((float) nt)/2;
	nw = (int) ntfft/2 + 1;
	nkx = nmx > 1 ? padx*nmx : 1;
	nky = nmy > 1 ? padx*nmy : 1;
	dkx = 2*PI/((float) nkx)/dmx;
	dky = 2*PI/((float) nky)/dmy;
	dw = 2*PI/((float) ntfft)/dt;
	if(fmax*dt*ntfft+1<nw) ifmax = (int) fmax*dt*ntfft + 1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = (int) fmin*dt*ntfft;
	else ifmin = 0;
	d_g_wx = alloc2complex(nw,nmx*nmy);
	d_s_wx = alloc2complex(nw,nmx*nmy);
	d_t = alloc1float(nt);
	d_w = alloc1complex(nw);
	for (it=0;it<nt;it++)  d_t[it] = 0.;  
	for (iw=0;iw<nw;iw++)  d_w[iw] = 0.; 

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
	igx = (int) (sx - omx)/dmx; /*position to inject source in x-dir*/
	igy = (int) (sy - omy)/dmy; /*position to inject source in y-dir*/
	/* source wavefield*/
	for (ix=0;ix<nmx*nmy;ix++) for (iw=0;iw<nw;iw++) d_s_wx[ix][iw] = 0.;
	for (it=0;it<nt;it++) d_t[it] = wav[0][it];

	f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
	for (iw=0;iw<nw;iw++) d_s_wx[igx*nmy + igy][iw] = d_w[iw];

	/* receiver wavefield*/
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++){
			for (it=0;it<nt;it++) d_t[it] = d[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) d_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++){ 
				w = iw*dw;
				d_g_wx[ix][iw] = w*w*d_w[iw];
			}
			for (iw=ifmax;iw<nw;iw++) d_g_wx[ix][iw] = 0.;
		}
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++){
			for (iw=0;iw<nw;iw++){
				d_g_wx[ix][iw] = 0.;
			}
		}
	}
	nthread = omp_thread_count();
	if (adj){
		m_threads = alloc2float(nz,nmx*nmy*nthread);
		for (imx=0;imx<nmx;imx++){
			for (imy=0;imy<nmy;imy++){
				for (ithread=0;ithread<nthread;ithread++){
					for (iz=0;iz<nz;iz++){
						m_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
					}
				}
			}
		}
	}
	else{
		m_threads = alloc2float(nz,nmx*nmy);
		for (imx=0;imx<nmx;imx++){
			for (imy=0;imy<nmy;imy++){
				for (iz=0;iz<nz;iz++){
					m_threads[imx*nmy + imy][iz] = m[imx*nmy + imy][iz];
				}
			}
		}
	}

	progress = 0.;
#pragma omp parallel for private(iw) shared(m_threads,d_g_wx,d_s_wx)
	for (iw=ifmin;iw<ifmax;iw++){ 
		progress += 1./((float) ifmax - ifmin);
		if (verbose) progress_msg(progress);
		extrap1f(m_threads,d_g_wx,d_s_wx,iw,ifmax,nw,ifmax,ntfft,dw,dkx,dky,nkx,nky,nz,oz,dz,gz,sz,nmx,omx,dmx,nmy,omy,dmy,nthread,vel,po,pd,vref,iref1,iref2,nref,p1,p2,adj,pspi,verbose);
	}
	if (verbose) fprintf(stderr,"\n");
	if (adj){
		// reduction over parallel axis 
		for (imx=0;imx<nmx;imx++) for (imy=0;imy<nmy;imy++) for (ithread=0;ithread<nthread;ithread++) for (iz=0;iz<nz;iz++) m[imx*nmy + imy][iz] += m_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++){
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++){ 
				w = iw*dw;
				d_w[iw] = w*w*d_g_wx[ix][iw];
			}
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) d[ix][it] = d_t[it];
		}
	}

	free1int(n); 
	fftwf_free(a);
	fftwf_free(b);  
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);  
	free1float(d_t);
	free1complex(d_w);
	free2float(m_threads);
	free2complex(d_g_wx);
	free2complex(d_s_wx);
	free2float(vref);
	free2int(iref1);
	free2int(iref2);
	return;
} 

void extrap1f(float **m,complex **d_g_wx, complex **d_s_wx,
		int iw, int ang_iw_max, int nw,int ifmax,int ntfft,float dw,float dkx,float dky,int nkx,int nky,
		int nz, float oz, float dz, float gz, float sz,
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		int nthread,
		float **v,float *po,float **pd,
		float **vref, int **iref1, int **iref2, int nref,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, bool pspi, bool verbose)   
/*< extrapolate 1 frequency >*/
{
	float w,factor,z;
	int iz,ix,imx,imy,ithread;
	complex *d_xg,*d_xs,**smig;
	ithread = omp_get_thread_num(); 
	d_xg = alloc1complex(nmx*nmy);
	d_xs = alloc1complex(nmx*nmy);
	for (ix=0;ix<nmx*nmy;ix++) d_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) d_xs[ix] = 0.;
	if (iw==0) factor = 1.;
	else factor = 2.;
	w = iw*dw;
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++){ 
			d_xs[ix] = d_s_wx[ix][iw]/sqrtf((float) ntfft);
			d_xg[ix] = d_g_wx[ix][iw]/sqrtf((float) ntfft);
		}
		for (iz=0;iz<nz;iz++){ // extrapolate source and receiver wavefields
			z = oz + dz*iz;
			if (z >= sz){
				if (pspi) pspiop(d_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,true,true,verbose);
				else ssop(d_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v,po,pd,p1,p2,true,true,verbose);
			} 
			if (z >= gz){
				if (pspi) pspiop(d_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,true,false,verbose);
				else ssop(d_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,v,po,pd,p1,p2,true,false,verbose);
				for (imx=0;imx<nmx;imx++){ 
					for (imy=0;imy<nmy;imy++){
						m[imx*nmy*nthread + imy*nthread + ithread][iz] += factor*crealf(d_xg[imx*nmy + imy]*conjf(d_xs[imx*nmy + imy]));
					}
				}
			}
		}
	}

	else{
		smig = alloc2complex(nz,nmx*nmy);
		for (ix=0;ix<nmx*nmy;ix++) d_xs[ix] = d_s_wx[ix][iw]/sqrtf((float) ntfft);
		for (iz=0;iz<nz;iz++){ // extrapolate source wavefield 
			z = oz + dz*iz;
			if (z >= sz){
				if (pspi) pspiop(d_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,true,true,verbose);
				else ssop(d_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v,po,pd,p1,p2,true,true,verbose);
				for (ix=0;ix<nmx*nmy;ix++) smig[ix][iz] = d_xs[ix];
			}
			else{
				for (ix=0;ix<nmx*nmy;ix++) smig[ix][iz] = 0.;
			}
		}
		for (ix=0;ix<nmx*nmy;ix++) d_xg[ix] = 0.;
		for (iz=nz-1;iz>=0;iz--){ // extrapolate receiver wavefield 
			z = oz + dz*iz;
			if (z >= gz){
				for (ix=0;ix<nmx*nmy;ix++){ 
					d_xg[ix] += m[ix][iz]*smig[ix][iz];
				}
				if (pspi) pspiop(d_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,false,false,verbose);
				else ssop(d_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v,po,pd,p1,p2,false,false,verbose);
			}
		}
		for (ix=0;ix<nmx*nmy;ix++){
			d_g_wx[ix][iw] = d_xg[ix]/sqrtf((float) ntfft);
		}
		free2complex(smig);
	}

	free1complex(d_xs);
	free1complex(d_xg);

	return;
}

void ssop(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
		float **v,float *po,float **pd,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
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
		if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
		else                           kx = -((float) dkx*nkx - dkx*ikx);
		for (iky=0;iky<nky;iky++){
			if (iky<= (int) nky/2) ky = (float) dky*iky;
			else                           ky = -((float) dky*nky - dky*iky);
			s = (w*w)*(po[iz]*po[iz]) - (kx*kx) - (ky*ky);
			if (s>=0) L = cexp(I*sqrtf(s)*dz); 
			else L = cexp(-0.2*sqrtf(fabs(s))*fabs(dz));
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

void pspiop(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,
		int nmx,float omx,float dmx,int nmy,float omy,float dmy,
		float dz,int iz,
		float **vel,float *po,float **pd,
		float **vref, int **iref1, int **iref2, int nref,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool src,
		bool verbose)
{

	float kx,ky,s;
	complex L;
	int ik,ikx,iky,imx,imy; 
	complex *d_k;
	fftwf_complex *a,*b;
	int lmx,lmy;
	complex **dref;
	int iref;
	float v,vref1,vref2;
	if (nmx>100) lmx=30;
	else lmx=0;
	if (nmy>100) lmy=30;
	else lmy=0;
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	dref = alloc2complex(nref,nmx*nmy);
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
	//	boundary_condition(d_x,nmx,lmx,nmy,lmy);    
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) a[imx*nky + imy] = d_x[imx*nmy + imy];
				else a[imx*nky + imy] = 0.;
			}
		}
	}
	fftwf_execute_dft(p1,a,a);
	for (iref=0;iref<nref;iref++){
		for (ikx=0;ikx<nkx;ikx++){
			if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
			else                           kx = -((float) dkx*nkx - dkx*ikx);
			for (iky=0;iky<nky;iky++){
				if (iky<= (int) nky/2) ky = (float) dky*iky;
				else                           ky = -((float) dky*nky - dky*iky);
				s = (w*w)*(1/(vref[iref][iz]*vref[iref][iz])) - (kx*kx) - (ky*ky);
				if (s>=0) L = cexpf(I*sqrtf(s)*dz);
				else L = cexpf(-0.2*sqrtf(fabs(s))*fabs(dz));
				d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);
			}
		}
		for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){
					dref[imx*nmy + imy][iref] = ((complex) b[imx*nky + imy])/sqrtf((float) nkx*nky);
				}
			}
		}
	}
	for (imx=0;imx<nmx*nmy;imx++){
		v = vel[imx][iz];
		vref1 = vref[iref1[imx][iz]][iz];
		vref2 = vref[iref2[imx][iz]][iz];
		if (vref2 - vref1 > 10.0){
			d_x[imx] = linear_interp(dref[imx][iref1[imx][iz]],dref[imx][iref2[imx][iz]],vref1,vref2,v);
		}
		else{
			d_x[imx] = dref[imx][iref1[imx][iz]];
		}
	}
	//if (adj){
	//	boundary_condition(d_x,nmx,lmx,nmy,lmy);    
	//}

	free1complex(d_k);
	fftwf_free(a);
	fftwf_free(b);
	free2complex(dref);

	return;
}

float linear_interp(complex y1,complex y2,float x1,float x2,float x)
	/*< linear interpolation between two points. x2-x1 must be nonzero. >*/
{
	//return  y1 + (y2-y1)*(x-x1)/(x2-x1);
	return  y1;
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

void compute_angles(float **angx, float **angy, float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		float sx,float sy,
		int nz, float oz, float dz, float sz,
		float **vel_p, float fmin, float fmax,
		bool verbose)
/*< source side incidence angle computation using the one way wave equation. >*/
{
	int iz,ix,igx,igy,ik,iw,it,nw,nkx,nky,nkz,ntfft; 
	float dw,dkx,dky,dkz;
	int ifmin,ifmax;
	float *d_t;
	complex *d_w,**d_s_wx;
	fftwf_complex *a1,*b1;
	fftwf_complex *a2,*b2;
	int *n1,*n2;
	fftwf_plan p1,p2,p3,p4;
	float *po_p,**pd_p;
	float progress;
	float **angx_sign,**angy_sign;

	ntfft = nt;
	nw = (int) ntfft/2 + 1;
	nkx = nmx;
	nky = nmy;
	nkz = nz;
	dkx = 2*PI/((float) nkx)/dmx;
	dky = 2*PI/((float) nky)/dmy;
	dkz = 2*PI/((float) nkz)/dz;
	dw = 2*PI/((float) ntfft)/dt;
	if(fmax*dt*ntfft+1<nw) ifmax = (int) fmax*dt*ntfft + 1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = (int) fmin*dt*ntfft;
	else ifmin = 0;
	d_s_wx = alloc2complex(nw,nmx*nmy);
	d_t = alloc1float(nt);
	d_w = alloc1complex(nw);
	for (it=0;it<nt;it++)  d_t[it] = 0.;  
	for (iw=0;iw<nw;iw++)  d_w[iw] = 0.; 
	/* decompose slowness into layer average, and layer purturbation */
	po_p = alloc1float(nz); 
	pd_p = alloc2float(nz,nmx*nmy); 
	for (iz=0;iz<nz;iz++){
		po_p[iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) po_p[iz] += vel_p[ix][iz];
		po_p[iz] /= (float) nmx*nmy;
		po_p[iz]  = 1./po_p[iz];
		for (ix=0;ix<nmx*nmy;ix++) pd_p[ix][iz] = 1.0/vel_p[ix][iz] - po_p[iz];
	}

	// set up fftw plans and pass them to the OMP region of the code
	// plan for extrapolation
	a1 = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b1 = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	n1 = alloc1int(2); n1[0] = nkx; n1[1] = nky;
	p1 = fftwf_plan_dft(2, n1, a1, a1, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftwf_plan_dft(2, n1, b1, b1, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (ik=0;ik<nkx*nky;ik++) a1[ik] = 0.;
	for (ik=0;ik<nkx*nky;ik++) b1[ik] = 0.;
	fftwf_execute_dft(p1,a1,a1);
	fftwf_execute_dft(p2,b1,b1);
	// plan for calculation of spatial derivatives 
	a2 = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky*nkz);
	b2 = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky*nkz);
	n2 = alloc1int(3); n2[0] = nkx; n2[1] = nky; n2[2] = nz;
	p3 = fftwf_plan_dft(3, n2, a2, a2, FFTW_FORWARD, FFTW_ESTIMATE);
	p4 = fftwf_plan_dft(3, n2, b2, b2, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (ik=0;ik<nkx*nky*nkz;ik++) a2[ik] = 0.;
	for (ik=0;ik<nkx*nky*nkz;ik++) b2[ik] = 0.;
	fftwf_execute_dft(p3,a2,a2);
	fftwf_execute_dft(p4,b2,b2);

	/**********************************************************************/
	igx = (int) (sx - omx)/dmx; /*position to inject source in x-dir*/
	igy = (int) (sy - omy)/dmy; /*position to inject source in y-dir*/

	/* source wavefield*/
	for (ix=0;ix<nmx*nmy;ix++) for (iw=0;iw<nw;iw++) d_s_wx[ix][iw] = 0.;
	for (it=0;it<nt;it++) d_t[it] = wav[0][it];

	f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
	for (iw=0;iw<nw;iw++) d_s_wx[igx*nmy + igy][iw] = d_w[iw];   

	angx_sign = alloc2float(nz,nmx*nmy);
	angy_sign = alloc2float(nz,nmx*nmy);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angx_sign[ix][iz] = 0.0;
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angy_sign[ix][iz] = 0.0;

	progress = 0.;
#pragma omp parallel for private(iw) shared(d_s_wx,angx,angy)
	for (iw=ifmin;iw<ifmax;iw++){ 
		progress += 1./((float) ifmax - ifmin);
		if (verbose) progress_msg(progress);
		extrapolate_source(d_s_wx,angx,angy,angx_sign,angy_sign,
				iw,nw,ifmin,ifmax,ntfft,dw,
				dkx,nkx,nmx,omx,dmx,
				dky,nky,nmy,omy,dmy,
				dkz,nkz,nz,oz,dz,
				sz,
				vel_p,po_p,pd_p,
				p1,p2,p3,p4,verbose);                       
	}
	if (verbose) fprintf(stderr,"\n");
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angx[ix][iz] *= signf1(angx_sign[ix][iz])/(float) (ifmax - ifmin + 1);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angy[ix][iz] *= signf1(angy_sign[ix][iz])/(float) (ifmax - ifmin + 1);

	free1int(n1); 
	fftwf_free(a1);
	fftwf_free(b1);
	free1int(n2); 
	fftwf_free(a2);
	fftwf_free(b2);
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);
	fftwf_destroy_plan(p3);
	fftwf_destroy_plan(p4);
	free1float(d_t);
	free1complex(d_w);
	free2complex(d_s_wx);
	free1float(po_p);
	free2float(pd_p);
	free2float(angx_sign);
	free2float(angy_sign);
	return;
} 

void extrapolate_source(complex **d_s_wx, float **angx, float **angy, float **angx_sign, float **angy_sign,
		int iw,int nw,int ifmin,int ifmax,int ntfft,float dw,
		float dkx,int nkx,int nmx,float omx,float dmx,
		float dky,int nky,int nmy,float omy,float dmy,
		float dkz,int nkz,int nz,float oz,float dz,
		float sz,
		float **v_p,float *po_p,float **pd_p,
		fftwf_plan p1,fftwf_plan p2,fftwf_plan p3,fftwf_plan p4,
		bool verbose)              
/*< extrapolate 1 frequency >*/
{
	float w,z;
	int iz,ix,imx,imy;
	complex *d_xs,*u_s,*u_sx,*u_sy,*u_sz;
	float *u_s_signx,*u_s_signy;

	d_xs = alloc1complex(nmx*nmy);
	u_s = alloc1complex(nmx*nmy*nz);
	for (ix=0;ix<nmx*nmy*nz;ix++) u_s[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) d_xs[ix] = 0.;
	w = iw*dw;
	for (ix=0;ix<nmx*nmy;ix++) d_xs[ix] = d_s_wx[ix][iw]/sqrtf((float) ntfft);
	for (iz=0;iz<nz;iz++){ // extrapolate source wavefield 
		z = oz + dz*iz;
		if (z >= sz){
			ssop_source(d_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v_p,po_p,pd_p,p1,p2,verbose);
			for (imx=0;imx<nmx;imx++) for (imy=0;imy<nmy;imy++) u_s[imx*nmy*nz + imy*nz + iz] = d_xs[imx*nmy + imy];
		}
	}

	u_sx = alloc1complex(nmx*nmy*nz);
	u_sy = alloc1complex(nmx*nmy*nz);
	u_sz = alloc1complex(nmx*nmy*nz);
	u_s_signx = alloc1float(nmx*nmy*nz);
	u_s_signy = alloc1float(nmx*nmy*nz);
	for (ix=0;ix<nmx*nmy*nz;ix++) u_sx[ix] = 0.;
	for (ix=0;ix<nmx*nmy*nz;ix++) u_sy[ix] = 0.;
	for (ix=0;ix<nmx*nmy*nz;ix++) u_sz[ix] = 0.;
	for (ix=0;ix<nmx*nmy*nz;ix++) u_s_signx[ix] = 1.;
	for (ix=0;ix<nmx*nmy*nz;ix++) u_s_signy[ix] = 1.;
	calculate_derivatives(u_s,u_sx,u_sy,u_sz,u_s_signx,u_s_signy,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,dkz,nkz,nz,oz,dz,p3,p4);

	for(imx=0;imx<nmx;imx++){ 
		for(imy=0;imy<nmy;imy++){
			for(iz=0;iz<nz;iz++){
				angx_sign[imx*nmy + imy][iz] += u_s_signx[imx*nmy*nz + imy*nz + iz]/ (float) nw; 
				angy_sign[imx*nmy + imy][iz] += u_s_signy[imx*nmy*nz + imy*nz + iz]/ (float) nw; 
			}
		}
	}
	for(imx=0;imx<nmx;imx++){ 
		for(imy=0;imy<nmy;imy++){
			for(iz=0;iz<nz;iz++){
				angx[imx*nmy + imy][iz] += (180/PI)*atanf(cabs(u_sx[imx*nmy*nz + imy*nz + iz])/(cabs(u_sz[imx*nmy*nz + imy*nz + iz]) + 1e-10)); 
				angy[imx*nmy + imy][iz] += (180/PI)*atanf(cabs(u_sy[imx*nmy*nz + imy*nz + iz])/(cabs(u_sz[imx*nmy*nz + imy*nz + iz]) + 1e-10)); 
			}
		}
	}
	free1complex(d_xs);
	free1complex(u_s);
	free1complex(u_sx);
	free1complex(u_sy);
	free1complex(u_sz);
	free1float(u_s_signx);
	free1float(u_s_signy);
	return;
}

void ssop_source(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
		float **v,float *po,float **pd,
		fftwf_plan p1,fftwf_plan p2, 
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


	boundary_condition(d_x,nmx,lmx,nmy,lmy);    
	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){ 
			if (imx < nmx && imy < nmy){
				L = cexpf(I*w*pd[imx*nmy + imy][iz]*dz);  
				a[imx*nky + imy] = d_x[imx*nmy + imy]*L; // SS operator
			}
			else a[imx*nky + imy] = 0.;
		}
	}


	fftwf_execute_dft(p1,a,a); 
	for (ikx=0;ikx<nkx;ikx++){
		if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
		else                           kx = -((float) dkx*nkx - dkx*ikx);
		for (iky=0;iky<nky;iky++){
			if (iky<= (int) nky/2) ky = (float) dky*iky;
			else                           ky = -((float) dky*nky - dky*iky);
			s = (w*w)*(po[iz]*po[iz]) - (kx*kx) - (ky*ky);
			if (s>=0) L = cexpf(I*sqrtf(s)*dz); 
			else L = cexpf(-0.2*sqrtf(fabs(s))*fabs(dz));
			d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);        
		}
	}
	for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
	fftwf_execute_dft(p2,b,b);

	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){ 
			if (imx < nmx && imy < nmy){ 
				d_x[imx*nmy + imy] = ((complex) b[imx*nky + imy])/sqrtf((float) nkx*nky);
			}
		}
	}

	free1complex(d_k);
	fftwf_free(a);
	fftwf_free(b);

	return;
}

void calculate_derivatives(complex *u_s, complex *u_sx, complex *u_sy, complex *u_sz,
		float *u_s_signx, float *u_s_signy,
		float dkx, int nkx, int nmx, float omx, float dmx,
		float dky, int nky, int nmy, float omy, float dmy,
		float dkz, int nkz, int nz, float oz, float dz, 
		fftwf_plan p3,fftwf_plan p4)              
/*< calculate spatial derivatives of the source wavefield >*/
{
	fftwf_complex *a,*b,*u_left,*u_right;
	int imx,imy,iz;
	float kx,ky,kz,kx_left,kx_right,ky_left,ky_right;

	/* set up fftw plans */
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky*nkz);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky*nkz);
	u_left  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky*nkz);
	u_right  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky*nkz);

	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){
			for(iz=0; iz<nkz;iz++){ 
				if (imx < nmx && imy < nmy && iz < nz){
					a[imx*nky*nkz + imy*nkz + iz] = u_s[imx*nmy*nz + imy*nz + iz]; 
				}
				else a[imx*nky*nkz + imy*nkz + iz] = 0.;
			}
		}
	}
	fftwf_execute_dft(p3,a,a);
	// x component
	for(imx=0;imx<nmx;imx++){ 
		if (imx<= (int) nkx/2) kx = (float) dkx*imx;
		else                           kx = -((float) dkx*nkx - dkx*imx);
		for(imy=0;imy<nky;imy++){
			for(iz=0;iz<nkz;iz++){ 
				b[imx*nky*nkz + imy*nkz + iz] = I*kx*a[imx*nky*nkz + imy*nkz + iz];
			}
		}
	}
	fftwf_execute_dft(p4,b,b);
	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){
			for(iz=0; iz<nkz;iz++){ 
				if (imx < nmx && imy < nmy && iz < nz){
					u_sx[imx*nmy*nz + imy*nz + iz] = b[imx*nky*nkz + imy*nkz + iz]/((float) nkx*nky*nkz); 
				}
			}
		}
	}
	// y component
	for(imx=0;imx<nmx;imx++){ 
		for(imy=0;imy<nky;imy++){
			if (imy<= (int) nky/2) ky = (float) dky*imy;
			else                           ky = -((float) dky*nky - dky*imy);
			for(iz=0;iz<nkz;iz++){ 
				b[imx*nky*nkz + imy*nkz + iz] = I*ky*a[imx*nky*nkz + imy*nkz + iz];
			}
		}
	}
	fftwf_execute_dft(p4,b,b);
	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){
			for(iz=0; iz<nkz;iz++){ 
				if (imx < nmx && imy < nmy && iz < nz){
					u_sy[imx*nmy*nz + imy*nz + iz] = b[imx*nky*nkz + imy*nkz + iz]/((float) nkx*nky*nkz); 
				}
			}
		}
	}
	// z component
	for(imx=0;imx<nmx;imx++){ 
		for(imy=0;imy<nky;imy++){
			for(iz=0;iz<nkz;iz++){ 
				if (iz<= (int) nkz/2) kz = (float) dkz*iz;
				else                          kz = -((float) dkz*nkz - dkz*iz);
				b[imx*nky*nkz + imy*nkz + iz] = I*kz*a[imx*nky*nkz + imy*nkz + iz];
			}
		}
	}
	fftwf_execute_dft(p4,b,b);
	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){
			for(iz=0; iz<nkz;iz++){ 
				if (imx < nmx && imy < nmy && iz < nz){
					u_sz[imx*nmy*nz + imy*nz + iz] = b[imx*nky*nkz + imy*nkz + iz]/((float) nkx*nky*nkz); 
				}
			}
		}
	}
	// sign for x direction
	for(imx=0;imx<nmx;imx++){ 
		if (imx<= (int) nkx/2) kx_left = 0.;
		else                           kx_left = -((float) dkx*nkx - dkx*imx);
		if (imx<= (int) nkx/2) kx_right = (float) dkx*imx;
		else                           kx_right = 0.;
		for(imy=0;imy<nky;imy++){
			for(iz=0;iz<nkz;iz++){ 
				u_left[imx*nky*nkz + imy*nkz + iz]  = I*kx_left*a[imx*nky*nkz + imy*nkz + iz];
				u_right[imx*nky*nkz + imy*nkz + iz] = I*kx_right*a[imx*nky*nkz + imy*nkz + iz];
			}
		}
	}
	fftwf_execute_dft(p4,u_left,u_left);
	fftwf_execute_dft(p4,u_right,u_right);
	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){
			for(iz=0; iz<nkz;iz++){ 
				if (imx < nmx && imy < nmy && iz < nz){
					if (cabs(u_left[imx*nky*nkz + imy*nkz + iz]) >= cabs(u_right[imx*nky*nkz + imy*nkz + iz])){
						u_s_signx[imx*nmy*nz + imy*nz + iz] = 1.*cabs(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
					else{
						u_s_signx[imx*nmy*nz + imy*nz + iz] = -1.*cabs(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
				}
			}
		}
	}
	// sign for y direction
	for(imx=0;imx<nmx;imx++){ 
		for(imy=0;imy<nky;imy++){
			if (imy<= (int) nky/2) ky_left = 0.;
			else                           ky_left = -((float) dky*nky - dky*imy);
			if (imy<= (int) nky/2) ky_right = (float) dky*imy;
			else                           ky_right = 0.;
			for(iz=0;iz<nkz;iz++){ 
				u_left[imx*nky*nkz + imy*nkz + iz]  = I*ky_left*a[imx*nky*nkz + imy*nkz + iz];
				u_right[imx*nky*nkz + imy*nkz + iz] = I*ky_right*a[imx*nky*nkz + imy*nkz + iz];
			}
		}
	}
	fftwf_execute_dft(p4,u_left,u_left);
	fftwf_execute_dft(p4,u_right,u_right);
	for(imx=0; imx<nkx;imx++){ 
		for(imy=0; imy<nky;imy++){
			for(iz=0; iz<nkz;iz++){ 
				if (imx < nmx && imy < nmy && iz < nz){
					if (cabs(u_left[imx*nky*nkz + imy*nkz + iz]) >= cabs(u_right[imx*nky*nkz + imy*nkz + iz])){
						u_s_signy[imx*nmy*nz + imy*nz + iz] = 1.*cabs(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
					else{
						u_s_signy[imx*nmy*nz + imy*nz + iz] = -1.*cabs(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
				}
			}
		}
	}

	fftwf_free(a);
	fftwf_free(b);
	fftwf_free(u_left);
	fftwf_free(u_right);

	return;
}

float signf1(float a)
	/*< sign of a float, always has an amplitude of 1 >*/
{
	float b;
	if (a>=0.) b = 1.;
	else b =-1.;
	return b;
}
