#include "seismic.h"
#include "wem.h"

void wesep2dop(float **ux,float **uy,
		float **uz,float **up,float **us,
		int nt,float ot,float dt,
		int nx,float ox,float dx,
		int ny,float oy,float dy,
		float vp,float vs,
		float fmin,float fmax,
		bool decomp, bool H,
		bool verbose);

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot Helmholtz decomposition for 3D-3C elastic, isotropic     \n"
		" data.                                                                \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char ux_name[512],uy_name[512],uz_name[512],up_name[512],us_name[512];
	struct SeisHeader *h;
	int nx,ny,nt,ix,it,ntraces;
	float **ux,**uy,**uz,**up,**us;
	float ox,dx,oy,dy,ot,dt,fmin,fmax,vp,vs;
	bool decomp,verbose,H;
	struct SeisFileHeader fh;

	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_bool(argc,argv,"H",&H)) H = true; // flag for exact Helmholtz decomp op (y), or adjoint of inverse Helmholtz decomp op (n)
	if (!par_read_bool(argc,argv,"decomp",&decomp)) decomp = true; // flag for decomposition from data components to potentials (y), or recomposition from potentials to data components (n)
	if (!par_read_string(argc,argv,"ux", ux_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"uy", uy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"uz", uz_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"up", up_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"us", us_name)) { docs (); exit (1); }
	if (!par_read_float(argc,argv,"vp",&vp)) vp = 2000;
	if (!par_read_float(argc,argv,"vs",&vs)) vs = 1000;

	if (!par_read_float(argc,argv,"fmin",&fmin)) fmin = 0;
	if (!par_read_float(argc,argv,"fmax",&fmax)) fmax = 80;
	// get dimensions from files
	
	
	InitFileHeader(&fh);
	if (decomp) ReadFileHeader(ux_name,&fh);
	else        ReadFileHeader(up_name,&fh);
	nt = fh.n1;
	ntraces = fh.n2*fh.n3*fh.n4*fh.n5;
	h = allocSeisHeader(ntraces);
	if (decomp){
		ux = alloc2float(nt,ntraces);
		SeisRead(ux_name,ux,h,&fh);
		uy = alloc2float(nt,ntraces);
		SeisRead(uy_name,uy,h,&fh);
		uz = alloc2float(nt,ntraces);
		SeisRead(uz_name,uz,h,&fh);
		up = alloc2float(nt,ntraces);
		us = alloc2float(nt,ntraces); 
		for (ix=0;ix<ntraces;ix++) for (it=0;it<nt;it++) up[ix][it] = 0.; 
		for (ix=0;ix<ntraces;ix++) for (it=0;it<nt;it++) us[ix][it] = 0.; 
	}        
	else{
		h = allocSeisHeader(ntraces);
		up = alloc2float(nt,ntraces);
		SeisRead(up_name,up,h,&fh);
		us = alloc2float(nt,ntraces);
		SeisRead(us_name,us,h,&fh);
		ux = alloc2float(nt,ntraces);
		uy = alloc2float(nt,ntraces); 
		uz = alloc2float(nt,ntraces);
		for (ix=0;ix<ntraces;ix++) for (it=0;it<nt;it++) ux[ix][it] = 0.; 
		for (ix=0;ix<ntraces;ix++) for (it=0;it<nt;it++) uy[ix][it] = 0.; 
		for (ix=0;ix<ntraces;ix++) for (it=0;it<nt;it++) uz[ix][it] = 0.; 
	}
	ot = h[0].o1;
	dt = h[0].d1;	
	nx = 827 ;//(int) (abs(h[ntraces-1].imx - h[0].imx) + 1)/4.;
	ny = 1; //(int) (abs(h[ntraces-1].imy - h[0].imy) + 1)/4.;
	dx = fabs(h[ntraces-1].gx - h[0].gx)/(float) nx;
	dy = fabs(h[ntraces-1].gy - h[0].gy)/(float) ny;
	if (nx*ny != ntraces) { fprintf(stderr,"Error, nx*ny=%d while ntraces=%d ...look in wavesep.c (or check ordering of input traces) \n",nx*ny,ntraces); exit (1); }
	ox = h[0].gx;
	oy = h[0].gy;

	wesep2dop(ux,uy,uz,up,us,nt,ot,dt,nx,ox,dx,ny,oy,dy,vp,vs,fmin,fmax,decomp,H,verbose);

	if (decomp){
		SeisWrite(up_name,up,h,&fh);
		SeisWrite(us_name,us,h,&fh);
	}
	else {
		SeisWrite(ux_name,ux,h,&fh);
		SeisWrite(uy_name,uy,h,&fh);  	
		SeisWrite(uz_name,uz,h,&fh);  	
	}         

	free2float(up); 
	free2float(us); 
	free2float(ux); 
	free2float(uy); 
	free2float(uz); 
	freeSeisHeader(h);

	return 0;
}

void wesep2dop(float **ux,float **uy,float **uz,
		float **up,float **us,
		int nt,float ot,float dt,
		int nx,float ox,float dx,
		int ny,float oy,float dy,
		float vp,float vs,
		float fmin,float fmax,
		bool decomp, bool H,
		bool verbose)
{
	int ix,iy,ik,ikx,iky,iw,it,nw,nkx,nky,padt,padx,ntfft;
	float dw,dkx,dky,w,kx,ky,s1,s2;
	float kzp,kzs,k,norm;
	int ifmin,ifmax;
	float *d_t;
	complex *d_w,*d_x,*d_y,*d_z,*d_p,*d_s;
	complex **dx_g_wx,**dy_g_wx,**dz_g_wx;
	complex **dp_g_wx,**ds_g_wx;
	fftwf_complex *a,*b;
	int *n;
	fftwf_plan p1,p2;
	padt = 2;
	padx = 2;
	ntfft = padt*nt;
	nw=ntfft/2+1;
	nkx = nx > 1 ? padx*nx : nx;
	dkx = nkx > 1 ? 2*PI/nkx/dx : 0;
	nky = ny > 1 ? padx*ny : ny;
	dky = nky > 1 ? 2*PI/nky/dy : 0;
	dw = 2*PI/ntfft/dt;
	if(fmax*dt*ntfft+1<nw) ifmax = truncf(fmax*dt*ntfft)+1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = truncf(fmin*dt*ntfft);
	else ifmin = 0;
	dx_g_wx = alloc2complex(nw,nx*ny);
	dy_g_wx = alloc2complex(nw,nx*ny);
	dz_g_wx = alloc2complex(nw,nx*ny);
	dp_g_wx = alloc2complex(nw,nx*ny);
	ds_g_wx = alloc2complex(nw,nx*ny);
	d_t = alloc1float(nt);
	d_w = alloc1complex(nw);
	d_x = alloc1complex(nkx*nky);
	d_y = alloc1complex(nkx*nky);
	d_z = alloc1complex(nkx*nky);
	d_p = alloc1complex(nkx*nky);
	d_s = alloc1complex(nkx*nky);
	/* set up fftw plans */
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
	if (decomp){
		for (ix=0;ix<nx*ny;ix++){
			// x component 
			for (it=0;it<nt;it++)        d_t[it] = ux[ix][it];
			f_op(d_w,d_t,nw,nt,1);       /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++)     dx_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) dx_g_wx[ix][iw] = d_w[iw]/sqrtf((float) ntfft);
			for (iw=ifmax;iw<nw;iw++)    dx_g_wx[ix][iw] = 0.;
			// y component 
			for (it=0;it<nt;it++)        d_t[it] = uy[ix][it];
			f_op(d_w,d_t,nw,nt,1);       /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++)     dy_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) dy_g_wx[ix][iw] = d_w[iw]/sqrtf((float) ntfft);
			for (iw=ifmax;iw<nw;iw++)    dy_g_wx[ix][iw] = 0.;
			// z component 
			for (it=0;it<nt;it++)        d_t[it] = uz[ix][it];
			f_op(d_w,d_t,nw,nt,1);       /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++)     dz_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) dz_g_wx[ix][iw] = d_w[iw]/sqrtf((float) ntfft);
			for (iw=ifmax;iw<nw;iw++)    dz_g_wx[ix][iw] = 0.;
		}
		for (iw=0;iw<ifmax;iw++){
			w = iw*dw;
			// x component 
			for(ix=0; ix<nkx;ix++){ 
				for(iy=0; iy<nky;iy++){ 
					if (ix < nx && iy < ny) a[ix*nky + iy] = dx_g_wx[ix*ny + iy][iw];
					else a[ix*nky + iy] = 0.;
				}
			}      
			fftwf_execute_dft(p1,a,a);
			for(ik=0; ik<nkx*nky;ik++) d_x[ik] = a[ik]/sqrtf(nkx*nky);
			// y component 
			for(ix=0; ix<nkx;ix++){ 
				for(iy=0; iy<nky;iy++){ 
					if (ix < nx && iy < ny) a[ix*nky + iy] = dy_g_wx[ix*ny + iy][iw];
					else a[ix*nky + iy] = 0.;
				}
			}      
			fftwf_execute_dft(p1,a,a);
			for(ik=0; ik<nkx*nky;ik++) d_y[ik] = a[ik]/sqrtf(nkx*nky);
			// z component 
			for(ix=0; ix<nkx;ix++){ 
				for(iy=0; iy<nky;iy++){ 
					if (ix < nx && iy < ny) a[ix*nky + iy] = dz_g_wx[ix*ny + iy][iw];
					else a[ix*nky + iy] = 0.;
				}
			}      
			fftwf_execute_dft(p1,a,a);
			for(ik=0; ik<nkx*nky;ik++) d_z[ik] = a[ik]/sqrtf(nkx*nky);
			for (ikx=0;ikx<nkx;ikx++){
				for (iky=0;iky<nky;iky++){
					if (ikx<nkx/2.) kx = dkx*ikx;
					else         kx = -(dkx*nkx - dkx*ikx);
					if (iky<nky/2.) ky = dky*iky;
					else         ky = -(dky*nky - dky*iky);
					s1 = w*w/(vp*vp) - kx*kx - ky*ky;
					s2 = w*w/(vs*vs) - kx*kx - ky*ky;
					if (s1>0) kzp = sqrtf(s1);
					else kzp = 0.;
					if (s2>0) kzs = sqrtf(s2);
					else kzs = 0.;
					// normalize the wavenumbers
					norm = sqrtf(kx*kx + ky*ky + kzp*kzs);
					//fprintf(stderr,"norm=%f\n",norm);
					kx = kx/(norm+1.e-7);
					ky = ky/(norm+1.e-7);
					kzp = kzp/(norm+1.e-7);
					kzs = kzs/(norm+1.e-7);
					k = sqrt(kx*kx + ky*ky); 
					//fprintf(stderr,"denom=%f, kzs=%f \n",denom, kzs);
					if (H){
						d_p[ikx*nky + iky] =  I*kx*d_x[ikx*nky + iky]  + I*kzs*d_z[ikx*nky + iky];  
						d_s[ikx*nky + iky] = -I*kzp*d_x[ikx*nky + iky] + I*kx*d_z[ikx*nky + iky]; 
					}  
					else {
						d_p[ikx*nky + iky] = (  kx*d_x[ikx*nky + iky] + kzs*d_z[ikx*nky + iky]  );  
						d_s[ikx*nky + iky] = (  (kx*kzp/(k+1.e-7))*d_x[ikx*nky + iky] - k*d_z[ikx*nky + iky]); 
					}  
				}
			}
			for (ik=0;ik<nkx*nky;ik++)    b[ik] = d_p[ik];
			fftwf_execute_dft(p2,b,b);
			for(ix=0; ix<nkx;ix++){ 
				for(iy=0; iy<nky;iy++){ 
					if (ix < nx && iy < ny) dp_g_wx[ix*ny + iy][iw] = b[ix*nky + iy];
				}
			}      
			for (ik=0;ik<nkx*nky;ik++)    b[ik] = d_s[ik];
			fftwf_execute_dft(p2,b,b);
			for(ix=0; ix<nkx;ix++){ 
				for(iy=0; iy<nky;iy++){ 
					if (ix < nx && iy < ny) ds_g_wx[ix*ny + iy][iw] = b[ix*nky + iy];
				}
			}      
		} 
		for (ix=0;ix<nx*ny;ix++){
			// p component 
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++){ 
				w = iw*dw;
				d_w[iw] = dp_g_wx[ix][iw];
			}
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) up[ix][it] = d_t[it]/sqrtf((float) ntfft);
			// s component 
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++){ 
				w = iw*dw;
				d_w[iw] = ds_g_wx[ix][iw];
			}
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) us[ix][it] = d_t[it]/sqrtf((float) ntfft);
		}
	}
	else{
		fprintf(stderr,"not developed yet. why not write it now? \n");
	}

	fftwf_destroy_plan(p1);fftwf_destroy_plan(p2);
	free1float(d_t);
	free1complex(d_w);
	free2complex(dx_g_wx);
	free2complex(dy_g_wx);
	free2complex(dz_g_wx);
	free2complex(dp_g_wx);
	free2complex(ds_g_wx);
	free1complex(d_x);
	free1complex(d_z);
	free1complex(d_p);
	free1complex(d_s);

	return;
}
