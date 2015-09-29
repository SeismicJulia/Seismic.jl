#include "seismic.h"
#include "ewem.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation angle of incidence computation for 3C-3D elastic, isotropic data. \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

void compute_eangles(float **angpx, float **angpy, float **angsx, float **angsy, float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		float sx,float sy,
		int nz, float oz, float dz, float sz,
		float **vel_p, float **vel_s, float fmin, float fmax,
		bool pade_flag, bool verbose);
void extrapolate_source(complex **d_s_wx, float **angx, float **angy, float **angx_sign, float **angy_sign,
		int iw,int nw,int ifmin,int ifmax,int ntfft,float dw,
		float dkx,int nkx,int nmx,float omx,float dmx,
		float dky,int nky,int nmy,float omy,float dmy,
		float dkz,int nkz,int nz,float oz,float dz,
		float sz,
		float **v_p,float *po_p,float **pd_p,
		fftwf_plan p1,fftwf_plan p2,fftwf_plan p3,fftwf_plan p4,
		bool pade_flag, bool verbose);              
void ssop_source(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
		float **v,float *po,float **pd,
		fftwf_plan p1,fftwf_plan p2, 
		bool pade_flag,
		bool verbose);
void calculate_derivatives(complex *u_s, complex *u_sx, complex *u_sy, complex *u_sz,
		float *u_s_signx, float *u_s_signy,
		float dkx, int nkx, int nmx, float omx, float dmx,
		float dky, int nky, int nmy, float omy, float dmy,
		float dkz, int nkz, int nz, float oz, float dz, 
		fftwf_plan p3,fftwf_plan p4);
float signf1(float a);

int main (int argc, char *argv[])
{
	char angpx_name[512],angpy_name[512],angsx_name[512],angsy_name[512],velp_name[512],vels_name[512],wav_name[512],dipx_name[512],dipy_name[512];
	struct SeisHeader *h_ang;
	struct SeisHeader *h_vel;
	struct SeisHeader *h_vel_aperture;
	struct SeisHeader *h_dipx,*h_dipy;
	struct SeisHeader *h_wav;
	int nx,ny,nz,nt,ix,iz,nxa,nya,imx,imy;
	int ixmin_aperture,ixmax_aperture,iymin_aperture,iymax_aperture;
	float **angpx,**angpy,**angsx,**angsy,**vp,**vp_aperture,**vs,**vs_aperture,**dipx,**dipy,**wav,sx,sy,sz;
	float ox,dx,oy,dy,oz,dz,ot,dt,fmin,fmax,xmin,xmax,ymin,ymax;
	float ohx,dhx,ohy,dhy;
	int nhx,nhy,ntraces;
	bool pade_flag,verbose,dip_flag;
	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_bool(argc,argv,"dip_flag",&dip_flag)) dip_flag = false; // flag to correct angles relative to reflector normal
	if (!par_read_bool(argc,argv,"pade_flag",&pade_flag)) pade_flag = false; // flag for Pade Fourier correction
	if (!par_read_string(argc,argv,"angpx", angpx_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"angpy", angpy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"angsx", angsx_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"angsy", angsy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vp", velp_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vs", vels_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"dipx", dipx_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"dipy", dipy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"wav", wav_name)) { docs (); exit (1); }
	if (!par_read_float(argc,argv,"sz",&sz)) sz = 0;
	if (!par_read_float(argc,argv,"sx",&sx)) sx = 0;
	if (!par_read_float(argc,argv,"sy",&sy)) sy = 0;	
	if (!par_read_float(argc,argv,"ohx",&ohx)) ohx = -1000;
	if (!par_read_float(argc,argv,"dhx",&dhx)) dhx = 10;
	if (!par_read_int(argc,argv,"nhx",&nhx)) nhx = 201;
	if (!par_read_float(argc,argv,"ohy",&ohy)) ohy = -1000;
	if (!par_read_float(argc,argv,"dhy",&dhy)) dhy = 10;
	if (!par_read_int(argc,argv,"nhy",&nhy)) nhy = 201;
	if (!par_read_float(argc,argv,"fmin",&fmin)) fmin = 0;
	if (!par_read_float(argc,argv,"fmax",&fmax)) fmax = 80;
	// get dimensions from velocity (nz,oz,dz,nx,ox,dx) and wavelet (nt,sx) files
	SeisDim(wav_name,&nt,&ntraces);
	SeisDim(velp_name,&nz,&ntraces);
	//fprintf(stderr,"nx=%d\n",nx);
	//fprintf(stderr,"nz=%d\n",nz);
	//fprintf(stderr,"nt=%d\n",nt);
	h_wav = allocSeisHeader(1);
	wav = alloc2float(nt,1); 		
	SeisRead(wav_name,wav,h_wav,nt,1);
	ot = h_wav[0].o1;
	dt = h_wav[0].d1;
	h_vel = allocSeisHeader(ntraces);
	vp = alloc2float(nz,ntraces);
	vs = alloc2float(nz,ntraces);
	SeisRead(velp_name,vp,h_vel,nz,ntraces); 	
	SeisRead(vels_name,vs,h_vel,nz,ntraces);
	nx = h_vel[ntraces-1].imx - h_vel[0].imx + 1;
	ny = h_vel[ntraces-1].imy - h_vel[0].imy + 1;
	oz = h_vel[0].o1;
	dz = h_vel[0].d1;
	ox = h_vel[0].mx;
	dx = dhx;
	oy = h_vel[0].my;
	dy = dhy;
	// calculate number of x points within the shot aperture
	xmin = ox > sx + ohx ? ox : sx + ohx;
	xmax = ox + (nx-1)*dx < sx + ohx + (nhx-1)*dhx ? ox + (nx-1)*dx : sx + ohx + (nhx-1)*dhx;
	ixmin_aperture = (int) truncf((xmin - ox)/dx);
	ixmax_aperture = (int) truncf((xmax - ox)/dx);
	nxa = ixmax_aperture - ixmin_aperture + 1;
	// calculate number of x points within the shot aperture
	ymin = oy > sy + ohy ? oy : sy + ohy;
	ymax = oy + (ny-1)*dy < sy + ohy + (nhy-1)*dhy ? oy + (ny-1)*dy : sy + ohy + (nhy-1)*dhy;
	iymin_aperture = (int) truncf((ymin - oy)/dy);
	iymax_aperture = (int) truncf((ymax - oy)/dy);
	nya = iymax_aperture - iymin_aperture + 1;
	// get part of velocity, data, model that fall within the shot aperture
	vp_aperture = alloc2float(nz,nxa*nya);
	vs_aperture = alloc2float(nz,nxa*nya);
	h_vel_aperture = allocSeisHeader(nxa*nya);

	for (imx=0;imx<nxa;imx++){
		for (imy=0;imy<nya;imy++){
			for (iz=0;iz<nz;iz++) vp_aperture[imx*nya + imy][iz] = vp[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
			for (iz=0;iz<nz;iz++) vs_aperture[imx*nya + imy][iz] = vs[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
			h_vel_aperture[imx*nya + imy] = h_vel[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)];
		}
	}
	free2float(vp); 
	free2float(vs); 

	angpx = alloc2float(nz,nxa*nya);
	angpy = alloc2float(nz,nxa*nya);
	angsx = alloc2float(nz,nxa*nya);
	angsy = alloc2float(nz,nxa*nya);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) angpx[ix][iz] = 0.0;
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) angpy[ix][iz] = 0.0;
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) angsx[ix][iz] = 0.0;
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) angsy[ix][iz] = 0.0;
	h_ang = allocSeisHeader(nxa*nya);
	for (ix=0;ix<nxa*nya;ix++) h_ang[ix] = h_vel_aperture[ix];
	for (ix=0;ix<nxa*nya;ix++) h_ang[ix].sx = sx;
	compute_eangles(angpx,angpy,angsx,angsy,
			wav,
			nt,ot,dt, 
			nxa,xmin,dx,
			nya,ymin,dy,
			sx,sy,
			nz,oz,dz,sz,
			vp_aperture,vs_aperture,fmin,fmax,
			pade_flag,verbose);

	if (dip_flag){
		dipx = alloc2float(nz,ntraces);
		dipy = alloc2float(nz,ntraces);
		h_dipx = allocSeisHeader(ntraces);
		h_dipy = allocSeisHeader(ntraces);
		SeisRead(dipx_name,dipx,h_dipx,nz,ntraces); 	
		SeisRead(dipy_name,dipy,h_dipy,nz,ntraces); 	
		for (imx=0;imx<nxa;imx++){
			for (imy=0;imy<nya;imy++){
				for (iz=0;iz<nz;iz++) angpx[imx*nya + imy][iz] += dipx[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
				for (iz=0;iz<nz;iz++) angpy[imx*nya + imy][iz] += dipy[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
				for (iz=0;iz<nz;iz++) angsx[imx*nya + imy][iz] += dipx[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
				for (iz=0;iz<nz;iz++) angsy[imx*nya + imy][iz] += dipy[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
			}
		}
		free2float(dipx); 
		free2float(dipy);
		freeSeisHeader(h_dipx);
		freeSeisHeader(h_dipy);
	}

	SeisWrite(angpx_name,angpx,h_ang,nz,nxa*nya);
	SeisWrite(angpy_name,angpy,h_ang,nz,nxa*nya);
	SeisWrite(angsx_name,angsx,h_ang,nz,nxa*nya);
	SeisWrite(angsy_name,angsy,h_ang,nz,nxa*nya);
	free2float(angpx); 
	free2float(angpy); 
	free2float(angsx); 
	free2float(angsy); 
	freeSeisHeader(h_ang);
	free2float(vp_aperture);
	free2float(vs_aperture);
	free2float(wav); 
	freeSeisHeader(h_vel);
	freeSeisHeader(h_wav);
	return 0;
}

void compute_eangles(float **angpx, float **angpy, float **angsx, float **angsy, float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		float sx,float sy,
		int nz, float oz, float dz, float sz,
		float **vel_p, float **vel_s, float fmin, float fmax,
		bool pade_flag, bool verbose)
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
	float **angpx_sign,**angpy_sign;

	ntfft = (int) 2*truncf(1*((float) nt)/2);
	nw = (int) truncf(ntfft/2)+1;
	nkx = nmx;
	nky = nmy;
	nkz = nz;
	dkx = 2*PI/((float) nkx)/dmx;
	dky = 2*PI/((float) nky)/dmy;
	dkz = 2*PI/((float) nkz)/dz;

	//fprintf(stderr,"dkx=%f dky=%f dkz=%f \n",dkx,dky,dkz);
	//fprintf(stderr,"nkx=%d nky=%d nkz=%d \n",nkx,nky,nkz);

	dw = 2*PI/((float) ntfft)/dt;
	if(fmax*dt*ntfft+1<nw) ifmax = trunc(fmax*dt*ntfft)+1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = trunc(fmin*dt*ntfft);
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
	igx = (int) truncf((sx - omx)/dmx); /*position to inject source in x-dir*/
	igy = (int) truncf((sy - omy)/dmy); /*position to inject source in y-dir*/

	/* source wavefield*/
	for (ix=0;ix<nmx*nmy;ix++) for (iw=0;iw<nw;iw++) d_s_wx[ix][iw] = 0.;
	for (it=0;it<nt;it++) d_t[it] = wav[0][it];

	f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
	for (iw=0;iw<nw;iw++) d_s_wx[igx*nmy + igy][iw] = d_w[iw];   

	angpx_sign = alloc2float(nz,nmx*nmy);
	angpy_sign = alloc2float(nz,nmx*nmy);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angpx_sign[ix][iz] = 0.0;
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angpy_sign[ix][iz] = 0.0;

	progress = 0.;
#pragma omp parallel for private(iw) shared(d_s_wx,angpx,angpy)
	for (iw=ifmin;iw<ifmax;iw++){ 
		progress += 1./((float) ifmax - ifmin);
		if (verbose) progress_msg(progress);
		extrapolate_source(d_s_wx,angpx,angpy,angpx_sign,angpy_sign,
				iw,nw,ifmin,ifmax,ntfft,dw,
				dkx,nkx,nmx,omx,dmx,
				dky,nky,nmy,omy,dmy,
				dkz,nkz,nz,oz,dz,
				sz,
				vel_p,po_p,pd_p,
				p1,p2,p3,p4,pade_flag,verbose);                       
	}
	if (verbose) fprintf(stderr,"\n");
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angpx[ix][iz] *= signf1(angpx_sign[ix][iz])/(float) (ifmax - ifmin + 1);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angpy[ix][iz] *= signf1(angpy_sign[ix][iz])/(float) (ifmax - ifmin + 1);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angsx[ix][iz] = (180./PI)*asinf(sinf((PI/180.)*angpx[ix][iz])*vel_s[ix][iz]/vel_p[ix][iz]);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nmx*nmy;ix++) angsy[ix][iz] = (180./PI)*asinf(sinf((PI/180.)*angpy[ix][iz])*vel_s[ix][iz]/vel_p[ix][iz]);

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
	free2float(angpx_sign);
	free2float(angpy_sign);
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
		bool pade_flag, bool verbose)              
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
			ssop_source(d_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,-dz,iz,v_p,po_p,pd_p,p1,p2,pade_flag,verbose);
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
				angx[imx*nmy + imy][iz] += (180./PI)*atanf(cabsf(u_sx[imx*nmy*nz + imy*nz + iz])/(cabsf(u_sz[imx*nmy*nz + imy*nz + iz]) + 1.e-10)); 
				angy[imx*nmy + imy][iz] += (180./PI)*atanf(cabsf(u_sy[imx*nmy*nz + imy*nz + iz])/(cabsf(u_sz[imx*nmy*nz + imy*nz + iz]) + 1.e-10)); 
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
		bool pade_flag,
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
	if (pade_flag) pade(d_x,nmx,omx,dmx,nmy,omy,dmy,dz,w,iz,v,po,pd,false,true,verbose);
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
		if (ikx<= (int) truncf(nkx/2)) kx = (float) dkx*ikx;
		else                           kx = -((float) dkx*nkx - dkx*ikx);
		for (iky=0;iky<nky;iky++){
			if (iky<= (int) truncf(nky/2)) ky = (float) dky*iky;
			else                           ky = -((float) dky*nky - dky*iky);
			s = (w*w)*(po[iz]*po[iz]) - (kx*kx) - (ky*ky);
			if (s>=0) L = cexpf(I*sqrtf(s)*dz); 
			else L = cexpf(-0.2*sqrtf(fabsf(s))*fabsf(dz));
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
		if (imx<= (int) truncf(nkx/2)) kx = (float) dkx*imx;
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
			if (imy<= (int) truncf(nky/2)) ky = (float) dky*imy;
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
				if (iz<= (int) truncf(nkz/2)) kz = (float) dkz*iz;
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
		if (imx<= (int) truncf(nkx/2)) kx_left = 0.;
		else                           kx_left = -((float) dkx*nkx - dkx*imx);
		if (imx<= (int) truncf(nkx/2)) kx_right = (float) dkx*imx;
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
					if (cabsf(u_left[imx*nky*nkz + imy*nkz + iz]) >= cabsf(u_right[imx*nky*nkz + imy*nkz + iz])){
						u_s_signx[imx*nmy*nz + imy*nz + iz] = 1.*cabsf(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
					else{
						u_s_signx[imx*nmy*nz + imy*nz + iz] = -1.*cabsf(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
				}
			}
		}
	}
	// sign for y direction
	for(imx=0;imx<nmx;imx++){ 
		for(imy=0;imy<nky;imy++){
			if (imy<= (int) truncf(nky/2)) ky_left = 0.;
			else                           ky_left = -((float) dky*nky - dky*imy);
			if (imy<= (int) truncf(nky/2)) ky_right = (float) dky*imy;
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
					if (cabsf(u_left[imx*nky*nkz + imy*nkz + iz]) >= cabsf(u_right[imx*nky*nkz + imy*nkz + iz])){
						u_s_signy[imx*nmy*nz + imy*nz + iz] = 1.*cabsf(u_s[imx*nmy*nz + imy*nz + iz]);
					} 
					else{
						u_s_signy[imx*nmy*nz + imy*nz + iz] = -1.*cabsf(u_s[imx*nmy*nz + imy*nz + iz]);
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
