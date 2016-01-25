#include "seismic.h"
#include "wem.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation migration for 3D acoustic, isotropic data. \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char d_name[512],m_name[512],vel_name[512],wav_name[512];
	struct SeisHeader *h_d;
	struct SeisHeader *h_m;
	struct SeisHeader *h_vel;
	struct SeisHeader *h_vel_aperture;
	struct SeisHeader *h_wav;
	int nx,ny,nz,nt,ix,iz,it,nxa,nya,imx,imy;
	int ixmin_aperture,ixmax_aperture,iymin_aperture,iymax_aperture;
	float **d,**m,**vel,**vel_aperture,**wav,sx,sy,sz,gz;
	float ox,dx,oy,dy,oz,dz,ot,dt,fmin,fmax,xmin,xmax,ymin,ymax;
	float ohx,dhx,ohy,dhy;
	float damping;
	int nhx,nhy,ntraces;
	int padt,padx;
	bool adj,pade_flag,verbose;
	struct SeisFileHeader fh;
	
	if (!par_read_bool(argc,argv,"adj",&adj)) adj = true;
	if (!par_read_float(argc,argv,"damping",&damping)) damping = 1000.;
	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_bool(argc,argv,"pade_flag",&pade_flag)) pade_flag= false; // flag for Pade Fourier correction
	if (!par_read_string(argc,argv,"d", d_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"m", m_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vel", vel_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"wav", wav_name)) { docs (); exit (1); }
	if (!par_read_float(argc,argv,"sz",&sz)) sz = 0;
	if (!par_read_float(argc,argv,"gz",&gz)) gz = 0;
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
	if (!par_read_int(argc,argv,"padt",&padt)) padt = 2;
	if (!par_read_int(argc,argv,"padx",&padx)) padx = 2;
	// get dimensions from velocity (nz,oz,dz,nx,ox,dx) and wavelet (nt,sx) files
	InitFileHeader(&fh);
	ReadFileHeader(wav_name,&fh);
	nt = fh.n1;
	ReadFileHeader(vel_name,&fh);
	nz = fh.n1; ntraces = fh.n2*fh.n3*fh.n4*fh.n5;
	
	//fprintf(stderr,"nx=%d\n",nx);
	//fprintf(stderr,"nz=%d\n",nz);
	//fprintf(stderr,"nt=%d\n",nt);
	h_wav = allocSeisHeader(1);
	wav = alloc2float(nt,1); 		
	SeisRead(wav_name,wav,h_wav,&fh);
	ot = fh.o1;
	dt = fh.d1;
	h_vel = allocSeisHeader(ntraces);
	vel = alloc2float(nz,ntraces);
	SeisRead(vel_name,vel,h_vel,&fh); 	
	nx = h_vel[ntraces-1].imx - h_vel[0].imx + 1;
	ny = h_vel[ntraces-1].imy - h_vel[0].imy + 1;
	oz = fh.o1;
	dz = fh.d1;
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
	vel_aperture = alloc2float(nz,nxa*nya);
	h_vel_aperture = allocSeisHeader(nxa*nya);

	for (imx=0;imx<nxa;imx++){
		for (imy=0;imy<nya;imy++){
			for (iz=0;iz<nz;iz++) vel_aperture[imx*nya + imy][iz] = vel[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
			h_vel_aperture[imx*nya + imy] = h_vel[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)];
		}
	}
	free2float(vel); 
	if (adj){
		h_d = allocSeisHeader(nxa*nya);
		d = alloc2float(nt,nxa*nya);
		SeisRead(d_name,d,h_d,&fh);
	}        
	else{
		h_m = allocSeisHeader(nxa*nya);
		m = alloc2float(nz,nxa*nya);
		SeisRead(m_name,m,h_m,&fh);
	}
	if (adj){
		m = alloc2float(nz,nxa*nya);
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) m[ix][iz] = 0.0;
		h_m = allocSeisHeader(nxa*nya);
		for (ix=0;ix<nxa*nya;ix++) h_m[ix] = h_vel[ix];	
	}       
	else{  	
		h_d = allocSeisHeader(nxa*nya);
		d = alloc2float(nt,nxa*nya);
		for (it=0;it<nt;it++) for (ix=0;ix<nxa*nya;ix++) d[ix][it] = 0.0;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix] = h_m[ix];
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].n1 = nt;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].d1 = dt;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].o1 = ot;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].gx = h_m[ix].mx;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].gy = h_m[ix].my;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].hx = h_m[ix].mx - sx;
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].hy = h_m[ix].my - sy;  		
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].ihx = (int) truncf((h_d[ix].hx - ohx)/dhx);
		for (ix=0;ix<nxa*nya;ix++) h_d[ix].ihy = (int) truncf((h_d[ix].hy - ohy)/dhy);
	}
	wem(d,
			m,
			wav,
			nt,ot,dt,
			nxa,xmin,dx,
			nya,ymin,dy,
			sx,sy,
			nz,oz,dz,
			gz,sz,
			vel_aperture,
			fmin,fmax,
			padt,padx,
			damping,
			adj,
			pade_flag,
			verbose);
	if (adj){
		for (ix=0;ix<nxa*nya;ix++) h_m[ix].mx = h_d[ix].gx;
		for (ix=0;ix<nxa*nya;ix++) h_m[ix].my = h_d[ix].gy;
		for (ix=0;ix<nxa*nya;ix++) h_m[ix].imx = (int) truncf((h_m[ix].mx - ox)/dx);
		for (ix=0;ix<nxa*nya;ix++) h_m[ix].imy = (int) truncf((h_m[ix].my - oy)/dy);
		InitFileHeader(&fh);
		fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
		fh.n2 = nxa; fh.o2 = xmin; fh.d2 = dx;
		fh.n3 = nya; fh.o3 = ymin; fh.d3 = dy;
		SeisWrite(m_name,m,h_m,&fh);
	}
	else{
		InitFileHeader(&fh);
		fh.n1 = nt; fh.o1 = ot; fh.d1 = dt;
		fh.n2 = nxa; fh.o2 = xmin; fh.d2 = dx;
		fh.n3 = nya; fh.o3 = ymin; fh.d3 = dy;
		SeisWrite(d_name,d,h_d,&fh);
	}
	free2float(d); 
	free2float(m); 
	free2float(vel_aperture); 
	free2float(wav); 
	freeSeisHeader(h_d);
	freeSeisHeader(h_m);
	freeSeisHeader(h_vel);
	freeSeisHeader(h_wav);
	return 0;
}
