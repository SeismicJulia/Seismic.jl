#include "seismic.h"
#include "ewem.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation migration for 3D-3C elastic, isotropic     \n"
		" data.                                                                \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char ux_name[512],uy_name[512],uz_name[512];
	char mpp_name[512],mps1_name[512],mps2_name[512];
	char vp_name[512],vs_name[512],wav_name[512];
	struct SeisHeader *h_ux;
	struct SeisHeader *h_uy;
	struct SeisHeader *h_uz;
	struct SeisHeader *h_mpp;
	struct SeisHeader *h_mps1;
	struct SeisHeader *h_mps2;
	struct SeisHeader *h_vp;
	struct SeisHeader *h_vp_aperture;
	struct SeisHeader *h_vs;
	struct SeisHeader *h_wav;
	int nx,ny,nz,nt,ix,iz,it,nxa,nya,imx,imy;
	int ixmin_aperture,ixmax_aperture,iymin_aperture,iymax_aperture;
	float **ux,**uy,**uz,**mpp,**mps1,**mps2,**vp,**vs,**vp_aperture,**vs_aperture,**wav,sx,sy,sz,gz;
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
	if (!par_read_string(argc,argv,"ux", ux_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"uy", uy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"uz", uz_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"mpp", mpp_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"mps1", mps1_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"mps2", mps2_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vp", vp_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vs", vs_name)) { docs (); exit (1); }
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
	ReadFileHeader(vp_name,&fh);
	nz = fh.n1; ntraces = fh.n2*fh.n3*fh.n4*fh.n5;

	//fprintf(stderr,"nx=%d\n",nx);
	//fprintf(stderr,"nz=%d\n",nz);
	//fprintf(stderr,"nt=%d\n",nt);
	h_wav = allocSeisHeader(1);
	wav = alloc2float(nt,1); 		
	SeisRead(wav_name,wav,h_wav,&fh);
	ot = fh.o1;
	dt = fh.d1;
	h_vp = allocSeisHeader(ntraces);
	vp = alloc2float(nz,ntraces);
	SeisRead(vp_name,vp,h_vp,&fh); 	
	h_vs = allocSeisHeader(ntraces); 	
	vs = alloc2float(nz,ntraces);
	SeisRead(vs_name,vs,h_vs,&fh); 	
	nx = h_vp[ntraces-1].imx - h_vp[0].imx + 1;
	ny = h_vp[ntraces-1].imy - h_vp[0].imy + 1;
	oz = fh.o1;
	dz = fh.d1;
	ox = h_vp[0].mx;
	dx = dhx;
	oy = h_vp[0].my;
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
	h_vp_aperture = allocSeisHeader(nxa*nya);

	for (imx=0;imx<nxa;imx++){
		for (imy=0;imy<nya;imy++){
			for (iz=0;iz<nz;iz++) vp_aperture[imx*nya + imy][iz] = vp[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
			for (iz=0;iz<nz;iz++) vs_aperture[imx*nya + imy][iz] = vs[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)][iz];
			h_vp_aperture[imx*nya + imy] = h_vp[(imx + ixmin_aperture)*ny + (imy + iymin_aperture)];
		}
	}
	//SeisWrite("vel_aperture",vel_aperture,h_vel_aperture,nz,nxa*nya);
	free2float(vp); 
	free2float(vs); 
	if (adj){
		h_ux = allocSeisHeader(nxa*nya);
		ux = alloc2float(nt,nxa*nya);
		SeisRead(ux_name,ux,h_ux,&fh);
		h_uy = allocSeisHeader(nxa*nya);
		uy = alloc2float(nt,nxa*nya);
		SeisRead(uy_name,uy,h_uy,&fh);
		h_uz = allocSeisHeader(nxa*nya);
		uz = alloc2float(nt,nxa*nya);
		SeisRead(uz_name,uz,h_uz,&fh);
	}        
	else{
		h_mpp = allocSeisHeader(nxa*nya);
		mpp = alloc2float(nz,nxa*nya);
		SeisRead(mpp_name,mpp,h_mpp,&fh);
		h_mps1 = allocSeisHeader(nxa*nya);
		mps1 = alloc2float(nz,nxa*nya);
		SeisRead(mps1_name,mps1,h_mps1,&fh);
		h_mps2 = allocSeisHeader(nxa*nya);
		mps2 = alloc2float(nz,nxa*nya);
		SeisRead(mps2_name,mps2,h_mps2,&fh);
	}
	if (adj){
		mpp = alloc2float(nz,nxa*nya);
		mps1 = alloc2float(nz,nxa*nya);
		mps2 = alloc2float(nz,nxa*nya);
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) mpp[ix][iz] = 0.0;		
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) mps1[ix][iz] = 0.0;
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nxa*nya;ix++) mps2[ix][iz] = 0.0;
		h_mpp = allocSeisHeader(nxa*nya);
		h_mps1 = allocSeisHeader(nxa*nya);
		h_mps2 = allocSeisHeader(nxa*nya);
		for (ix=0;ix<nxa*nya;ix++) h_mpp[ix] = h_vp_aperture[ix];
		for (ix=0;ix<nxa*nya;ix++) h_mps1[ix] = h_vp_aperture[ix];
		for (ix=0;ix<nxa*nya;ix++) h_mps2[ix] = h_vp_aperture[ix];
	}       
	else{
		h_ux = allocSeisHeader(nxa*nya);
		ux = alloc2float(nt,nxa*nya);
		for (it=0;it<nt;it++) for (ix=0;ix<nxa*nya;ix++) ux[ix][it] = 0.0;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix] = h_mpp[ix];
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].n1 = nt;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].d1 = dt;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].o1 = ot;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].gx = h_mpp[ix].mx;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].gy = h_mpp[ix].my;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].hx = h_mpp[ix].mx - sx;
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].hy = h_mpp[ix].my - sy;  		
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].ihx = (int) truncf((h_ux[ix].hx - ohx)/dhx);
		for (ix=0;ix<nxa*nya;ix++) h_ux[ix].ihy = (int) truncf((h_ux[ix].hy - ohy)/dhy);
		h_uy = allocSeisHeader(nxa*nya);
		uy = alloc2float(nt,nxa*nya);
		for (it=0;it<nt;it++) for (ix=0;ix<nxa*nya;ix++) uy[ix][it] = 0.0;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix] = h_mpp[ix];
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].n1 = nt;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].d1 = dt;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].o1 = ot;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].gx = h_mpp[ix].mx;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].gy = h_mpp[ix].my;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].hx = h_mpp[ix].mx - sx;
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].hy = h_mpp[ix].my - sy;  		
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].ihx = (int) truncf((h_ux[ix].hx - ohx)/dhx);
		for (ix=0;ix<nxa*nya;ix++) h_uy[ix].ihy = (int) truncf((h_ux[ix].hx - ohy)/dhy);
		h_uz = allocSeisHeader(nxa*nya);
		uz = alloc2float(nt,nxa*nya);
		for (it=0;it<nt;it++) for (ix=0;ix<nxa*nya;ix++) uz[ix][it] = 0.0;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix] = h_mpp[ix];
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].n1 = nt;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].d1 = dt;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].o1 = ot;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].gx = h_mpp[ix].mx;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].gy = h_mpp[ix].my;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].hx = h_mpp[ix].mx - sx;
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].hy = h_mpp[ix].my - sy;  		
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].ihx = (int) truncf((h_uz[ix].hx - ohx)/dhx);
		for (ix=0;ix<nxa*nya;ix++) h_uz[ix].ihy = (int) truncf((h_uz[ix].hy - ohy)/dhy);
	}
	ewem(ux,uy,uz,
		mpp,mps1,mps2, 
		wav,
		nt,ot,dt, 
		nxa,xmin,dx,
		nya,ymin,dy,
		sx,sy,
		nz,oz,dz,gz,sz,
		vp_aperture,vs_aperture, 
		fmin,fmax,
		padt,padx,
		damping,
		adj,pade_flag,verbose);

	if (adj){
		InitFileHeader(&fh);
		fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
		fh.n2 = nxa; fh.o2 = xmin; fh.d2 = dx;
		fh.n3 = nya; fh.o3 = ymin; fh.d3 = dy;
		for (ix=0;ix<nxa*nya;ix++) h_mpp[ix].trid = 1;
		for (ix=0;ix<nxa*nya;ix++) h_mpp[ix].mx = h_ux[ix].gx;
		for (ix=0;ix<nxa*nya;ix++) h_mpp[ix].my = h_ux[ix].gy;
		for (ix=0;ix<nxa*nya;ix++) h_mpp[ix].imx = (int) truncf((h_mpp[ix].mx - ox)/dx);
		for (ix=0;ix<nxa*nya;ix++) h_mpp[ix].imy = (int) truncf((h_mpp[ix].my - oy)/dy);
		SeisWrite(mpp_name,mpp,h_mpp,&fh);
		for (ix=0;ix<nxa*nya;ix++) h_mps1[ix].trid = 2;
		for (ix=0;ix<nxa*nya;ix++) h_mps1[ix].mx = h_ux[ix].gx;
		for (ix=0;ix<nxa*nya;ix++) h_mps1[ix].my = h_ux[ix].gy;
		for (ix=0;ix<nxa*nya;ix++) h_mps1[ix].imx = (int) truncf((h_mps1[ix].mx - ox)/dx);
		for (ix=0;ix<nxa*nya;ix++) h_mps1[ix].imy = (int) truncf((h_mps1[ix].my - oy)/dy);
		SeisWrite(mps1_name,mps1,h_mps1,&fh);
		for (ix=0;ix<nxa*nya;ix++) h_mps2[ix].trid = 2;
		for (ix=0;ix<nxa*nya;ix++) h_mps2[ix].mx = h_ux[ix].gx;
		for (ix=0;ix<nxa*nya;ix++) h_mps2[ix].my = h_ux[ix].gy;
		for (ix=0;ix<nxa*nya;ix++) h_mps2[ix].imx = (int) truncf((h_mps2[ix].mx - ox)/dx);
		for (ix=0;ix<nxa*nya;ix++) h_mps2[ix].imy = (int) truncf((h_mps2[ix].my - oy)/dy);
		SeisWrite(mps2_name,mps2,h_mps2,&fh);
	}
	else{
		InitFileHeader(&fh);
		fh.n1 = nt; fh.o1 = ot; fh.d1 = dt;
		fh.n2 = nxa; fh.o2 = xmin; fh.d2 = dx;
		fh.n3 = nya; fh.o3 = ymin; fh.d3 = dy;
		sprintf(fh.dname,"%s%s",getenv("DATAPATH"); check if datapath is NULL and just use current directory if so.
		SeisWrite(ux_name,ux,h_ux,&fh);
		SeisWrite(uy_name,uy,h_uy,&fh);
		SeisWrite(uz_name,uz,h_uz,&fh);
	}

	free2float(ux); 
	free2float(uy); 
	free2float(uz); 
	free2float(mpp); 
	free2float(mps1); 
	free2float(mps2); 
	free2float(vp_aperture); 
	free2float(vs_aperture); 
	free2float(wav); 
	freeSeisHeader(h_ux);
	freeSeisHeader(h_uy);
	freeSeisHeader(h_uz);
	freeSeisHeader(h_mpp);
	freeSeisHeader(h_mps1);
	freeSeisHeader(h_mps2);
	freeSeisHeader(h_vp);
	freeSeisHeader(h_vs);
	freeSeisHeader(h_wav);

	return 0;
}
