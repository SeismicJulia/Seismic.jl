#include "seismic.h"
#include "wem.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation angle of incidence computation for 3D acoustic, isotropic data. \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char angx_name[512],angy_name[512],vel_name[512],wav_name[512],dipx_name[512],dipy_name[512];
	struct SeisHeader *h_angx;
	struct SeisHeader *h_angy;
	struct SeisHeader *h_vel;
	struct SeisHeader *h_dipx,*h_dipy;
	struct SeisHeader *h_wav;
	int nx,ny,nz,nt,ix,iz,imx,imy;
	float **angx,**angy,**vel,**dipx,**dipy,**wav,sx,sy,sz;
	float ox,dx,oy,dy,oz,dz,ot,dt,fmin,fmax;
	int ntraces;
	bool verbose,dip_flag;
	struct SeisFileHeader fh;
	
	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_bool(argc,argv,"dip_flag",&dip_flag)) dip_flag = false; // flag to correct angles relative to reflector normal
	if (!par_read_string(argc,argv,"angx", angx_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"angy", angy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vel", vel_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"dipx", dipx_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"dipy", dipy_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"wav", wav_name)) { docs (); exit (1); }
	if (!par_read_float(argc,argv,"sz",&sz)) sz = 0;
	if (!par_read_float(argc,argv,"sx",&sx)) sx = 0;
	if (!par_read_float(argc,argv,"sy",&sy)) sy = 0;	
	if (!par_read_float(argc,argv,"fmin",&fmin)) fmin = 0;
	if (!par_read_float(argc,argv,"fmax",&fmax)) fmax = 80;
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
	ot = h_wav[0].o1;
	dt = h_wav[0].d1;
	h_vel = allocSeisHeader(ntraces);
	vel = alloc2float(nz,ntraces);
	SeisRead(vel_name,vel,h_vel,&fh); 	
	nx = h_vel[ntraces-1].imx - h_vel[0].imx + 1;
	ny = h_vel[ntraces-1].imy - h_vel[0].imy + 1;
	oz = h_vel[0].o1;
	dz = h_vel[0].d1;
	ox = h_vel[0].mx;
	dx = h_vel[1].mx - h_vel[0].mx;
	oy = h_vel[0].my;
	dy = ntraces > nx ? h_vel[nx].my - h_vel[nx-1].my : dx;
	angx = alloc2float(nz,nx*ny);
	angy = alloc2float(nz,nx*ny);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nx*ny;ix++) angx[ix][iz] = 0.0;
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nx*ny;ix++) angy[ix][iz] = 0.0;
	h_angx = allocSeisHeader(nx*ny);
	h_angy = allocSeisHeader(nx*ny);
	for (ix=0;ix<nx*ny;ix++) h_angx[ix] = h_vel[ix];
	for (ix=0;ix<nx*ny;ix++) h_angx[ix].sx = sx;
	for (ix=0;ix<nx*ny;ix++) h_angy[ix] = h_vel[ix];
	for (ix=0;ix<nx*ny;ix++) h_angy[ix].sx = sx;
	compute_angles(angx,angy,wav,
			nt,ot,dt, 
			nx,ox,dx,
			ny,oy,dy,
			sx,sy,
			nz,oz,dz,sz,
			vel,fmin,fmax,
			verbose);

	if (dip_flag){
		dipx = alloc2float(nz,ntraces);
		dipy = alloc2float(nz,ntraces);
		h_dipx = allocSeisHeader(ntraces);
		h_dipy = allocSeisHeader(ntraces);
		SeisRead(dipx_name,dipx,h_dipx,&fh); 	
		SeisRead(dipy_name,dipy,h_dipy,&fh); 	
		for (imx=0;imx<nx;imx++){
			for (imy=0;imy<ny;imy++){
				for (iz=0;iz<nz;iz++) angx[imx*ny + imy][iz] += dipx[imx*ny + imy][iz];
				for (iz=0;iz<nz;iz++) angy[imx*ny + imy][iz] += dipy[imx*ny + imy][iz];
			}
		}
		free2float(dipx); 
		free2float(dipy);
		freeSeisHeader(h_dipx);
		freeSeisHeader(h_dipy);
	}

	InitFileHeader(&fh);
	fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
	fh.n2 = nx; fh.o2 = ox; fh.d2 = dx;
	fh.n3 = ny; fh.o3 = oy; fh.d3 = dy;
	SeisWrite(angx_name,angx,h_angx,&fh);
	freeSeisHeader(h_angx);
	SeisWrite(angy_name,angy,h_angy,&fh);
	freeSeisHeader(h_angy);
	free2float(angx); 
	free2float(angy); 
	free2float(vel);
	free2float(wav); 
	freeSeisHeader(h_vel);
	freeSeisHeader(h_wav);
	return 0;
}
