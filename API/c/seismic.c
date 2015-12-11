#include "seismic.h"
int par_read_name(int argc, char *argv[], int arg_count, char *arg, int *n) {
  char *buf;
  buf = argv[arg_count];
  for (*n = 0; (buf[*n] != '=') && (buf[*n] != '\0'); (*n)++);
  strncpy (arg, buf, (*n)+1);
  arg[*n] = '\0';
  return 0;
}

int par_read_string(int argc, char *argv[], char *arg, char *value) {
  char name[256], *buf;
  int i, n;
  for (i = 1; i < argc; i++) {
    par_read_name(argc,argv, i, name, &n);
    if (strcmp (arg, name) == 0) {
      buf = argv[i];
      strcpy(value, buf+n+1);
      buf = 0;
      return 1;
    }
  }
  return 0;
}

int par_read_int(int argc, char *argv[], char *arg, int *value) {
  char name[512], *buf;
  int i, n;
  for (i = 1; i < argc; i++) {
    par_read_name(argc,argv, i, name, &n);

    if (strcmp (arg, name) == 0) {
      buf = argv[i];
      *value = atoi (buf+n+1);
      return 1;
    }
  }
  return 0;
}

int par_read_bool(int argc, char *argv[], char *arg, bool *value) {
  char name[512], *buf;
  int i, n;

  for (i = 1; i < argc; i++) {
    par_read_name(argc,argv, i, name, &n);

    if (strcmp (arg, name) == 0) {
      buf = argv[i]+n+1;
      *value = (bool) (('y' == (int) buf[0]));
      return 1;
    }
  }
  return 0;
}


int par_read_float(int argc, char *argv[], char *arg, float *value) {
  char name[512], *buf;
  int i, n;
  for (i = 1; i < argc; i++) {
    par_read_name(argc,argv, i, name, &n);

    if (strcmp (arg, name) == 0) {
      buf = argv[i];
      *value = (float) atof(buf+n+1);
      return 1;
    }
  }
  return 0;
}

/* allocation and deallocation routines from cwp */

/* allocate a 1-d array */
void *alloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}

/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
	size_t i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}

/* allocate a 3-d array */
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
	size_t i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* allocate a 1-d array of ints */
int *alloc1int(size_t n1)
{
	return (int*)alloc1(n1,sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(size_t n1, size_t n2)
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(size_t n1, size_t n2, size_t n3)
{
	return (int***)alloc3(n1,n2,n3,sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(size_t n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of complexs */
complex *alloc1complex(size_t n1)
{
	return (complex*)alloc1(n1,sizeof(complex));
}

/* free a 1-d array of complexs */
void free1complex(complex *p)
{
	free1(p);
}

/* allocate a 2-d array of complexs */
complex **alloc2complex(size_t n1, size_t n2)
{
	return (complex**)alloc2(n1,n2,sizeof(complex));
}

/* free a 2-d array of complexs */
void free2complex(complex **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of complexs */
complex ***alloc3complex(size_t n1, size_t n2, size_t n3)
{
	return (complex***)alloc3(n1,n2,n3,sizeof(complex));
}

/* free a 3-d array of complexs */
void free3complex(complex ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of SeisHeaders */
struct SeisHeader *allocSeisHeader(size_t n1)
{
	return (struct SeisHeader*)alloc1(n1,sizeof(struct SeisHeader));
}

/* free a 1-d array of SeisHeaders */
void freeSeisHeader(struct SeisHeader *p)
{
	free1(p);
}

/* read a seis file */
void SeisRead(char *in,float **d,struct SeisHeader *h, int nt, int nx)
{
	FILE *fp_d,*fp_h;
	char in_d[512],in_h[512];
	int ix;
	sprintf(in_d, "%s.seisd",in);
	sprintf(in_h, "%s.seish",in);  	
  
	fp_h = fopen(in_h,"rb");
	if (!fp_h){
		printf("Unable to open %s\n",in_h);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fread(&h[ix],sizeof(struct SeisHeader),1,fp_h);
	}
	fclose(fp_h);

	fp_d=fopen(in_d,"rb");
	if (!fp_d){
		printf("Unable to open %s\n",in_d);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fread(d[ix],sizeof(float),nt,fp_d);
	}	
	fclose(fp_d);
	
	return;
}

/* write a seis file */
void SeisWrite(char *out,float **d,struct SeisHeader *h, int nt, int nx)
{
	FILE *fp_d,*fp_h;
	char out_d[512],out_h[512];
	int ix;
	sprintf(out_d, "%s.seisd",out);
	sprintf(out_h, "%s.seish",out);  	
  
	fp_h = fopen(out_h,"wb");
	if (!fp_h){
		printf("Unable to open %s\n",out_h);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fwrite(&h[ix],sizeof(struct SeisHeader),1,fp_h);
	}
	fclose(fp_h);

	fp_d=fopen(out_d,"wb");
	if (!fp_d){
		printf("Unable to open %s\n",out_d);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fwrite(d[ix],sizeof(float),nt,fp_d);
	}
	fclose(fp_d);
	
	return;
}

/* return the dimensions (nt and nx) of a seis file */
void SeisDim(char *in,int *nt, int *nx)
{
	FILE *fp_h;
	char in_h[512];
	struct SeisHeader h;
	struct stat st;
	
	sprintf(in_h, "%s.seish",in);  	
	fp_h = fopen(in_h,"rb");
	if (!fp_h){
		printf("Unable to open %s\n",in_h);
		return;
	}
	fread(&h,sizeof(struct SeisHeader),1,fp_h);
	fclose(fp_h);
	stat(in_h, &st);
    *nt = h.n1;
    *nx = (int) truncf((float) st.st_size/124);
	return;
}

void InitSeisHeader(struct SeisHeader *h, int nx)
{
	int ix;
	for (ix=0;ix<nx;ix++){
    	h[ix].tracenum = 0;
    	h[ix].o1    = 0;
    	h[ix].n1    = 0;
    	h[ix].d1    = 0;
    	h[ix].sx    = 0;
    	h[ix].sy    = 0;
    	h[ix].gx    = 0;
    	h[ix].gy    = 0;
    	h[ix].mx    = 0;
    	h[ix].my    = 0;
    	h[ix].hx    = 0;
    	h[ix].hy    = 0;
    	h[ix].h     = 0;
    	h[ix].az    = 0;
    	h[ix].ang   = 0;
    	h[ix].isx   = 0;
    	h[ix].isy   = 0;
    	h[ix].igx   = 0;
    	h[ix].igy   = 0;
    	h[ix].imx   = 0;
    	h[ix].imy   = 0;
    	h[ix].ihx   = 0;
    	h[ix].ihy   = 0;
    	h[ix].ih    = 0;
    	h[ix].iaz   = 0;
    	h[ix].iang  = 0;
    	h[ix].selev = 0;
    	h[ix].gelev = 0;
    	h[ix].sstat = 0;
    	h[ix].gstat = 0;
    	h[ix].trid  = 0;
    }
};



