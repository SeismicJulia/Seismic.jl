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
void SeisRead(char *in,float **d,struct SeisHeader *h, struct SeisFileHeader *fh)
{
	int ix,nx;
	FILE *fp_d,*fp_h;
	InitFileHeader(fh);
	ReadFileHeader(in,fh);
	nx = fh->n2*fh->n3*fh->n4*fh->n5;
	fp_h = fopen(fh->hname,"rb");
	if (!fp_h){
		printf("SeisRead: Unable to open %s\n",fh->hname);
		return;
	}
	for (ix=0; ix<nx; ix++){
		fread(&h[ix],sizeof(struct SeisHeader),1,fp_h);
	}
	fclose(fp_h);
	fp_d=fopen(fh->dname,"rb");
	if (!fp_d){
		printf("SeisRead: Unable to open %s\n",fh->dname);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fread(d[ix],sizeof(float),fh->n1,fp_d);
	}	
	fclose(fp_d);
	
	return;
}

/* write a seis file */
void SeisWrite(char *out,float **d,struct SeisHeader *h, struct SeisFileHeader *fh)
{
	int ix,nx;
	FILE *fp_d,*fp_h;
	const char* datapath;
	// get DATAPATH environmental variable if it exists, and prepend it to dname and hname. 
	datapath = getenv("DATAPATH");
	if (datapath != NULL){
		strcpy (fh->dname,datapath);
		strcat (fh->dname,out);
		strcpy (fh->hname,datapath);
		strcat (fh->hname,out);
	}
	else{
		strcpy (fh->dname,out);
		strcpy (fh->hname,out);
	}	
	strcat (fh->dname,"@data@");
	strcat (fh->hname,"@headers@");
	
	WriteFileHeader(out,fh);
	nx = fh->n2*fh->n3*fh->n4*fh->n5;
  
	fp_h = fopen(fh->hname,"wb");
	if (!fp_h){
		printf("SeisWrite: Unable to open %s\n",fh->hname);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fwrite(&h[ix],sizeof(struct SeisHeader),1,fp_h);
	}
	fclose(fp_h);

	fp_d=fopen(fh->dname,"wb");
	if (!fp_d){
		printf("SeisWrite: Unable to open %s\n",fh->dname);
		return;
	}
	for (ix=0; ix<nx; ix++){	
		fwrite(d[ix],sizeof(float),fh->n1,fp_d);
	}
	fclose(fp_d);
	
	return;
}

void InitFileHeader(struct SeisFileHeader *fh)
{
	fh->n1 = 1;
	fh->n2 = 1;
	fh->n3 = 1;
	fh->n4 = 1;
	fh->n5 = 1;
	fh->o1 = 0;
	fh->o2 = 1;
	fh->o3 = 1;
	fh->o4 = 1;
	fh->o5 = 1;
	fh->d1 = 0;
	fh->d2 = 1;
	fh->d3 = 1;
	fh->d4 = 1;
	fh->d5 = 1;
	strcpy(fh->label1,"");
	strcpy(fh->label2,"");
	strcpy(fh->label3,"");
	strcpy(fh->label4,"");
	strcpy(fh->label5,"");
	strcpy(fh->unit1,"");
	strcpy(fh->unit2,"");
	strcpy(fh->unit3,"");
	strcpy(fh->unit4,"");
	strcpy(fh->unit5,"");
	strcpy(fh->title,"");
	strcpy(fh->data_format,"native_float");
	fh->esize = 4;
	strcpy(fh->dname,"NULL");
	strcpy(fh->hname,"NULL");
}

void ReadFileHeader(char *filename,struct SeisFileHeader *fh)
{
	FILE *fp;
	char fline[100];
	char *newline = NULL;
	char tmp[100];
	fp = fopen(filename, "r");
	while (fgets(fline, 100, fp) != NULL) {
		if (newline == strchr(fline, '\n')) *newline = '\0';
		if (strstr(fline, "in=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->dname);
		if (strstr(fline, "headers=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->hname);
		if (strstr(fline, "title=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->title);
		if (strstr(fline, "label1=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->label1);
		if (strstr(fline, "unit1=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->unit1);
		if (strstr(fline, "n1=") != NULL) sscanf(fline,"%[^= ]=%d", tmp, &fh->n1);
		if (strstr(fline, "o1=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->o1);
		if (strstr(fline, "d1=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->d1);
		if (strstr(fline, "label2=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->label2);
		if (strstr(fline, "unit2=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->unit2);
		if (strstr(fline, "n2=") != NULL) sscanf(fline,"%[^= ]=%d", tmp, &fh->n2);
		if (strstr(fline, "o2=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->o2);
		if (strstr(fline, "d2=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->d2);
		if (strstr(fline, "label3=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->label3);
		if (strstr(fline, "unit3=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->unit3);
		if (strstr(fline, "n3=") != NULL) sscanf(fline,"%[^= ]=%d", tmp, &fh->n3);
		if (strstr(fline, "o3=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->o3);
		if (strstr(fline, "d3=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->d3);
		if (strstr(fline, "label4=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->label4);
		if (strstr(fline, "unit4=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->unit4);
		if (strstr(fline, "n4=") != NULL) sscanf(fline,"%[^= ]=%d", tmp, &fh->n4);
		if (strstr(fline, "o4=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->o4);
		if (strstr(fline, "d4=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->d4);
		if (strstr(fline, "label5=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->label5);
		if (strstr(fline, "unit5=") != NULL) sscanf(fline,"%[^=\" ]=\"%99[^\"]", tmp, fh->unit5);
		if (strstr(fline, "n5=") != NULL) sscanf(fline,"%[^= ]=%d", tmp, &fh->n5);
		if (strstr(fline, "o5=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->o5);
		if (strstr(fline, "d5=") != NULL) sscanf(fline,"%[^= ]=%f", tmp, &fh->d5);
	}
	fclose(fp);
	return;

}

void WriteFileHeader(char *filename,struct SeisFileHeader *fh)
{
	FILE *fp;
	fp = fopen(filename, "w");
	fprintf(fp, "	n1=%d\n",fh->n1);
	fprintf(fp, "	n2=%d\n",fh->n2);
	fprintf(fp, "	n3=%d\n",fh->n3);
	fprintf(fp, "	n4=%d\n",fh->n4);
	fprintf(fp, "	n5=%d\n",fh->n5);
	fprintf(fp, "	o1=%f\n",fh->o1);
	fprintf(fp, "	o2=%f\n",fh->o2);
	fprintf(fp, "	o3=%f\n",fh->o3);
	fprintf(fp, "	o4=%f\n",fh->o4);
	fprintf(fp, "	o5=%f\n",fh->o5);
	fprintf(fp, "	d1=%f\n",fh->d1);
	fprintf(fp, "	d2=%f\n",fh->d2);
	fprintf(fp, "	d3=%f\n",fh->d3);
	fprintf(fp, "	d4=%f\n",fh->d4);
	fprintf(fp, "	d5=%f\n",fh->d5);
	fprintf(fp, "	label1=\"%s\"\n",fh->label1);
	fprintf(fp, "	label2=\"%s\"\n",fh->label2);
	fprintf(fp, "	label3=\"%s\"\n",fh->label3);
	fprintf(fp, "	label4=\"%s\"\n",fh->label4);
	fprintf(fp, "	label5=\"%s\"\n",fh->label5);
	fprintf(fp, "	unit1=\"%s\"\n",fh->unit1);
	fprintf(fp, "	unit2=\"%s\"\n",fh->unit2);
	fprintf(fp, "	unit3=\"%s\"\n",fh->unit3);
	fprintf(fp, "	unit4=\"%s\"\n",fh->unit4);
	fprintf(fp, "	unit5=\"%s\"\n",fh->unit5);
	fprintf(fp, "	unit5=\"%s\"\n",fh->unit5);
	fprintf(fp, "	title=\"%s\"\n",fh->title);
	fprintf(fp, "	data_format=\"%s\"\n",fh->data_format);
	fprintf(fp, "	esize=%d\n",fh->esize);
	fprintf(fp, "	in=\"%s\"\n",fh->dname);
	fprintf(fp, "	headers=\"%s\"\n",fh->hname);
	fclose(fp);
	return;

}
