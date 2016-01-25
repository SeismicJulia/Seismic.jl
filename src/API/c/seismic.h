#ifndef _seismic_h_
#define _seismic_h_
#include <math.h> 
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stddef.h>
#include <sys/stat.h>
#include <complex.h> 
#include <omp.h>
#include <fftw3.h>
#ifndef MARK 
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif 
#ifndef PI
#define PI (3.141592653589793)
#endif
struct SeisHeader
{
    int tracenum;
    float o1;
    int n1;
    float d1;
    float sx;
    float sy;
    float gx;
    float gy;
    float mx;
    float my;
    float hx;
    float hy;
    float h;
    float az;
    float ang;
    int isx;
    int isy;
    int igx;
    int igy;
    int imx;
    int imy;
    int ihx;
    int ihy;
    int ih;
    int iaz;
    int iang;
    float selev;
    float gelev;
    float sstat;
    float gstat;
    int trid;      
};
struct SeisFileHeader
{
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    float o1;
    float o2;
    float o3;
    float o4;
    float o5;
    float d1;
    float d2;
    float d3;
    float d4;
    float d5;
	char label1[100];   
	char label2[100];   
	char label3[100];   
	char label4[100];   
	char label5[100];   
	char unit1[100];   
	char unit2[100];   
	char unit3[100];   
	char unit4[100];   
	char unit5[100];   
    char title[100];
    char data_format[100];
    int esize;
    char dname[100];
    char hname[100];
};
int par_read_name(int argc, char *argv[], int arg_count, char *arg, int *n);
int par_read_string(int argc, char *argv[], char *arg, char *value);
int par_read_int(int argc, char *argv[], char *arg, int *value);
int par_read_bool(int argc, char *argv[], char *arg, bool *value);
int par_read_float(int argc, char *argv[], char *arg, float *value);
void *alloc1 (size_t n1, size_t size);
void free1 (void *p);
void **alloc2 (size_t n1, size_t n2, size_t size);
void free2 (void **p);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void free3 (void ***p);
int *alloc1int(size_t n1);
void free1int(int *p);
int **alloc2int(size_t n1, size_t n2);
void free2int(int **p);
int ***alloc3int(size_t n1, size_t n2, size_t n3);
void free3int(int ***p);
float *alloc1float(size_t n1);
void free1float(float *p);
float **alloc2float(size_t n1, size_t n2);
void free2float(float **p);
float ***alloc3float(size_t n1, size_t n2, size_t n3);
void free3float(float ***p);
complex *alloc1complex(size_t n1);
void free1complex(complex *p);
complex **alloc2complex(size_t n1, size_t n2);
void free2complex(complex **p);
complex ***alloc3complex(size_t n1, size_t n2, size_t n3);
void free3complex(complex ***p);
struct SeisHeader *allocSeisHeader(size_t n1);
void freeSeisHeader(struct SeisHeader *p);
void SeisRead(char *in,float **d,struct SeisHeader *h, struct SeisFileHeader *fh);
void SeisWrite(char *out,float **d,struct SeisHeader *h, struct SeisFileHeader *fh);
void InitSeisHeader(struct SeisHeader *h, int nx);
void InitFileHeader(struct SeisFileHeader *fh);
void ReadFileHeader(char *filename,struct SeisFileHeader *fh);
void WriteFileHeader(char *filename,struct SeisFileHeader *fh);
#endif
