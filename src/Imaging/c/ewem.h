#ifndef _wem_h_
#define _wem_h_
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
          bool adj, bool verbose);
void elastic_extrap1f(float **mpp, float **mps1, float **mps2,
                      complex **ux_g_wx, complex **uy_g_wx, complex **uz_g_wx, 
                      complex **u_s_wx,
                      float max_source, int iw, int nw,int ifmax,int ntfft,float dw,float dkx,float dky,int nkx,int nky,
                      int nz, float oz, float dz, float gz, float sz,
                      int nmx,float omx, float dmx,
                      int nmy,float omy, float dmy,
                      int nthread,
                      float **vp,float *po_p,float **pd_p,
                      float **vs,float *po_s,float **pd_s,
                      fftwf_plan p1,fftwf_plan p2,
                      bool adj, bool verbose);
void ssop(complex *d_x,
          float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
          float **v,float *po,float **pd,
          fftwf_plan p1,fftwf_plan p2,
          bool adj,
          bool src,
          bool verbose);
void f_op(complex *m,float *d,int nw,int nt,bool adj);
void progress_msg(float progress);
float signf(float a);
float signfnonzero(float a);
int compare (const void * a, const void * b);
int omp_thread_count();
void boundary_condition(complex *d_x,int nmx,int lmx,int nmy,int lmy);
void elastic_separate_3d(complex *ux, complex *uy, complex *uz,
                      complex *up, complex *us1, complex *us2,
                      float w, 
                      float dkx, int nkx, int nmx, float omx, float dmx,
                      float dky, int nky, int nmy, float omy, float dmy,
                      float vp,float vs,
                      fftwf_plan p1,fftwf_plan p2,
                      bool inv, bool adj);
#endif
