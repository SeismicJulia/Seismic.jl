#ifndef _wem_h_
#define _wem_h_
void wem(float **d, float **m, float **wav,
         int nt, float ot, float dt, 
         int nmx,float omx, float dmx,
         int nmy,float omy, float dmy,
         float sx,float sy,
         int nz, float oz, float dz, float gz, float sz,
         float **vel_src, float fmin, float fmax,
         int padt, int padx,
	 float damping,
         bool adj, bool pade_flag, bool verbose);
void extrap1f(float **m,complex **d_g_wx, complex **d_s_wx,float sigma,
              int iw, int ang_iw_max, int nw,int ifmax,int ntfft,float dw,float dkx,float dky,int nkx,int nky,
              int nz, float oz, float dz, float gz, float sz,
              int nmx,float omx, float dmx,
              int nmy,float omy, float dmy,
              int nthread,
              float **v,float *po,float **pd,
              fftwf_plan p1,fftwf_plan p2,
              bool adj, bool pade_flag, bool verbose);
void ssop(complex *d_x,
          float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
          float **v,float *po,float **pd,
          fftwf_plan p1,fftwf_plan p2,
          bool adj, 
          bool pade_flag,bool src,
          bool verbose);        
void pade(complex *d,
           int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,
           float w, int iz,
           float **v,float *po,float **pd,
           bool adj,bool src,
           bool verbose);
void fdop(complex *X,int nx,float dx, complex *a, bool adj);
void f_op(complex *m,float *d,int nw,int nt,bool adj);
void progress_msg(float progress);
float signf(float a);
int compare (const void * a, const void * b);
int omp_thread_count();
void boundary_condition(complex *d_x,int nmx,int lmx,int nmy,int lmy);
void compute_angles(float **angx, float **angy, float **wav,
         			int nt, float ot, float dt, 
         			int nmx,float omx, float dmx,
                    int nmy,float omy, float dmy,
                    float sx,float sy,
                    int nz, float oz, float dz, float sz,
                    float **vel_p, float fmin, float fmax,
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
#endif
