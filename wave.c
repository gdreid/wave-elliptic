//=============================================================================
// headers
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd_w.h>
#include <math.h>
#include <mpi.h>
#include <sys/timeb.h>

#include "init_f.h"
#include "init_f_t.h"

#include "res_f.h"
#include "res_f_t.h"
#include "res_psi_t0.h"

#include "u_f.h"
#include "u_f_t.h"

#include "lop_psi_t0.h"
#include "relax_psi_t0.h"
#include "res_mg_psi_t0.h"

#include "num.h"


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//=============================================================================
// function prototypes
//=============================================================================
void set_gfns(void);
void ldptr(void);
void const_f(real *f, real c);
void zero(real *f);
void copy_gf(real *f, real*g);

int wave_id(void);
void wave_var_pre_init(char *pfile);
void wave_var_post_init(char *pfile);
void wave_AMRH_var_clear(void);
void wave_free_data(void);
void wave_t0_cnst_data(void);
void wave_pre_io_calc(void);
real wave_evo_residual(void);
real wave_MG_residual(void);
void wave_evolve(int iter, int *ifc_mask);
void wave_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised);
void wave_fill_bh_bboxes(real *bbox, int *num, int max_num);
void wave_post_tstep(int L);
real wave_MG_relax(void);
void wave_L_op(void);
void wave_scale_tre(void);
void wave_post_regrid(void);
void elapsed_time(void);

real nbs_MG_residual(void);
real nbs_MG_relax(void);
void nbs_L_op(void);

//=============================================================================
// parameters
//=============================================================================

real famp, r0, z0, delr, delz, idsignum;

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *R, *Z;
real res, myzero;

//=============================================================================
// multigrid variables
//=============================================================================
real *f,     *f_t,     *psi;
int   f_gfn,  f_t_gfn,  psi_gfn;

real *psi_t0,     *psi_t0_res,     *psi_t0_lop,     *psi_t0_rhs;
int   psi_t0_gfn,  psi_t0_res_gfn,  psi_t0_lop_gfn,  psi_t0_rhs_gfn; 

real *mg_w1;
int  mg_w1_gfn;

//=============================================================================
// amr variables
//=============================================================================

real *n_f,    *np1_f;
int   n_f_gfn, np1_f_gfn;

real *n_f_t,    *np1_f_t;
int   n_f_t_gfn, np1_f_t_gfn;

real *n_psi,    *np1_psi;
int   n_psi_gfn, np1_psi_gfn;

real *n_psi_t0,    *np1_psi_t0;
int   n_psi_t0_gfn, np1_psi_t0_gfn;

real *w1,     *w2;
int   w1_gfn,  w2_gfn;

real *mask,    *mask_mg;
int  mask_gfn, mask_mg_gfn;

//=============================================================================
// other variables
//=============================================================================

real res_f, res_f_t, res_psi, res_psi_t0;

int shape[3], ghost_width[6], NR, NZ, phys_bdy[6], size, g_rank, dim, g_L, ngfs;

real base_bbox[6], bbox[6], hR, hZ, ht, t;

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((f_gfn     = PAMR_get_gfn("f",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error f",    0);
    if ((n_f_gfn   = PAMR_get_gfn("f",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error n_f",  0);
    if ((np1_f_gfn = PAMR_get_gfn("f",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error np1_f",0);
    
    if ((f_t_gfn     = PAMR_get_gfn("f_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error f_t",    0);
    if ((n_f_t_gfn   = PAMR_get_gfn("f_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error n_f_t",  0);
    if ((np1_f_t_gfn = PAMR_get_gfn("f_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error np1_f_t",0);
    
    if ((w1_gfn  = PAMR_get_gfn("w1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error w1",0);
    if ((w2_gfn  = PAMR_get_gfn("w2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error w2",0);
    
    if ((psi_gfn     = PAMR_get_gfn("psi",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error psi",    0);
    if ((n_psi_gfn   = PAMR_get_gfn("psi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error n_psi",  0);
    if ((np1_psi_gfn = PAMR_get_gfn("psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error np1_psi",0);
    
    if ((psi_t0_gfn     = PAMR_get_gfn("psi_t0",    PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error psi_t0",    0);
    if ((psi_t0_res_gfn = PAMR_get_gfn("psi_t0_res",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error psi_t0_res",0);
    if ((psi_t0_lop_gfn = PAMR_get_gfn("psi_t0_lop",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error psi_t0_lop",0);
    if ((psi_t0_rhs_gfn = PAMR_get_gfn("psi_t0_rhs",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error psi_t0_rhs",0);
    
    if ((mask_mg_gfn=PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn=PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    
    if ((mg_w1_gfn  = PAMR_get_gfn("mg_w1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error mg_w1",0);
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr(void)
{
   real dx0[3];
   real *x0[3];
   real *gfs[PAMR_MAX_GFNS];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

/* 
   if (!(PAMR_get_g_attribs(&g_rank,&dim,shape,bbox,ghost_width,&t,&ngfs,x0,gfs))) 
      AMRD_stop("ldptr: PAMR_get_g_attribs failed\n","");
*/

   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_dim(&dim);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_t(&t);
   PAMR_get_g_ngfs(&ngfs);
   PAMR_get_g_x(x0);
   PAMR_get_g_gfs(gfs);

   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&ht);
   hZ=dx0[0];
   hR=dx0[1];
   
   Z=x0[0];
   R=x0[1];  
   hZ = Z[1] - Z[0];
   hR = R[1] - R[0];
   
   NZ=shape[0];
   NR=shape[1]; 
   size = NZ*NR;

   if ((bbox[0]-base_bbox[0])<hZ/2.0) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<hZ/2.0) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<hR/2.0) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<hR/2.0) phys_bdy[3]=1; else phys_bdy[3]=0;

   PAMR_get_g_gfs(gfs);

   f     = gfs[f_gfn     -1];
   n_f   = gfs[n_f_gfn   -1];
   np1_f = gfs[np1_f_gfn -1];
   
   f_t     = gfs[f_t_gfn     -1];
   n_f_t   = gfs[n_f_t_gfn   -1];
   np1_f_t = gfs[np1_f_t_gfn -1];
   
   psi     = gfs[psi_gfn     -1];
   n_psi   = gfs[n_psi_gfn   -1];
   np1_psi = gfs[np1_psi_gfn -1];
   
   psi_t0     = gfs[psi_t0_gfn     -1]; 
   psi_t0_res = gfs[psi_t0_res_gfn -1]; 
   psi_t0_lop = gfs[psi_t0_lop_gfn -1]; 
   psi_t0_rhs = gfs[psi_t0_rhs_gfn -1]; 
   
   mask=gfs[mask_gfn-1];
   mask_mg=gfs[mask_mg_gfn-1];

   w1    = gfs[w1_gfn    -1];
   w2    = gfs[w2_gfn    -1];
   mg_w1 = gfs[mg_w1_gfn -1];
}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<NR*NZ; i++) f[i]=c;
}

void zero(real *f)
{
   const_f(f,0);
}

double norm(real *f)
{
   int i;
   real norm = 0;
   for (i = 0; i < NZ*NR; i++) norm += f[i]*f[i];
   return sqrt(norm/(double)(NZ*NR)); 
}

void copy_gf(real *f, real*g)
{
   int i;
   for (i = 0; i < NZ*NR; i++) g[i] = f[i];
}


//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int wave_id(void)
{
   if( g_rank == 0 ) elapsed_time();
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void wave_var_pre_init(char *pfile)
{
   return;
}

void wave_var_post_init(char *pfile)
{
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
      system("date > date.dat");
      printf("===================================================================\n");
      printf("Reading wave parameters:\n\n");
   }

   famp = r0 = z0 = delr = delz = idsignum = 0;
   myzero = 0;

   AMRD_real_param(pfile,"famp",     &famp,1);
   AMRD_real_param(pfile,"r0",       &r0,1);
   AMRD_real_param(pfile,"z0",       &z0,1);
   AMRD_real_param(pfile,"delr",     &delr,1);
   AMRD_real_param(pfile,"delz",     &delz,1);
   AMRD_real_param(pfile,"idsignum", &idsignum,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void wave_AMRH_var_clear(void)
{
   ldptr();

   zero(n_f); 
   zero(np1_f);
   
   zero(n_f_t);
   zero(np1_f_t);
   
   zero(n_psi);
   zero(np1_psi);

   zero(w1);
   zero(w2);

   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2)
//=============================================================================
void wave_free_data(void)
{
   ldptr();
   
   init_f_(R,Z,&NZ,&NR,&delr,&delz,&famp,&r0,&z0,n_f);
   init_f_t_(R,Z,&NZ,&NR,&delr,&delz,&famp,&idsignum,&r0,&z0,n_f_t);

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//=============================================================================
void wave_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void wave_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration.
//=============================================================================
real wave_evo_residual(void)
{
   ldptr();

   res = 0;
   res_f_(R,Z,n_f,n_f_t,np1_f,np1_f_t,&NZ,&NR,&ht,&myzero,phys_bdy,&res_f);
   res = MAX(res_f,res);
   
   res_f_t_(R,Z,n_f,n_f_t,np1_f,np1_f_t,&NZ,&NR,&hR,&hZ,&ht,&myzero,phys_bdy,&res_f_t);
   res = MAX(res_f_t,res);
   
   return res;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real wave_MG_residual(void)
{
   ldptr();
   
   zero(psi_t0_res);
   res_mg_psi_t0_(R,Z,f,psi_t0,psi_t0_rhs,&NZ,&NR,&hR,&hZ,&myzero,phys_bdy,mask_mg,psi_t0_res);   
   
   return norm(psi_t0_res);
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real wave_MG_relax(void)
{
   int i, j;
   ldptr();
   
   copy_gf(psi_t0,mg_w1);
   relax_psi_t0_(R,Z,f,psi_t0,psi_t0_rhs,&NZ,&NR,&hR,&hZ,&myzero,phys_bdy,mask_mg,mg_w1);
   copy_gf(mg_w1,psi_t0);
   
   zero(mg_w1);
   res_mg_psi_t0_(R,Z,f,psi_t0,psi_t0_rhs,&NZ,&NR,&hR,&hZ,&myzero,phys_bdy,mask_mg,mg_w1);
   
   if (NZ < 10 && NR < 10) {
   //printf("%d %d\n", NZ,NR);
   }
   
   return norm(mg_w1);
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void wave_L_op(void)
{
   ldptr();
   
   zero(psi_t0_lop);
   lop_psi_t0_(R,Z,f,psi_t0,&NZ,&NR,&hR,&hZ,&myzero,phys_bdy,mask_mg,psi_t0_lop);

   return;
}

//=============================================================================
// Performs 1 iteration of the evolution equations (new version of AMRD 
// adds additional parameter, *ifc_mask).
//============================================================================= 
void wave_evolve(int iter, int *ifc_mask)
{
   ldptr();
   
   u_f_(R,Z,n_f,n_f_t,np1_f,np1_f_t,&NZ,&NR,&ht,&myzero,phys_bdy,np1_f);
   u_f_t_(R,Z,n_f,n_f_t,np1_f,np1_f_t,&NZ,&NR,&hR,&hZ,&ht,&myzero,phys_bdy,np1_f_t);

   return;
}


//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real nbs_MG_residual(void)
{
   real norm;

   ldptr();

   residual_(psi_t0_res,psi_t0_rhs,psi_t0,f,f,mask_mg,Z,R,&norm,&NZ,&NR);

   return norm;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real nbs_MG_relax(void)
{
   real normal;

   ldptr();

   relax_(psi_t0,psi_t0_rhs,f,f,mask_mg,phys_bdy,Z,R,&normal,&NZ,&NR);
   
   //printf("%5.5e\n",(double) (norm(mask_mg)));

   return normal;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void nbs_L_op(void)
{
   ldptr();

   lop_(psi_t0_lop,psi_t0,f,f,mask_mg,Z,R,&NZ,&NR);

   return;
}


//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
//=============================================================================
void wave_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
   return;
}

//=============================================================================
void wave_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
   return;
}

//=============================================================================
void wave_post_tstep(int L)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void wave_scale_tre(void)
{
   return;
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void wave_post_regrid(void)
{
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd(argc,argv,&wave_id,&wave_var_pre_init,
        &wave_var_post_init, &wave_AMRH_var_clear,
        &wave_free_data, &wave_t0_cnst_data,
        &wave_evo_residual, &wave_MG_residual,
        &wave_evolve, &wave_MG_relax, &wave_L_op, 
        &wave_pre_io_calc, &wave_scale_tre, 
        &wave_post_regrid, &wave_post_tstep,
        &wave_fill_ex_mask, &wave_fill_bh_bboxes);
   if (my_rank==0) elapsed_time();
}

//=============================================================================
// Maintains/reports elapsed wall-clock time.
//=============================================================================
void elapsed_time(void) {
   static int    first = 1;
   struct        timeb t;
   static double msinit;
   double        mscurr, mselapsed;

   ftime(&t);
   mscurr = 1000.0 * t.time + t.millitm;
   if( first ) {
      msinit = mscurr;
      first = 0;
   }
	mselapsed = mscurr - msinit;
   printf("elapsed_time: Seconds since initial call: %12.3f\n",
         mselapsed / 1000.0);
}
