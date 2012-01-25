

#include "ceigs.h"

#include <suitesparse/cs.h>

#include <math.h>
#include <assert.h>
#include <string.h>


/**
 * @brief Data structure to conserve LU factorization information.
 */
typedef struct cs_lu_data_s {
   css *S;     /**< Symoblic information. */
   csn *N;     /**< Numeric LU factorization. */
   double *x;  /**< Workspace. */
   int n;      /**< Dimension. */
} cs_lu_data_t;


/* LU Factorization and solving. */
static cs_lu_data_t* cs_lu_data_init( const cs *A );
static void cs_lu_data_free( cs_lu_data_t *lud );
static void cs_lu_data_solve( double *b, cs_lu_data_t *lud );

/* Common .*/
static int eigs_w_Av_lu( double *w, int n, const cs *A, const double *v );
/* DSDRV1 */
static int eigs_dsdrv1_lu( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );
/* DSDRV2 */
static void* eigs_dsdrv2_init_lu( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts );
static void eigs_dsdrv2_free_lu( void* data, const EigsOpts_t *opts );
static int eigs_dsdrv2_lu( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );
/* DSDRV3 */
static void* eigs_dsdrv3_init_lu( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts );
static void eigs_dsdrv3_free_lu( void* data, const EigsOpts_t *opts );
static int eigs_dsdrv3_lu( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );
/* DSDRV4 */
static void* eigs_dsdrv4_init_lu( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts );
static void eigs_dsdrv4_free_lu( void* data, const EigsOpts_t *opts );
static int eigs_dsdrv4_lu( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );


const EigsDriverGroup_t eigs_drv_lu = {
   .driver1 = {
      .init  = NULL,
      .free  = NULL,
      .dsdrv = eigs_dsdrv1_lu,
   },
   .driver2 = {
      .init  = eigs_dsdrv2_init_lu,
      .free  = eigs_dsdrv2_free_lu,
      .dsdrv = eigs_dsdrv2_lu,
   },
   .driver3 = {
      .init  = eigs_dsdrv3_init_lu,
      .free  = eigs_dsdrv3_free_lu,
      .dsdrv = eigs_dsdrv3_lu,
   },
   .driver4 = {
      .init  = eigs_dsdrv4_init_lu,
      .free  = eigs_dsdrv4_free_lu,
      .dsdrv = eigs_dsdrv4_lu,
   },
   .driver5 = {
      .init  = NULL,
      .free  = NULL,
      .dsdrv = NULL,
   },
   .driver6 = {
      .init  = NULL,
      .free  = NULL,
      .dsdrv = NULL,
   },
};


/**
 * @brief Initializes the LU data.
 */
static cs_lu_data_t* cs_lu_data_init( const cs *A )
{
   int order;
   double tol;
   cs_lu_data_t *lud;

   order    = 2; /* order 0:natural, 1:Chol, 2:LU, 3:QR */
   tol      = 1e-16;

   lud      = malloc( sizeof(cs_lu_data_t) );
   assert( lud != NULL );
   lud->S   = cs_sqr( order, A, 0 );
   assert( lud->S != NULL );
   lud->N   = cs_lu( A, lud->S, tol );
   assert( lud->N != NULL );
   lud->x   = cs_malloc( A->n, sizeof(double) );
   assert( lud->x != NULL );
   lud->n   = A->n;
   return lud;
}
static void cs_lu_data_free( cs_lu_data_t *lud )
{
   cs_free(  lud->x );
   cs_sfree( lud->S );
   cs_nfree( lud->N );
   free(     lud );
}
static void cs_lu_data_solve( double *b, cs_lu_data_t *lud )
{
   cs_ipvec(  lud->N->pinv, b, lud->x, lud->n );
   cs_lsolve( lud->N->L, lud->x );
   cs_usolve( lud->N->U, lud->x );
   cs_ipvec(  lud->S->q, lud->x, b, lud->n );
}


/**
 * @brief Performs w<-Av where A is a sparse matrix and w,v are both vectors.
 *    @param[out] w Output vector.
 *    @param[in] n Length of vectors.
 *    @param[in] A Sparse matrix to multiply v by.
 *    @param[in] v Vector to multiply the sparse matrix by.
 *    @return 0 on success.
 */
static int eigs_w_Av_lu( double *w, int n, const cs *A, const double *v )
{
   int err;
   memset( w, 0, n*sizeof(double) );
   err = cs_gaxpy( A, v, w );
   if (err != 1)
      fprintf( stderr, "error while running cs_gaxpy\n" );
   return !err;
}


/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_I_REGULAR.
 */
static int eigs_dsdrv1_lu( int ido, int n, double *workd, const int *ipntr, const void *data_A, const void *data_M, void *extra )
{
   const cs *A =  (const cs*) data_A;
   (void) data_M;
   (void) extra;

   /* Matrix vector multiplication ${\bf w}\leftarrow{\bf A}{\bf v}$.
    * The vector is in workd(ipntr(1)).
    * The result vector must be returned in the array workd(ipntr(2)). */
   if (ido == 1)
      eigs_w_Av_lu( &workd[ ipntr[1]-1 ], n, A, &workd[ ipntr[0]-1 ] );

   return 0;
}
static void* eigs_dsdrv2_init_lu( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts )
{
   const cs *A = (const cs*) data_A;
   (void) data_M;
   cs *C, *B, *T;
   int i;
   cs_lu_data_t *lud;

   /* Create Temoporary B matrix. */
   B = cs_spalloc( n, n, n, 1, 1 );
   for (i=0; i<n; i++)
      cs_entry( B, i, i, 1 );
   T = B;
   B = cs_compress( T );
   cs_spfree( T );

   /* C = A - sigma B */
   C = cs_add( A, B, 1.0, -opts->sigma );
   cs_spfree( B );

   /* Factorize C and keep factorization. */
   lud = cs_lu_data_init( C );
   cs_spfree( C );
   return lud;
}
static void eigs_dsdrv2_free_lu( void* data, const EigsOpts_t *opts )
{
   (void) opts;
   cs_lu_data_t *lud = (cs_lu_data_t*) data;
   cs_lu_data_free( lud );
}
static int eigs_dsdrv2_lu( int ido, int n, double *workd, const int *ipntr, const void *data_A, const void *data_M, void *extra )
{
   (void) data_A;
   (void) data_M;
   cs_lu_data_t *lud = (cs_lu_data_t*) extra;

   /* Matrix vector multiplication ${\bf w}\leftarrow{\bf A}{\bf v}$.
    * The vector is in workd(ipntr(1)).
    * The result vector must be returned in the array workd(ipntr(2)). */
   if ((ido == -1) || (ido == 1)) {
      memcpy( &workd[ ipntr[1]-1 ], &workd[ ipntr[0]-1 ], n*sizeof(double) );
      cs_lu_data_solve( &workd[ ipntr[1]-1 ], lud );
   }

   return 0;
}
/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_G_REGINVERSE.
 */
static void* eigs_dsdrv3_init_lu( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts )
{
   (void) data_A;
   const cs *M = (const cs*) data_M;
   (void) n;
   (void) opts;
   cs_lu_data_t *lud;
   lud = cs_lu_data_init( M );
   return lud;
}
static void eigs_dsdrv3_free_lu( void* data, const EigsOpts_t *opts )
{
   (void) opts;
   cs_lu_data_t *lud = (cs_lu_data_t*) data;
   cs_lu_data_free( lud );
}
static int eigs_dsdrv3_lu( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra )
{
   const cs *A = (const cs*) data_A;
   const cs *M = (const cs*) data_M;
   cs_lu_data_t *lud = (cs_lu_data_t*) extra;

   if ((ido == -1) || (ido == 1)) {
      /* y <-- A*x */
      eigs_w_Av_lu( &workd[ ipntr[1]-1 ], n, A, &workd[ ipntr[0]-1 ] );
      memcpy( &workd[ ipntr[0]-1 ], &workd[ ipntr[1]-1 ], n*sizeof(double) );

      /* M*y = A*x, solve for y */
      cs_lu_data_solve( &workd[ ipntr[1]-1 ], lud );
   }
   else if (ido == 2)
      /* y <-- M*x */
      eigs_w_Av_lu( &workd[ ipntr[1]-1 ], n, M, &workd[ ipntr[0]-1 ] );
   return 0;
}
static void* eigs_dsdrv4_init_lu( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts )
{
   const cs *A = (const cs*) data_A;
   const cs *M = (const cs*) data_M;
   (void) n;
   cs *C;
   cs_lu_data_t *lud;

   /* C = A - sigma M */
   C   = cs_add( A, M, 1.0, -opts->sigma );
   lud = cs_lu_data_init( C ); /* Only care about factorization. */
   cs_spfree( C );
   return lud;
}
static void eigs_dsdrv4_free_lu( void* data, const EigsOpts_t *opts )
{
   (void) opts;
   cs_lu_data_t *lud = (cs_lu_data_t*) data;
   cs_lu_data_free( lud );
}
/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_G_SHIFTINVERT.
 * 
 * Solving:  A*x = lambda*M*x
 */
static int eigs_dsdrv4_lu( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra )
{
   (void) data_A;
   const cs *M  = (const cs*) data_M;
   cs_lu_data_t *lud = (cs_lu_data_t*) extra;

   if (ido == -1) {
      /* y <-- inv[ A - sigma*M ] * M*x
       * input:  workd(ipntr(1))
       * output: workd(ipntr(2)) */
      eigs_w_Av_lu( &workd[ ipntr[1]-1 ], n, M, &workd[ ipntr[0]-1 ] );
      cs_lu_data_solve( &workd[ ipntr[1]-1 ], lud );
   }
   else if (ido == 1) {
      /* y <-- inv[ A - sigma*M ] * M*x
       * input:  M*x in workd(ipntr(3))
       * output: workd(ipntr(2)) */
      memcpy( &workd[ ipntr[1]-1 ], &workd[ ipntr[2]-1 ], n*sizeof(double) );
      cs_lu_data_solve( &workd[ ipntr[1]-1 ], lud );
   }
   else if (ido == 2) {
      /* y <-- M*x 
       * input:  workd(ipntr(1))
       * output: workd(ipntr(2)) */
      eigs_w_Av_lu( &workd[ ipntr[1]-1 ], n, M, &workd[ ipntr[0]-1 ] );
   }
   return 0;
}



