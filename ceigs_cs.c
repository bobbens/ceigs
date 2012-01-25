

#include "ceigs.h"
#include "ceigs_cs.h"

#include <suitesparse/cs.h>

#include <math.h>
#include <assert.h>
#include <string.h>


/* Common .*/
static int eigs_w_Av_cs( double *w, int n, const cs *A, const double *v );


/**
 * @brief Performs w<-Av where A is a sparse matrix and w,v are both vectors.
 *    @param[out] w Output vector.
 *    @param[in] n Length of vectors.
 *    @param[in] A Sparse matrix to multiply v by.
 *    @param[in] v Vector to multiply the sparse matrix by.
 *    @return 0 on success.
 */
static int eigs_w_Av_cs( double *w, int n, const cs *A, const double *v )
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
int eigs_dsdrv1_cs( int ido, int n, double *workd, const int *ipntr, const void *data_A, const void *data_M, void *extra )
{
   const cs *A =  (const cs*) data_A;
   (void) data_M;
   (void) extra;

   /* Matrix vector multiplication ${\bf w}\leftarrow{\bf A}{\bf v}$.
    * The vector is in workd(ipntr(1)).
    * The result vector must be returned in the array workd(ipntr(2)). */
   if (ido == 1)
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, A, &workd[ ipntr[0]-1 ] );

   return 0;
}
void* eigs_dsdrv2_init_cs( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts, cs_fact_type_t type )
{
   const cs *A = (const cs*) data_A;
   (void) data_M;
   cs *C, *B, *T;
   int i;
   cs_fact_t *fact;

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
   fact = cs_fact_init_type( C, type );
   cs_spfree( C );
   return fact;
}
void eigs_dsdrv2_free_cs( void* data, const EigsOpts_t *opts )
{
   (void) opts;
   cs_fact_t *fact = (cs_fact_t*) data;
   cs_fact_free( fact );
}
int eigs_dsdrv2_cs( int ido, int n, double *workd, const int *ipntr, const void *data_A, const void *data_M, void *extra )
{
   (void) data_A;
   (void) data_M;
   cs_fact_t *fact = (cs_fact_t*) extra;

   /* Matrix vector multiplication ${\bf w}\leftarrow{\bf A}{\bf v}$.
    * The vector is in workd(ipntr(1)).
    * The result vector must be returned in the array workd(ipntr(2)). */
   if ((ido == -1) || (ido == 1)) {
      memcpy( &workd[ ipntr[1]-1 ], &workd[ ipntr[0]-1 ], n*sizeof(double) );
      cs_fact_solve( &workd[ ipntr[1]-1 ], fact );
   }

   return 0;
}
/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_G_REGINVERSE.
 */
void* eigs_dsdrv3_init_cs( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts, cs_fact_type_t type )
{
   (void) data_A;
   const cs *M = (const cs*) data_M;
   (void) n;
   (void) opts;
   cs_fact_t *fact;
   fact = cs_fact_init_type( M, type );
   return fact;
}
void eigs_dsdrv3_free_cs( void* data, const EigsOpts_t *opts )
{
   (void) opts;
   cs_fact_t *fact = (cs_fact_t*) data;
   cs_fact_free( fact );
}
int eigs_dsdrv3_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra )
{
   const cs *A = (const cs*) data_A;
   const cs *M = (const cs*) data_M;
   cs_fact_t *fact = (cs_fact_t*) extra;

   if ((ido == -1) || (ido == 1)) {
      /* y <-- A*x */
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, A, &workd[ ipntr[0]-1 ] );
      memcpy( &workd[ ipntr[0]-1 ], &workd[ ipntr[1]-1 ], n*sizeof(double) );

      /* M*y = A*x, solve for y */
      cs_fact_solve( &workd[ ipntr[1]-1 ], fact );
   }
   else if (ido == 2)
      /* y <-- M*x */
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, M, &workd[ ipntr[0]-1 ] );
   return 0;
}
void* eigs_dsdrv4_init_cs( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts, cs_fact_type_t type )
{
   const cs *A = (const cs*) data_A;
   const cs *M = (const cs*) data_M;
   (void) n;
   cs *C;
   cs_fact_t *fact;

   /* C = A - sigma M */
   C   = cs_add( A, M, 1.0, -opts->sigma );
   fact = cs_fact_init_type( C, type ); /* Only care about factorization. */
   cs_spfree( C );
   return fact;
}
void eigs_dsdrv4_free_cs( void* data, const EigsOpts_t *opts )
{
   (void) opts;
   cs_fact_t *fact = (cs_fact_t*) data;
   cs_fact_free( fact );
}
/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_G_SHIFTINVERT.
 * 
 * Solving:  A*x = lambda*M*x
 */
int eigs_dsdrv4_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra )
{
   (void) data_A;
   const cs *M  = (const cs*) data_M;
   cs_fact_t *fact = (cs_fact_t*) extra;

   if (ido == -1) {
      /* y <-- inv[ A - sigma*M ] * M*x
       * input:  workd(ipntr(1))
       * output: workd(ipntr(2)) */
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, M, &workd[ ipntr[0]-1 ] );
      cs_fact_solve( &workd[ ipntr[1]-1 ], fact );
   }
   else if (ido == 1) {
      /* y <-- inv[ A - sigma*M ] * M*x
       * input:  M*x in workd(ipntr(3))
       * output: workd(ipntr(2)) */
      memcpy( &workd[ ipntr[1]-1 ], &workd[ ipntr[2]-1 ], n*sizeof(double) );
      cs_fact_solve( &workd[ ipntr[1]-1 ], fact );
   }
   else if (ido == 2) {
      /* y <-- M*x 
       * input:  workd(ipntr(1))
       * output: workd(ipntr(2)) */
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, M, &workd[ ipntr[0]-1 ] );
   }
   return 0;
}



