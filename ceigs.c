

#include "ceigs.h"

#define F77_NAME(x) x ## _
#include <suitesparse/cs.h>

#include <math.h>
#include <stdio.h>
#include <string.h>



/**
 * @brief ARPACK header for dsaupd_.
 */
extern void 
F77_NAME(dsaupd)(int *ido, char *bmat, int *n, char *which,
      int *nev, double *tol, double *resid, int *ncv,
      double *v, int *ldv, int *iparam, int *ipntr,
      double *workd, double *workl, int *lworkl,
      int *info);
/**
 * @brief ARPACK header for dseupd_.
 */
extern void
F77_NAME(dseupd)(int *rvec, char *All, int *select, double *d,
      double *v, int *ldv, double *sigma, 
      char *bmat, int *n, char *which, int *nev,
      double *tol, double *resid, int *ncv, double *v2,
      int *ldv2, int *iparam, int *ipntr, double *workd,
      double *workl, int *lworkl, int *ierr);
/**
 * @brief LAPACK header for dgttrf_
 */
extern void
F77_NAME(dgttrf)(const int* n, double* dl, double* d,
       double* du, double* du2, int* ipiv, int* info);
/**
 * @brief LAPACK header for dgttrs_
 */
extern void
F77_NAME(dgttrs)(const char* trans, const int* n, const int* nrhs,
      double* dl, double* d, double* du, double* du2,
      int* ipiv, double* b, const int* ldb, int* info);


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
   memset( w, 0, n*sizeof(double) );
   return !cs_gaxpy( A, v, w );
}


/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_I_REGULAR.
 */
static int eigs_dsdrv1_cs( int ido, int n, double *workd, const int *ipntr, const void *data_A, const void *data_B, void *extra )
{
   const cs *A =  (const cs*) data_A;
   (void) data_B;
   (void) extra;

   /* Matrix vector multiplication ${\bf w}\leftarrow{\bf A}{\bf v}$.
    * The vector is in workd(ipntr(1)).
    * The result vector must be returned in the array workd(ipntr(2)). */
   if (ido == 1)
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, A, &workd[ ipntr[0]-1 ] );

   return 0;
}
/**
 * @brief Default driver for dsaupd_ using mode EIGS_MODE_G_REGINVERSE.
 */
static int eigs_dsdrv3_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_B, void *extra )
{
   const cs *A = (const cs*) data_A;
   const cs *B = (const cs*) data_B;
   int err;
   (void) extra;

   if ((ido == -1) || (ido == 1)) {
      /* y <-- Ax */
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, A, &workd[ ipntr[0]-1 ] );
      memcpy( &workd[ ipntr[0]-1 ], &workd[ ipntr[1]-1 ], n*sizeof(double) );

      /* My = Ax, solve for y */
      err = cs_lsolve( B, &workd[ ipntr[1]-1 ] );
      if (err != 1)
         fprintf( stderr, "Error while running cs_lsolve\n" );
   }
   else if (ido == 2)
      /* y <-- Bx */
      eigs_w_Av_cs( &workd[ ipntr[1]-1 ], n, B, &workd[ ipntr[0]-1 ] );
   return 0;
}


void eigs_optsDefault( EigsOpts_t *opts )
{
   memset( opts, 0, sizeof(EigsOpts_t) );

   opts->iters = 3000;
   opts->tol   = 0.0; /* Maximum computer precision. */
}


int eigs( int n, int nev, double *lambda, double *vec, const void *data_A, const void *data_B,
      EigsOrder_t order, EigsMode_t mode, const EigsDriverGroup_t *drvlist, const EigsOpts_t *opts )
{
   int i, j;
   int ido, ncv, ldv, lworkl, info, ierr, rvec, ret;
   double *resid, *v, *workd, *workl, *d;
   double tol, sigma;
   char *which, *bmat;
   int iparam[11], ipntr[11];
   int *sel;
   void *drv_data;
   const EigsOpts_t *opts_use;
   EigsOpts_t opts_default;
   const EigsDriver_t *drv;

   /* Choose options to use. */
   if (opts == NULL) {
      eigs_optsDefault( &opts_default );
      opts_use = &opts_default;
   }
   else
      opts_use = opts;

   EigsDriverGroup_t drv_default = {
      .driver1 = {
         .init  = NULL,
         .free  = NULL,
         .dsdrv = eigs_dsdrv1_cs,
      },
      .driver2 = {
         .init  = NULL,
         .free  = NULL,
         .dsdrv = NULL,
      },
      .driver3 = {
         .init  = NULL,
         .free  = NULL,
         .dsdrv = eigs_dsdrv3_cs,
      },
      .driver4 = {
         .init  = NULL,
         .free  = NULL,
         .dsdrv = NULL,
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

   /* Use default drivers if not found. */
   if (drvlist == NULL)
      drvlist = &drv_default;

   ido   = 0; /* Initialization of the reverse communication
                 parameter. */
   tol   = opts_use->tol; /* Sets the tolerance; tol<=0 specifies 
                             machine precision */
   resid = malloc( n * sizeof(double) );
   ncv   = 4*nev; /* The largest number of basis vectors that will
                     be used in the Implicitly Restarted Arnoldi
                     Process.  Work per major iteration is
                     proportional to N*NCV*NCV. */
   if (ncv>n)
      ncv = n;

   ldv = n;
   v   = malloc( ldv*ncv * sizeof(double) );

   iparam[0] = 1;   /* Specifies the shift strategy (1->exact). */
   iparam[2] = opts_use->iters; /* Maximum number of iterations. */

   /* Get Eigenvalue output order. */
   switch (order) {
      case EIGS_ORDER_LA: which = "LA"; break;
      case EIGS_ORDER_SA: which = "SA"; break;
      case EIGS_ORDER_LM: which = "LM"; break;
      case EIGS_ORDER_SM: which = "SM"; break;
      case EIGS_ORDER_BE: which = "BE"; break;
      default:
         fprintf( stderr, "Invalid eigenvalue order type.\n" );
         return -1;
   }

   /* Get operating mode. */
   switch (mode) {
      case EIGS_MODE_I_REGULAR:
         bmat      = "I";
         iparam[6] = 1;
         drv       = &drvlist->driver1;
         break;
      case EIGS_MODE_I_SHIFTINVERT:
         bmat      = "I";
         iparam[6] = 3;
         drv       = &drvlist->driver2;
         break;
      case EIGS_MODE_G_REGINVERSE:
         bmat      = "G";
         iparam[6] = 2;
         drv       = &drvlist->driver3;
         break;
      case EIGS_MODE_G_SHIFTINVERT:
         bmat      = "G";
         iparam[6] = 3;
         drv       = &drvlist->driver4;
         break;
      case EIGS_MODE_G_BUCKLING:
         bmat      = "G";
         iparam[6] = 4;
         drv       = &drvlist->driver5;
         break;
      case EIGS_MODE_G_CAYLEY:
         bmat      = "G";
         iparam[6] = 5;
         drv       = &drvlist->driver6;
         break;
      default:
         fprintf( stderr, "Invalid Eigs mode.\n" );
         return -1;
   }
   lworkl = ncv*(ncv+8); /* Length of the workl array */
   workd  = malloc( 3*n *    sizeof(double) );
   workl  = malloc( lworkl * sizeof(double) );

   info   = 0; /* Passes convergence information out of the iteration
                  routine. */

   sel = malloc( ncv   * sizeof(int) );
   d   = malloc( 2*ncv * sizeof(double) ); /* This vector will return the eigenvalues from
                             the second routine, dseupd. */
   rvec = (vec != NULL); /* Specifies that eigenvectors should not be calculated */
   ret  = 0; /* Default return value. */

   /* Make sure we have a driver. */
   if (drv->dsdrv == NULL) {
      fprintf( stderr, "Unable to find driver to use for ARPACK.\n" );
      ret = -1;
      goto err;
   }

   /* Initialize driver. */
   if (drv->init != NULL)
      drv_data = drv->init( n );
   else
      drv_data = NULL;

   /* Main loop using ARPACK. */
   do {
      /* Reverse Communication Interface of ARPACK. */
      dsaupd_( &ido, bmat, &n, which, &nev, &tol, resid, 
               &ncv, v, &ldv, iparam, ipntr, workd, workl,
               &lworkl, &info );

      /* We'll use the driver that got set up. */
      drv->dsdrv( ido, n, workd, ipntr, data_A, data_B, drv_data );

   } while (ido != 99); /* Finish condition. */

   /* From those results, the eigenvalues and vectors are
      extracted. */
   if (info<0) {
      fprintf( stderr, "Error with dsaupd, info = %d\nCheck the documentation of dsaupd.\n\n", info );
      ret = -1;
      goto err;
   }
   else {
      dseupd_( &rvec, "All", sel, d, v, &ldv, &sigma, bmat,
               &n, which, &nev, &tol, resid, &ncv, v, &ldv,
               iparam, ipntr, workd, workl, &lworkl, &ierr);

      if (ierr!=0) {
         ret = -1;
         fprintf( stderr, "Error with dseupd, info = %d\nCheck the documentation of dseupd.\n\n", info );
         goto err;
      }
      else if (info==1)
         printf( "Maximum number of iterations reached.\n\n" );
      else if (info==3)
         fprintf( stderr, "No shifts could be applied during implicit Arnoldi update, try increasing NCV.\n\n" );

      /* Eigenvalues. */
      for (i=0; i<nev; i++)
         lambda[i] = d[i];

      /* Eigenvectors. */
      if (rvec) {
         for (i=0; i<nev; i++)
            for (j=0; j<n; j++)
               vec[ i*n+j ] = v[ i*n+j ];
      }
   }

err:
   /* Clean up after the driver. */
   if (drv->free != NULL)
      drv->free( drv_data );
   free( resid );
   free( v );
   free( workd );
   free( workl );
   free( sel );
   free( d );
   return ret;
}


void eigs_version( int *major, int *minor )
{
   *major = EIGS_VERSION_MAJOR;
   *minor = EIGS_VERSION_MINOR;
}


