

#include "ceigs.h"

#include <cs.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define POW2(x)      ((x)*(x))


static int compute_error( const char *name, const char *str,
      double *derr, double *verr, int n, int nev,
      const double *d1, const double *v1,
      const double *d2, const double *v2, double tol, int def )
{
   int i, j;
   double d, v, nvec, sign;
   int ret = 0;
   d = 0.;
   v = 0.;
   nvec = (def) ? 1 : nev;
   for (i=0; i<nvec; i++) {
      d += POW2( d1[i]-d2[i] );
     
      sign = (POW2(v1[i*n+0]-v2[i*n+0]) < 1e-10) ? -1.0 : 1.0;
      for (j=0; j<n; j++)
         v += POW2( v1[i*n+j]+sign*v2[i*n+j] );
   }
   *derr = d;
   *verr = v;

   if ((d > tol) || (v > tol))
      fprintf( stderr, "%s: %s failed!\n", name, str );
   if (d > tol) {
      fprintf( stderr, "Eigenvalues do not match (err=%.3e)!\n", d );
      fprintf( stderr, "   gt:  " );
      for (i=0; i<nev; i++)
         fprintf( stderr, " %+.3e,", d2[i] );
      fprintf( stderr, "\n   res: " );
      for (i=0; i<nev; i++)
         fprintf( stderr, " %+.3e,", d1[i] );
      fprintf( stderr, "\n" );
      ret = 1;
   }
   if (v > tol) {
      fprintf( stderr, "Eigenvectors do not match (err=%.3e)!\n", v );
      for (j=0; j<nvec; j++) {
         d = 0.;
         sign = (POW2(v1[j*n+0]-v2[j*n+0]) < 1e-10) ? -1.0 : 1.0;
         for (i=0; i<n; i++)
            d += POW2( v1[j*n+i]+sign*v2[j*n+i] );
         if (d < tol)
            continue;

         fprintf( stderr, "   gt [%d]: ", j );
         for (i=0; i<n; i++)
            fprintf( stderr, " %+.3e,", v2[j*n+i] );
         fprintf( stderr, "\n   res[%d]: ", j );
         for (i=0; i<n; i++)
            fprintf( stderr, " %+.3e,", v1[j*n+i] );
         fprintf( stderr, "\n" );
      }
      ret = 1;
   }
   return ret;
}


static int eigs_test( const char *name, const char *str,
      const double *lambda_true, const double *vec_true, double tol, int def,
      int n, int nev, double *lambda, double *vec,
      const void *data_A, const void *data_M,
      EigsOrder_t order, EigsMode_t mode, const EigsDriverGroup_t *drvlist,
      const EigsOpts_t *opts )
{
   double derr, verr;
   int rret, ret = 0;

   /* Clear memory when initializing. */
   memset( lambda, 0, nev*sizeof(double) );
   memset( vec,    0, n*nev*sizeof(double) );

   /* Calculate eigenvectors. */
   rret = eigs( n, nev, lambda, vec, data_A, data_M, order, mode, drvlist, opts );
   if (rret != 0) {
      fprintf( stderr, "%s: %s failed to run!\n", name, str );
      return rret;
   }

   /* Compute error. */
   ret |= compute_error( name, str,
         &derr, &verr, n, nev, lambda, vec, lambda_true, vec_true, tol, def );
   return ret;
}


static int test_asym( const char *name, const EigsDriverGroup_t *drv, double tol )
{
   int n, nev, i, j;
   double *lambda, *vec;
   int ret = 0;

   n = 6; /* The order of the matrix */

   cs *T;
   cs *A = cs_spalloc( n, n, n*n, 1, 1 );
   for (i=0; i<n; i++)
      for (j=0; j<n; j++)
         if (!cs_entry( A, i, j, (pow(i+1, j+1) + pow(j+1, i+1)) ))
            fprintf( stderr, "Failed to add item to matrix.\n" );
   T = A;
   A = cs_compress( T );
   cs_spfree( T );

   cs *M = cs_spalloc( n, n, n, 1, 1 );
   for (i=0; i<n; i++)
      cs_entry( M, i, i, ((double) i+1) );
   T = M;
   M = cs_compress( T );
   cs_spfree( T );

   nev = 3; /* The number of values to calculate */

   /* Allocate eigenvalue and eigenvectors. */
   lambda = calloc( nev,   sizeof(double) );
   vec    = calloc( nev*n, sizeof(double) );

   /* We'll calculate Av=vd first. */
   const double l1[] = {
      9.953394836768425012e+04,
      6.141970787731389692e+02,
      -2.075006682779244827e+01
   };
   const double v1[] = {
      8.608327302135085704e-05, 1.132378335414539966e-03, 1.018141277439137994e-02, 5.679145464671306320e-02, 2.438958713544855661e-01, 9.680829426026675844e-01,
      1.144048975579629000e-02, 7.047849876936132518e-02, 2.608572146544452797e-01, 5.735162007046817889e-01, 7.404068272277498641e-01, -2.230074162262173088e-01,
      -6.001339802578492533e-02, -1.212174333083029382e-01, 2.355345514097250681e-01, 7.372730538926840493e-01, -6.091134665651405378e-01, 1.078769199806836054e-01
   };
   ret |= eigs_test( name, "Asym test Av=vd regular (dsdrv1)", l1, v1, tol, 0,
         n, nev, lambda, vec, A, NULL, EIGS_ORDER_LM, EIGS_MODE_I_REGULAR, drv, NULL );
   ret |= eigs_test( name, "Asym test Av=vd shiftinvert (dsdrv2)", l1, v1, tol, 0,
         n, nev, lambda, vec, A, NULL, EIGS_ORDER_SM, EIGS_MODE_I_SHIFTINVERT, drv, NULL );

   /* Now calculate Av = Mvd. */
   const double l2[] = {
      5.412473132453965441e-01,
      -2.984619993926789297e-02,
      -4.913844353579759350e+00
      };
   const double v2[] = {
      8.048873206468857289e-01, -2.585034763897858245e-01, -2.395618556517035702e-01, 1.064557896099580792e-01, 1.222991277965221855e-02, -6.576017055643864490e-03,
      1.726984061702139805e-01, -4.276016740271261218e-01, 4.046988721759349761e-01, -1.660254170588852940e-01, 2.403141630938674528e-02, -8.615466920500390005e-05,
      -8.976873680085523111e-02, -1.079480371895823465e-01, 1.522288623684176501e-01, 3.417208372352413259e-01, -2.885499220770642581e-01, 5.118158320569202169e-02
   };
   ret |= eigs_test( name, "Asym test Av=Mvd regular (dsdrv3)", l2, v2, tol, 0,
         n, nev, lambda, vec, A, M, EIGS_ORDER_SM, EIGS_MODE_G_REGINVERSE, drv, NULL );
   ret |= eigs_test( name, "Asym test Av=Mvd shiftinvert (dsdrv4)", l2, v2, tol, 0,
         n, nev, lambda, vec, A, M, EIGS_ORDER_LM, EIGS_MODE_G_SHIFTINVERT, drv, NULL );

   /* Clean up. */
   free( lambda );
   free( vec );
   cs_spfree( A );
   cs_spfree( M );

   return ret;
}


/**
 * @brief Tests the cholesky backend driver for symmetric positive definite matrices.
 */
static int test_sym( const char *name, const EigsDriverGroup_t *drv, double tol )
{
   int n, nev, i;
   double *lambda, *vec;
   int ret = 0;

   n = 5; /* The order of the matrix */

   cs *T;
   cs *A = cs_spalloc( n, n, n*n, 1, 1 );
   cs_entry( A, 0, 0, 1. );
   for (i=1; i<n; i++) {
      cs_entry( A, i, 0, -1. );
      cs_entry( A, 0, i, -1. );
      cs_entry( A, i, i, (double)(i+1) );
   }
   cs_entry( A, 3, 2, 1. );
   cs_entry( A, 4, 2, 1. );
   cs_entry( A, 2, 3, 1. );
   cs_entry( A, 2, 4, 1. );
   cs_entry( A, 3, 4, 2. );
   cs_entry( A, 4, 3, 2. );
   T = A;
   A = cs_compress( T );
   cs_spfree( T );

   /* This is just diagonl matrix. */
   cs *M = cs_spalloc( n, n, n, 1, 1 );
   for (i=0; i<n; i++)
      cs_entry( M, i, i, ((double) i+1) );
   T = M;
   M = cs_compress( T );
   cs_spfree( T );

   nev = 3; /* The number of values to calculate */

   /* Allocate eigenvalue and eigenvectors. */
   lambda = calloc( nev,   sizeof(double) );
   vec    = calloc( nev*n, sizeof(double) );

   /* We'll calculate Av=vd first. */
   const double l1[] = {
      7.487499930704983875,
      2.828934782420331917,
      2.396468531986137407
   };
   const double v1[] = {
       0.255653564112301834, -0.046588349401483709, -0.340340358589111491, -0.572071344927197578, -0.699552426545247852,
      -0.329043545565872841,  0.396947446945257187,  0.644073931218181128,  0.114210542611589028, -0.553432735358693084,
      -0.254721696675916931,  0.642476454309929612, -0.080823372416091277, -0.600264971072408837,  0.394322723000606445
   };
   ret |= eigs_test( name, "Sym test Av=vd regular (dsdrv1)", l1, v1, tol, 0,
         n, nev, lambda, vec, A, NULL, EIGS_ORDER_LM, EIGS_MODE_I_REGULAR, drv, NULL );
   ret |= eigs_test( name, "Sym test Av=vd shiftinvert (dsdrv2)", l1, v1, tol, 0,
         n, nev, lambda, vec, A, NULL, EIGS_ORDER_SM, EIGS_MODE_I_SHIFTINVERT, drv, NULL );

   /* Now calculate Av = Mvd. */
   const double l2[] = {
      0.787763146961615424,
      0.551298200038270458,
      0.006381569440378952
      };
   const double v2[] = {
       0.054618487217833120,  0.128673428850631683, -0.472625518114244914,  0.159718406288255427,  0.195825738820187861,
      -0.010411926777927509, -0.011602278817263223, -0.050107150349627372,  0.362488031421858536, -0.305450452541293993,
      -0.739873616712546589, -0.372312748011245365, -0.190339145144644672, -0.103589190311946536, -0.068910978382554500
   };
   ret |= eigs_test( name, "Sym test Av=Mvd regular (dsdrv3)", l2, v2, tol, 0,
         n, nev, lambda, vec, A, M, EIGS_ORDER_SM, EIGS_MODE_G_REGINVERSE, drv, NULL );
   ret |= eigs_test( name, "Sym test Av=Mvd shiftinvert (dsdrv4)", l2, v2, tol, 0,
         n, nev, lambda, vec, A, M, EIGS_ORDER_LM, EIGS_MODE_G_SHIFTINVERT, drv, NULL );

   /* Clean up. */
   free( lambda );
   free( vec );
   cs_spfree( A );
   cs_spfree( M );

   return ret;
}


/**
 * @brief Tests the cholesky backend driver for symmetric positive definite matrices.
 */
static int test_def( const char *name, const EigsDriverGroup_t *drv, double tol )
{
   int n, nev, i;
   double *lambda, *vec;
   EigsOpts_t opts;
   int ret = 0;

   n = 4; /* The order of the matrix */

   cs *T;
   cs *A = cs_spalloc( n, n, n*n, 1, 1 );
   for (i=0; i<4; i++)
      cs_entry( A, i, i, 4. );
   for (i=0; i<2; i++) {
      cs_entry( A, 0, 1+i, -2. );
      cs_entry( A, 1+i, 0, -2. );
      cs_entry( A, 3, 1+i, -2. );
      cs_entry( A, 1+i, 3, -2. );
   }
   T = A;
   A = cs_compress( T );
   cs_spfree( T );

   /* This is just diagonl matrix. */
   cs *M = cs_spalloc( n, n, n, 1, 1 );
   for (i=0; i<n; i++)
      cs_entry( M, i, i, ((double) i+1) );
   T = M;
   M = cs_compress( T );
   cs_spfree( T );

   nev = 3; /* The number of values to calculate */

   /* Options. */
   eigs_optsDefault( &opts );
   opts.sigma = 1e-5;

   /* Allocate eigenvalue and eigenvectors. */
   lambda = calloc( nev,   sizeof(double) );
   vec    = calloc( nev*n, sizeof(double) );

   /* We'll calculate Av=vd first. */
   const double l1[] = { 8.0, 4.0, 4.0 };
   double v1[] = {
       0.500000000000000000, -0.500000000000000000, -0.500000000000000000, 0.5000000000000000000,
       0.430981323461368637, -0.560584604522356922,  0.560584604522357033, -0.430981323461368526,
      -0.560584604522356700, -0.430981323461368304,  0.430981323461368582,  0.560584604522357033,
   };
   ret |= eigs_test( name, "Def test Av=vd regular (dsdrv1)", l1, v1, tol, 1,
         n, nev, lambda, vec, A, NULL, EIGS_ORDER_LM, EIGS_MODE_I_REGULAR, drv, NULL );
   ret |= eigs_test( name, "Def test Av=vd shiftinvert (dsdrv2)", l1, v1, tol, 1,
         n, nev, lambda, vec, A, NULL, EIGS_ORDER_SM, EIGS_MODE_I_SHIFTINVERT, drv, &opts );

   /* Now calculate Av = Mvd. */
   const double l2[] = {
      1.798798590723448054,
      1.460977999424757368,
      0.000000000000000036,

      };
   const double v2[] = {
      -0.347624864863106398, -0.537484369474989121,  0.154888198256868670,  0.239482252260619699,
       0.217227322371342069, -0.151922222596149653,  0.427694697909593136, -0.299116742726955376,
       0.316227766016837886,  0.316227766016837941,  0.316227766016837886,  0.316227766016837886,
   };
   ret |= eigs_test( name, "Def test Av=Mvd regular (dsdrv3)", l2, v2, tol, 0,
         n, nev, lambda, vec, A, M, EIGS_ORDER_SM, EIGS_MODE_G_REGINVERSE, drv, NULL );
   ret |= eigs_test( name, "Def test Av=Mvd shiftinvert (dsdrv4)", l2, v2, tol, 0,
         n, nev, lambda, vec, A, M, EIGS_ORDER_LM, EIGS_MODE_G_SHIFTINVERT, drv, &opts );

   /* Clean up. */
   free( lambda );
   free( vec );
   cs_spfree( A );
   cs_spfree( M );

   return ret;
}


int main( int argc, char *argv[] )
{
   (void) argc;
   (void) argv;
   int ret = 0;

   /* Symmetrical matrix tests. */
   ret |= test_sym( "default",  NULL, 1e-8 );
   ret |= test_sym( "cholesky", &eigs_drv_cholesky, 1e-10 );
   ret |= test_sym( "lu",       &eigs_drv_lu, 1e-8 );
   ret |= test_sym( "qr",       &eigs_drv_qr, 1e-5 );
#ifdef USE_UMFPACK
   ret |= test_sym( "umfpack",  &eigs_drv_umfpack, 1e-10 );
#endif /* USE_UMFPACK */

   /* Assymetrical matrix tests. */
   ret |= test_asym( "default",  NULL, 1e-8 );
   ret |= test_asym( "cholesky", &eigs_drv_cholesky, 1e-8 );
   ret |= test_asym( "lu",       &eigs_drv_lu, 1e-8 );
   ret |= test_asym( "qr",       &eigs_drv_qr, 1e-5 );
#ifdef USE_UMFPACK
   ret |= test_asym( "umfpack",  &eigs_drv_umfpack, 1e-10 );
#endif /* USE_UMFPACK */

   /* Rank deficient matrix tests. */
   ret |= test_def( "default",  NULL, 1e-8 );
   ret |= test_def( "cholesky", &eigs_drv_cholesky, 1e-8 );
   ret |= test_def( "lu",       &eigs_drv_lu, 1e-8 );
   ret |= test_def( "qr",       &eigs_drv_qr, 1e-5 );
#ifdef USE_UMFPACK
   ret |= test_def( "umfpack",  &eigs_drv_umfpack, 1e-10 );
#endif /* USE_UMFPACK */

   if (ret)
      fprintf( stderr, "ceigs test failed!!\n" );
   else
      printf( "ceigs test success!!\n" );
   return ret;
}



