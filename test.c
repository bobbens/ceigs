

#include "ceigs.h"

#include <suitesparse/cs.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define POW2(x)      ((x)*(x))

static int compute_error( double *derr, double *verr, int n, int nev,
      const double *d1, const double *v1,
      const double *d2, const double *v2 )
{
   int i, j;
   double d, v, sign;
   int ret = 0;
   d = 0.;
   v = 0.;
   for (i=0; i<nev; i++) {
      d += POW2( d1[i]-d2[i] );
     
      sign = (POW2(v1[i*n+0]-v2[i*n+0]) < 1e-10) ? -1.0 : 1.0;
      for (j=0; j<n; j++)
         v += POW2( v1[i*n+j]+sign*v2[i*n+j] );
   }
   *derr = d;
   *verr = v;

   if (d > 1e-10) {
      fprintf( stderr, "Eigenvalues do not match (err=%.3e)!\n", d );
      ret = 1;
   }
   if (v > 1e-10) {
      fprintf( stderr, "Eigenvectors do not match (err=%.3e)!\n", v );
      ret = 1;
   }
   return ret;
}

int main( int argc, char *argv[] )
{
   (void) argc;
   (void) argv;
   int n, nev, i, j;
   double *lambda, *vec;
   double derr, verr;
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

   cs *B = cs_spalloc( n, n, n, 1, 1 );
   for (i=0; i<n; i++)
      cs_entry( B, i, i, ((double) i+1) );
   T = B;
   B = cs_compress( T );
   cs_spfree( T );

   nev = 3; /* The number of values to calculate */

   /* Allocate eigenvalue and eigenvectors. */
   lambda = calloc( nev,   sizeof(double) );
   vec    = calloc( nev*n, sizeof(double) );

   /* We'll calculate Av=vd first. */
   eigs( n, nev, lambda, vec, A, NULL, EIGS_ORDER_LM, EIGS_MODE_I_REGULAR, NULL, NULL );
   const double l1[] = { -2.075006682779244827e+01, 6.141970787731389692e+02, 9.953394836768425012e+04 };
   const double v1[] = {
      -6.001339802578492533e-02, -1.212174333083029382e-01, 2.355345514097250681e-01, 7.372730538926840493e-01, -6.091134665651405378e-01, 1.078769199806836054e-01,
      1.144048975579629000e-02, 7.047849876936132518e-02, 2.608572146544452797e-01, 5.735162007046817889e-01, 7.404068272277498641e-01, -2.230074162262173088e-01,
      8.608327302135085704e-05, 1.132378335414539966e-03, 1.018141277439137994e-02, 5.679145464671306320e-02, 2.438958713544855661e-01, 9.680829426026675844e-01
   };
   ret |= compute_error( &derr, &verr, n, nev, lambda, vec, l1, v1 );

   /* Now calculate Av = Mvd. */
   eigs( n, nev, lambda, vec, A, B, EIGS_ORDER_SM, EIGS_MODE_G_REGINVERSE, NULL, NULL );
   const double l2[] = { -4.913844353579759350e+00, -2.984619993926789297e-02, 5.412473132453965441e-01 };
   const double v2[] = {
      -8.976873680085523111e-02, -1.079480371895823465e-01, 1.522288623684176501e-01, 3.417208372352413259e-01, -2.885499220770642581e-01, 5.118158320569202169e-02,
      1.726984061702139805e-01, -4.276016740271261218e-01, 4.046988721759349761e-01, -1.660254170588852940e-01, 2.403141630938674528e-02, -8.615466920500390005e-05,
      8.048873206468857289e-01, -2.585034763897858245e-01, -2.395618556517035702e-01, 1.064557896099580792e-01, 1.222991277965221855e-02, -6.576017055643864490e-03
   };
   ret |= compute_error( &derr, &verr, n, nev, lambda, vec, l2, v2 );

   /* Clean up. */
   free( lambda );
   free( vec );
   cs_spfree( A );
   cs_spfree( B );

   if (ret)
      fprintf( stderr, "ceigs test failed!!\n" );
   else
      printf( "ceigs test success!!\n" );
   return ret;
}



