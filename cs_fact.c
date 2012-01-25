

#include "cs_fact.h"

#include <math.h>
#include <string.h>


/**
 * @brief Data structure to conserve matrix factorization.
 */
struct cs_fact_s {
   cs_fact_type_t type; /**< Type of factorization used. */
   css *S;              /**< Symoblic information. */
   csn *N;              /**< Factorization information.. */
   double *x;           /**< Workspace. */
   int n;               /**< Dimension. */
};


/*
 * Prototypes.
 */
static int cs_fact_init_cholesky( cs_fact_t *chold, const cs *A );
static int cs_fact_init_lu( cs_fact_t *lud, const cs *A );
static int cs_fact_init_qr( cs_fact_t *qrd, const cs *A );


/**
 * @brief Performs Cholesky factorization on a matrix.
 *
 *    @param[in,out] chold Cholesky factorization of the matrix.
 *    @param[in] A Matrix to factorize.
 *    @return 0 on success.
 */
static int cs_fact_init_cholesky( cs_fact_t *chold, const cs *A )
{
   int order = 1; /* order 0:natural, 1:Chol, 2:LU, 3:QR */

   chold->S = cs_schol( order, A );
   if (chold->S == NULL)
      goto err_S;
   chold->N = cs_chol( A, chold->S );
   if (chold->N == NULL)
      goto err_N;
   chold->x = cs_malloc( A->n, sizeof(double) );
   if (chold->x == NULL)
      goto err_x;

   chold->type = CS_FACT_CHOLESKY;
   return 0;
err_x:
   cs_nfree( chold->N );
err_N:
   cs_sfree( chold->S );
err_S:
   return -1;
}


/**
 * @brief Performs LU factorization on a matrix.
 *
 *    @param[in,out] lud LU factorization of the matrix.
 *    @param[in] A Matrix to factorize.
 *    @return 0 on success.
 */
static int cs_fact_init_lu( cs_fact_t *lud, const cs *A )
{
   int order  = 2; /* order 0:natural, 1:Chol, 2:LU, 3:QR */
   double tol = 1e-16;

   lud->S   = cs_sqr( order, A, 0 );
   if (lud->S == NULL)
      goto err_S;
   lud->N   = cs_lu( A, lud->S, tol );
   if (lud->N == NULL)
      goto err_N;
   lud->x   = cs_malloc( A->n, sizeof(double) );
   if (lud->x == NULL)
      goto err_x;

   lud->type = CS_FACT_LU;
   return 0;
err_x:
   cs_nfree( lud->N );
err_N:
   cs_sfree( lud->S );
err_S:
   return -1;
}


/**
 * @brief Performs QR factorization on a matrix.
 *
 *    @param[in,out] qrd QR factorization of the matrix.
 *    @param[in] A Matrix to factorize.
 *    @return 0 on success.
 */
static int cs_fact_init_qr( cs_fact_t *qrd, const cs *A )
{
   int order  = 3; /* order 0:natural, 1:Chol, 2:LU, 3:QR */

   qrd->S   = cs_sqr( order, A, 1 );
   if (qrd->S == NULL)
      goto err_S;
   qrd->N   = cs_qr( A, qrd->S );
   if (qrd->N == NULL)
      goto err_N;
   qrd->x   = cs_malloc( qrd->S->m2, sizeof(double) );
   if (qrd->x == NULL)
      goto err_x;

   qrd->type = CS_FACT_QR;
   return 0;
err_x:
   cs_nfree( qrd->N );
err_N:
   cs_sfree( qrd->S );
err_S:
   return -1;
}


/**
 * @brief Attempts to factorize a sparse matrix.
 *
 *    @param[in] A Matrix to factorize.
 *    @return Factorization of the matrix or NULL on error.
 */
cs_fact_t* cs_fact_init_type( const cs *A, cs_fact_type_t type )
{
   cs_fact_t *fact;

   /* Set up. */
   fact        = malloc( sizeof(cs_fact_t) );
   fact->type  = CS_FACT_NULL;
   fact->n     = A->n;

   /* Try in order Cholesky, LU and QR. */
   switch (type) {
      case CS_FACT_NULL:
      case CS_FACT_CHOLESKY:
         if (cs_fact_init_cholesky( fact, A )==0)
            return fact;
      case CS_FACT_LU:
         if (cs_fact_init_lu( fact, A )==0)
            return fact;
      case CS_FACT_QR:
         if (cs_fact_init_qr( fact, A )==0)
            return fact;
   }

   /* Failed to generate anything. */
   free( fact );
   return NULL;
}


/**
 * @brief Cleans up a sparse matrix factorization.
 *
 *    @param[in] fact Factorization to clean up.
 */
void cs_fact_free( cs_fact_t *fact )
{
   cs_free(  fact->x );
   cs_sfree( fact->S );
   cs_nfree( fact->N );
   free(     fact );
}


/**
 * @brief Solves a sparse matrix factorization.
 *
 * This solves Ax = b.
 *
 *    @param[in,out] b Input right side vector that gets set to the output solution.
 *    @param[in] fact Factorization of the sparse matrix.
 */
void cs_fact_solve( double *b, cs_fact_t *fact )
{
   int k;
   switch (fact->type) {
      case CS_FACT_CHOLESKY:
         cs_ipvec(   fact->S->pinv, b, fact->x, fact->n );
         cs_lsolve(  fact->N->L, fact->x );
         cs_ltsolve( fact->N->L, fact->x );
         cs_pvec(    fact->S->pinv, fact->x, b, fact->n );
         break;
      case CS_FACT_LU:
         cs_ipvec(  fact->N->pinv, b, fact->x, fact->n );
         cs_lsolve( fact->N->L, fact->x );
         cs_usolve( fact->N->U, fact->x );
         cs_ipvec(  fact->S->q, fact->x, b, fact->n );
         break;
      case CS_FACT_QR:
         cs_ipvec(  fact->S->pinv, b, fact->x, fact->n );
         for (k=0; k<fact->n; k++)
            cs_happly( fact->N->L, k, fact->N->B[k], fact->x );
         cs_usolve( fact->N->U, fact->x );
         cs_ipvec(  fact->S->q, fact->x, b, fact->n );
         break;
      default:
         break;
   }
}



