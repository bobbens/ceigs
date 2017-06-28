

#include "cs_fact.h"

#ifdef USE_UMFPACK
#include <umfpack.h>
#endif /* USE_UMFPACK */

#include <math.h>
#include <string.h>


/**
 * @brief Data structure to conserve matrix factorization.
 */
struct cs_fact_s {
   /* General information. */
   int n;               /**< Dimension. */
   /* For CSPARSE routines. */
   cs_fact_type_t type; /**< Type of factorization used. */
   css *S;              /**< Symoblic information. */
   csn *N;              /**< Factorization information.. */
   double *x;           /**< Workspace. */
#ifdef USE_UMFPACK
   /* For UMFPACK routines. */
   cs *A;               /**< Matrix. */
   void *numeric;       /**< UMFPACK numeric information. */
   int *wi;             /**< First workspace. */
   double *w;           /**< Second workspace. */
#endif /* USE_UMFPACK */
};


/*
 * Prototypes.
 */
static int cs_fact_init_cholesky( cs_fact_t *chold, const cs *A );
static int cs_fact_init_lu( cs_fact_t *lud, const cs *A );
static int cs_fact_init_qr( cs_fact_t *qrd, const cs *A );
static int cs_fact_init_umfpack( cs_fact_t *umfd, const cs *A );


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


static int cs_fact_init_umfpack( cs_fact_t *umfd, const cs *A )
{
#ifdef USE_UMFPACK
   int ret;
   void *symbolic;
   cs *T;

   /* Create matrix and sort. */
   /* Non-symmetric case. */
   if (0) {
      T        = cs_transpose( A, 1 );
      umfd->A  = cs_transpose( T, 1 );
      cs_spfree( T );
   }
   /* Symmetric case. */
   else {
      umfd->A  = cs_transpose( A, 1 );
   }
   cs_dupl( umfd->A );

   /* Generate symbolic. */
   ret = umfpack_di_symbolic( umfd->n, umfd->n, umfd->A->p, umfd->A->i, umfd->A->x, &symbolic, NULL, NULL );
   if (ret < UMFPACK_OK)
      goto err_sym;

   /* Generate numeric. */
   ret = umfpack_di_numeric(  umfd->A->p, umfd->A->i, umfd->A->x, symbolic, &umfd->numeric, NULL, NULL );
   if (ret < UMFPACK_OK)
      goto err_num;
   else if (ret > UMFPACK_OK) {
      switch (ret) {
         case UMFPACK_WARNING_singular_matrix:
            fprintf( stderr, "UMFPACK: Matrix is singular.\n" );
            break;
         case UMFPACK_WARNING_determinant_underflow:
            fprintf( stderr, "UMFPACK: Determinant underflow.\n" );
            break;
         case UMFPACK_WARNING_determinant_overflow:
            fprintf( stderr, "UMFPACK: Determinant overflow.\n" );
            break;
      }
   }

   /* Clean up symbolic. */
   umfpack_di_free_symbolic( &symbolic );

   /* Allocate buffers. */
   umfd->wi = malloc( umfd->n * sizeof(int) );
   if (umfd->wi == NULL)
      goto err_wi;
   umfd->w  = malloc( umfd->n*5 * sizeof(double) ); /* We consider iteration refinement. */
   if (umfd->w == NULL)
      goto err_w;
   umfd->x  = calloc( umfd->n, sizeof(double) );
   if (umfd->x == NULL)
      goto err_x;

   umfd->type = CS_FACT_UMFPACK;
   return 0;
err_x:
   free( umfd->w );
err_w:
   free( umfd->wi );
err_wi:
   umfpack_di_free_numeric( &umfd->numeric );
err_num:
   umfpack_di_free_symbolic( &symbolic );
err_sym:
   cs_spfree( umfd->A );
   return -1;
#else /* USE_UMFPACK */
   (void) umfd;
   (void) A;
   return -1;
#endif /* USE_UMFPACK */
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
      /* More advanced UMFPACK routines. */
#ifdef USE_UMFPACK
      case CS_FACT_NULL:
#endif /* USE_UMFPACK */
      case CS_FACT_UMFPACK:
         if (cs_fact_init_umfpack( fact, A )==0)
            return fact;
         break;

      /* Default CSPARSE routines. */
#ifndef USE_UMFPACK
      case CS_FACT_NULL:
#endif /* USE_UMFPACK */
      case CS_FACT_CHOLESKY:
         if (cs_fact_init_cholesky( fact, A )==0)
            return fact;
         /* fallthrough */
      case CS_FACT_LU:
         if (cs_fact_init_lu( fact, A )==0)
            return fact;
         /* fallthrough */
      case CS_FACT_QR:
         if (cs_fact_init_qr( fact, A )==0)
            return fact;
         break;
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
   if (fact == NULL)
      return;
   switch (fact->type) {
      case CS_FACT_NULL:
         break;
      case CS_FACT_CHOLESKY:
      case CS_FACT_LU:
      case CS_FACT_QR:
         cs_free(  fact->x );
         cs_sfree( fact->S );
         cs_nfree( fact->N );
         break;
      case CS_FACT_UMFPACK:
#ifdef USE_UMFPACK
         cs_spfree( fact->A );
         umfpack_di_free_numeric( &fact->numeric );
         free( fact->x );
         free( fact->wi );
         free( fact->w );
#endif /* USE_UMFPACK */
         break;
   }
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
#ifdef USE_UMFPACK
   int ret;
   double info[ UMFPACK_INFO ];
   int fnan, finf;
#endif /* USE_UMFPACK */
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
      case CS_FACT_UMFPACK:
#ifdef USE_UMFPACK
         ret = umfpack_di_wsolve( UMFPACK_A, /* Solving Ax=b problem. */
               fact->A->p, fact->A->i, fact->A->x,
               fact->x, b, fact->numeric, NULL, info, fact->wi, fact->w );
         if (ret == UMFPACK_WARNING_singular_matrix) {
            fprintf( stderr, "UMFPACK: wsolver Matrix singular!\n" );
            fnan = 0;
            finf = 0;
            for (k=0; k<fact->n; k++) {
               if (isnan(fact->x[k])) {
                  fact->x[k] = 0.;
                  fnan       = 1;
               }
               else if (isinf(fact->x[k])) {
                  fact->x[k] = 0.;
                  finf       = 1;
               }
            }
            if (fnan)
               fprintf( stderr, "UMFPACK: NaN values detected!\n" );
            if (finf)
               fprintf( stderr, "UMFPACK: Infinity values detected!\n" );
         }
         memcpy( b, fact->x, fact->n*sizeof(double) );
#endif /* USE_UMFPACK */
      default:
         break;
   }
}



