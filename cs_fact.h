


#ifndef _CS_FACT_H
#  define _CS_FACT_H


#include <cs.h>


/**
 * @brief Enumeration of possible matrix factorization.
 */
typedef enum cs_fact_type_e {
   CS_FACT_NULL,     /**< No factorization. */
   CS_FACT_CHOLESKY, /**< Cholesky factorization. */
   CS_FACT_LU,       /**< LU factorization. */
   CS_FACT_QR,       /**< QR factorization. */
   CS_FACT_UMFPACK   /**< Using UMFPACK library. */
} cs_fact_type_t;


struct cs_fact_s;
typedef struct cs_fact_s cs_fact_t;


#define cs_fact_init( A ) \
cs_fact_init_type( A, CS_FACT_CHOLESKY )
cs_fact_t* cs_fact_init_type( const cs *A, cs_fact_type_t type );
void cs_fact_free( cs_fact_t *fact );
void cs_fact_solve( double *b, cs_fact_t *fact );


#endif /* _CS_FACT_H */


