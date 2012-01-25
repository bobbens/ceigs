

#ifndef _CEIGS_CS
#  define _CEIGS_CS


#include "cs_fact.h"


/* DSDRV1 */
int eigs_dsdrv1_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );
/* DSDRV2 */
void* eigs_dsdrv2_init_cs( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts, cs_fact_type_t type );
void eigs_dsdrv2_free_cs( void* data, const EigsOpts_t *opts );
int eigs_dsdrv2_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );
/* DSDRV3 */
void* eigs_dsdrv3_init_cs( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts, cs_fact_type_t type );
void eigs_dsdrv3_free_cs( void* data, const EigsOpts_t *opts );
int eigs_dsdrv3_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );
/* DSDRV4 */
void* eigs_dsdrv4_init_cs( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts, cs_fact_type_t type );
void eigs_dsdrv4_free_cs( void* data, const EigsOpts_t *opts );
int eigs_dsdrv4_cs( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );


#endif /* _CEIGS_CS */


