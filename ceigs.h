
#ifndef CEIG_SPARSE_H
#  define CEIG_SPARSE_H

/**
 * @mainpage ceigs doxygen documentation
 * @author Edgar Simo-Serra <esimo@iri.upc.edu>
 * @version 1.0
 * @date December 2011
 *
 * @section License
 *
 @verbatim
    Copyright 2011, 2012 Edgar Simo-Serra

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 @endverbatim
 *
 *
 * @section Overview
 *
 * This is a simple C frontend for ARPACK. This allows easy access to
 * calculating a subset of eigenvectors and eigenvalues of sparse matrices.
 * Specifically it can solve two problems:
 * - \f$Av = vd\f$
 * - \f$Av = Mvd\f$
 *
 * Where \f$A, M\f$ are sparse matrices, \f$v\f$ is the subset of
 * eigenvectors and \f$d\f$ is the diagonal matrix of eigenvalues.
 *
 * Currently only symmetric real matrices are supported.
 *
 *
 * @section Usage
 *
 * ARPACK uses a Reverse Communication Interface to be extremely
 * flexible. This wrapper aims to keep that flexibility while also making
 * it much easier to use. A small example of usage would be:
 *
 * @code
 * #include <ceigs.h>
 * #include <suitesparse/csparse.h>
 * #include <stdio.h>
 *
 * int main( int argc, char *argv[] )
 * {
 *    cs *A; // Matrix to work with
 *    int n; // Order of the matrix
 *    int nev; // Number of eigenvectors to calculate, nev < n+1
 *    A = cs_spalloc( ... );
 *    // Fill A here
 *
 *    // Output data
 *    double *v, *d;
 *    v = malloc( n * nev * sizeof(double) );
 *    d = malloc( nev * sizeof(double) );
 *
 *    // Actual algorithm
 *    int ret;
 *    ret = eigs( n, nev, d, v, // Output and problem information
 *          A, NULL, // We are solving Av = vd, so no need for second matrix
 *          EIGS_ORDER_SM, // Order to get eignvalues in
 *          EIGS_MODE_I_REGULAR, // Mode of operation of ARPACK
 *          NULL, NULL ); // Use default driver (csparse) and default options
 *    if (ret != 0)
 *       fprintf( stderr, "An error occurred while running eigs!\n" );
 *
 *    // Print some output
 *    int i;
 *    for (i=0; i<nev; i++)
 *       printf( "Eigenvalue  %d: %f\n"
 *               "Eigenvector %d: ( %f, %f, %f, %f, ... )\n",
 *               i, d[i], i, v[i*n+0], v[i*n+1], v[i*n+2], v[i*n+3] );
 *
 *    // Clean up
 *    cs_spfree( A );
 *    free( v );
 *    free( d );
 *
 *    return 0;
 * }
 * @endcode
 *
 * This example would calculate nev eigenvectors and eigenvalues of the matrix
 * A in the order of smallest magnitude first. To compile you would have to
 * do:
 *
 * @code
 * $ gcc -lceigs -larpack -lcxsparse ceigs_test.c -o ceigs_test
 * @endcode
 *
 *
 * @section Changelog
 *
 * - Version 1.1, (unreleased)
 *    - Invert the eigenvector/value order to match octave/matlab's eigs(...) function.
 *    - Support for EIGS_MODE_I_SHIFTINVERT with default driver backend.
 *    - Support for EIGS_MODE_G_SHIFTINVERT with default driver backend.
 *    - Added number of Lanczos vectors to use as a parameter.
 *    - Added driver that tries Cholesky factorization, then LU and finally QR (default).
 *    - Added driver that tries LU factorization then QR.
 *    - Added driver that tries QR factorization.
 * - Version 1.0, December 2011
 *    - Initial Revision.
 *    - Support for EIGS_MODE_I_REGULAR with default driver backend.
 *    - Support for EIGS_MODE_G_REGINVERSE with default driver backend.
 *
 *
 * @section References
 *
 * - Rich Lehoucq, Kristi Maschhoff, Danny Sorensen and Chao Yang. ARPACK. http://www.caam.rice.edu/software/ARPACK/
 *
 */


/**
 * @file ceigs.h
 * @brief The main include of the ceigs library.
 */

#define EIGS_VERSION_MAJOR   1 /**< Major version of the ceigs library. */
#define EIGS_VERSION_MINOR   0 /**< Minor version of the ceigs library. */

/**
 * @brief Ordering to return eigenvalues as.
 */
typedef enum EigsOrder_e {
   EIGS_ORDER_LA, /**< Largest algebraic eigenvalues first. */
   EIGS_ORDER_SA, /**< Smallest algebraic eigenvalues first. */
   EIGS_ORDER_LM, /**< Largest eigenvalues in magnitude first. */
   EIGS_ORDER_SM, /**< Smallest eigenvalues in magnitude first. */
   EIGS_ORDER_BE, /**< Compute nev eigenvalues, half from each end of the
                       spectrum. When nev is odd, compute one more from the
                       high end than from the low end. */
} EigsOrder_t;


/**
 * @brief Mode of operation of the ARPACK backend.
 *
 * @note Note that currently neither EIGS_MODE_G_BUCKLING nor EIGS_MODE_G_CAYLEY are
 * implemented yet.
 *
 * @note By default all the CSPARSE drivers use LU factorization.
 */
typedef enum EigsMode_e {
   /* For use in solving Av = vd */
   EIGS_MODE_I_REGULAR,     /**< For solving Av=vd in regular mode. */
   EIGS_MODE_I_SHIFTINVERT, /**< For solving Av=vd in shift-invert mode. */
   /* For use in solving Av = Mvd */
   EIGS_MODE_G_REGINVERSE,  /**< For solving Av=Mvd in regular inverse mode. */
   EIGS_MODE_G_SHIFTINVERT, /**< For solving Av=Mvd in shift-invert mode. */
   EIGS_MODE_G_BUCKLING,    /**< For solving Av=Mvd in Buckling mode. */
   EIGS_MODE_G_CAYLEY,      /**< For solving Av=Mvd in Cayley mode. */
} EigsMode_t;


/**
 * @brief Options to use.
 */
typedef struct EigsOpts_s {
   int iters;     /**< Maximum iterations during algorithm execution.
                       Default is 3000. */
   double tol;    /**< Tolerance to use. A value of 0.0 indicates to use
                       maximum machine precision.
                       Default is 0.0. */
   double sigma;  /**< Value used for the shift-invert modes to choose where to
                       calculate eigenvalues near.
                       Default is 0.0. */
   int ncv;       /**< Number of Lanczos vectors to computer (0 to autoset).
                       This value must be larger than one plus the number of
                       eigenvalues being calculated.
                       Default is 0. */
} EigsOpts_t;


/**
 * @brief Prototype for driver initialization.
 *
 * This function can generate auxiliary data structures for the backend driver
 * if necessary. The only input parameter is the size of the matrix.
 *
 * @sa EigsFreedrv_t
 */
typedef void* (*EigsInitdrv_t)( int n, const void *data_A, const void *data_M,
      const EigsOpts_t *opts );
/**
 * @brief Prototype for driver clean up.
 *
 * This function should free all data allocated by it's reciprocal driver
 * initialization function. The only input parameter is the pointer returned
 * by the driver initialization function.
 *
 * @sa EigsInitdrv_t
 */
typedef void (*EigsFreedrv_t)( void* data, const EigsOpts_t *opts );
/**
 * @brief Prototype for a dsaupd_ driver.
 *
 * These drivers are the core of the reverse communication interface used by
 * ARPACK. Their objective is to provide an implementation to manipulate the
 * data type used by the matrices with ARPACK. For exact consultation please
 * refer to the ARPACK user guide.
 */
typedef int (*EigsDsdrv_t)( int ido, int n, double *workd, const int *ipntr,
      const void *data_A, const void *data_M, void *extra );

/**
 * @brief Driver backend.
 */
typedef struct EigsDriver_s {
   EigsInitdrv_t init;  /**< Driver initialization function. Set to NULL if not needed. */
   EigsFreedrv_t free;  /**< Driver clean up function. Set to NULL if not needed. */
   EigsDsdrv_t   dsdrv; /**< Actual driver implementation. For exact specification
                             please refer to the ARPACK user guide. Set to NULL
                             to indicate this backend is not supported. */
} EigsDriver_t;

/**
 * @brief List of drivers that can be used. Any can be set to NULL to indicate unsupported.
 */
typedef struct EigsDriverGroup_s {
   EigsDriver_t driver1; /**< Driver for EIGS_MODE_I_REGULAR. */
   EigsDriver_t driver2; /**< Driver for EIGS_MODE_I_SHIFTINVERT. */
   EigsDriver_t driver3; /**< Driver for EIGS_MODE_G_REGINVERSE. */
   EigsDriver_t driver4; /**< Driver for EIGS_MODE_G_SHIFTINVERT. */
   EigsDriver_t driver5; /**< Driver for EIGS_MODE_G_BUCKLING. */
   EigsDriver_t driver6; /**< Driver for EIGS_MODE_G_CAYLEY. */
} EigsDriverGroup_t;


/**
 * @brief Driver group using Cholesky factorization.
 *
 * This driver will first attempt Cholesky factorization, then LU
 * factorization and finally QR factorization if the previous fails.
 * This driver is the default.
 *
 * @note For Cholesky factorization to work, the A matrix must be Hermitian,
 *       positive definite, and non-singular.
 */
extern const EigsDriverGroup_t eigs_drv_cholesky;
/**
 * @brief Driver group using LU factorization.
 *
 * This driver will fall back to QR factorization if LU factorization fails.
 *
 * @note For LU factorization to work, the A matrix must be non-singular.
 */
extern const EigsDriverGroup_t eigs_drv_lu; 
/**
 * @brief Driver group using QR factorization.
 *
 * @note This works with least-squares so the A matrix can be singular.
 */
extern const EigsDriverGroup_t eigs_drv_qr;


/**
 * @brief Sets the default parameters.
 *
 * This function should generally be used before manipulating the options of
 * the ARPACK for future compatibility.
 *
 *    @param[out] opts Default options for eigs.
 */
void eigs_optsDefault( EigsOpts_t *opts );


/**
 * @brief Main interface to the ARPACK eigen vector calculator.
 *
 * This either solves the problem "Ax=lx" or "Ax=Mlx".
 *
 * For the algorithm to work the following conditions must be met:
 *  - nev   <= n
 *  - nev+1 <= ncv
 *
 * Various alternative dirvers are provided, these include:
 *  - eigs_drv_lu: LU factorization driver. Works for non-singular matrices.
 *  - eigs_drv_cholesky: Cholesky factorization driver. Works for non-singular
 *                       symmetric positive definite matrices.
 *  - eigs_drv_qr: QR factorization driver. Always works, but has worse precision.
 *
 * @note When using the default drivers, matrix A and matrix B should be
 *       be compressed csparse matrices.
 *
 *    @param[in] n Order of the matrix.
 *    @param[in] nev Number of eigenvectors and eigenvalues to calculate.
 *    @param[out] lambda Eigenvalues calculated (should be of size nev).
 *    @param[out] vec Eigenvectors calculated (should be of size nev*n). If NULL
 *                    ARPACK will not calculate the eigenvectors and only compute
 *                    eigenvalues.
 *    @param[in] data_A A matrix data.
 *    @param[in] data_M M matrix data if applicable.
 *    @param[in] order Ordering to find eigenvalues in.
 *    @param[in] mode Mode of operation.
 *    @param[in] drvlist List of drivers to use or NULL to use defaults.
 *    @param[in] opts Options to use or NULL to use defaults.
 *    @return 0 on success.
 * @sa eigs_optsDefault
 */
int eigs( int n, int nev, double *lambda, double *vec,
      const void *data_A, const void *data_M,
      EigsOrder_t order, EigsMode_t mode, const EigsDriverGroup_t *drvlist,
      const EigsOpts_t *opts );

/**
 * @brief Gets the version of the library during runtime.
 *
 * The version takes the form of major.minor.
 *
 *    @param[out] major Major version of the library.
 *    @param[out] minor Minor version of the library.
 *
 * @sa EIGS_VERSION_MAJOR
 * @sa EIGS_VERSION_MINOR
 */
void eigs_version( int *major, int *minor );


#endif /* CEIG_SPARSE_H */


