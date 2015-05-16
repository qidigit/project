/*
 * =====================================================================================
 *
 *       Filename:  exec_baro.c
 *
 *    Description:  This is the main function to execute barotropic model
 *
 *        Version:  1.0
 *        Created:  05/15/2015 23:55:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include <mpi.h>
#include "../model_constants/model_constants.h"
#include "../transform_s2/gauss_legendre.h"
#include "../transform_s2/legendre_polynomial.h"
#include "../field_management/field_management.h"
#include "../transform_s2/transform_s2.h"
#include "../spectral_diff/spectral_diff.h"
#include "../time_inte/time_inte.h"
#include "../barotropic/barotropic.h"
#include "../time_inte/time_inte.h"

int main(int argc, const char *argv[])
{
    // initialize mpi
    int nproc, rproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);

    // set up model resolution, need to be power of 2, 3, 5
    NY = 256;
    initial_resolution();
    if (rproc == 0) {
        printf("Model resolution is set as: NY = %d, M_rsv = %d, M = %d.\n", rproc, NY, M_rsv, M);
        printf("Calculation up to %d days.\n", (int) (TT/24/60/60));
    }

    // total number of integration step
    delta_t = 1800.0;
    int nstep = (int) (TT / delta_t);

    // initialize model variables
    initialize_barotropic();

    // set up initial values
    int num = 1;
    set_initial_value();

    int step_count = 0;
    while (step_count < nstep) {
        // one step integration
        barotropic_integration();
        leapfrog_time_inte(curl_pre, delta_curl, delta_t, curl_post);
        barotropic_update();

        // diagnostic step
        
        // increment
        step_count++;
    }

    finalize_barotropic();
    MPI_Finalize();
    return 0;
}
