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
#include "../diagnostics/diagnostics.h"


int main(int argc, char **argv)
{
    // initialize mpi
    int nproc, rproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);

    // set up model resolution, need to be power of 2, 3, 5
    NY = 256;
    TT = 50;
    initial_resolution();
    if (rproc == 0) {
        printf("Model resolution is set as: NY = %d, M_rsv = %d, M = %d.\n", NY, M_rsv, M);
        printf("Calculation up to %d days.\n", (int) (TT));
    }

    // total number of integration step & diagnostic time step
    delta_t  = 600.0;
    delta_at = 120*600.0;
    int nstep = (int) round(TT * 24.0 * 60.0 * 60.0 / delta_t);
    int dstep = (int) round(delta_at / delta_t);
//    printf("dstep = %d\n", dstep);

    // initialize model variables
    initialize_barotropic();

    // set up initial values
    int num = 1;
    set_initial_value(num);


    // for test purpose
//        int start = 4;
//        int i, m, n, cur_ind;
//        if (rproc == 0) {
//        printf("curl pre\n");
//        for (i = 0; i < M+2; i++) {
//            printf("%12g+%12gi\n", creal(curl_pre.spec[start][i]), cimag(curl_pre.spec[start][i]));
//        }
//       }


    // diagnostic initialize
    char f_status = 'w';
    diag_initial(f_status);

    int step_count = 0;
    while (step_count < nstep) {
        // one step integration
        barotropic_integration();

        // for test purpose
//        if (rproc == 0) {
//        printf("time tendency\n");
//        for (i = 0; i < M+2; i++) {
//            printf("%12g+%12gi\n", creal(delta_curl.spec[start][i]), cimag(delta_curl.spec[start][i]));
//        }
//        for (m = 0; m < delta_curl.sind[2]; m++) {
//            cur_ind = delta_curl.sind[0]+m*delta_curl.sind[1];
//            for (n = 0; n < M+2-cur_ind; n++) {
//                if (isnan(creal(delta_curl.spec[m][n])) || isnan(cimag(delta_curl.spec[m][n])) ) {
//                    printf("delta_curl becomes NaN at m=%d, n=%d\n", m, n);
//                }        
//            }
//        }
//        }


        leapfrog_time_inte(&curl_pre, &delta_curl, delta_t, &curl_post);

        // for test purpose
//        if (rproc == 0) {
//        printf("curl post\n");
//        for (i = 0; i < M+2; i++) {
//            printf("%12g+%12gi\n", creal(curl_post.spec[start][i]), cimag(curl_post.spec[start][i]));
//        }
//        }

        barotropic_update();

//        for (m = 0; m < delta_curl.sind[2]; m++) {
//            cur_ind = delta_curl.sind[0]+m*delta_curl.sind[1];
//            for (n = 0; n < M+2-cur_ind; n++) {
//                if (delta_curl.spec[m][n] != delta_curl.spec[m][n] ) {
//                    printf("delta_curl becomes NaN at m=%d, n=%d\n", m, n);
//                }        
//            }
//        }

        // for test purpose
//        if (rproc == 0) {
//        printf("curl update\n");
//        for (i = 0; i < M+2; i++) {
//            printf("%12g+%12gi\n", creal(curl_gp.spec[start][i]), cimag(curl_gp.spec[start][i]));
//        }
//        printf("curl grid\n");
//        for (i = 0; i < NY/nproc; i+=1) {
//            printf("%12g\n", curl_gp.grid[i*2*NY+start]);
//        }
//        }

 
        // count increment
        step_count++;

        // diagnostic step
        if ((step_count % dstep) == 0) {
            diag_writeData();
        }
    }

    finalize_barotropic();
    MPI_Finalize();
    return 0;
}
