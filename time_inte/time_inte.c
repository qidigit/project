/*
 * =====================================================================================
 *
 *       Filename:  time_inte.c
 *
 *    Description:  Calculate time integration with leapfrog scheme
 *
 *        Version:  1.0
 *        Created:  05/15/2015 16:26:15
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
#include "../field_management/field_management.h"
#include "time_inte.h"


int leapfrog_time_inte(spec_field *xi_pre, spec_field *delta_xi, double dt, spec_field *xi_post)
{
    int m, n;
    int cur_ind;

    for (m = 0; m < xi_post->sind[2]; m++) {
        cur_ind = xi_post->sind[0]+m*xi_post->sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {
            xi_post->spec[m][n] = xi_pre->spec[m][n] + 2.0 * dt * delta_xi->spec[m][n];
        }

        // set unused mode to be zero
        xi_post->spec[m][M_rsv+1-cur_ind] = 0;
    }

    return 0;
}
