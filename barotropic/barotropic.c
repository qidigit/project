/*
 * =====================================================================================
 *
 *       Filename:  barotropic.c
 *
 *    Description:  model for barotropic vorticity flow
 *
 *        Version:  1.0
 *        Created:  05/15/2015 17:32:37
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
#include "barotropic.h"

// these variables cannot be used before initialization
trans_group u_gp, v_gp, curl_gp;
spec_field delta_curl, curl_pre, curl_post;
spec_field psi_s, pot_s;
grid_field freq;
legendre trans;

int initialize_barotropic()
{
    int m, n, i, j;
    int cur_ind;
    int ind[3], sind[3];

    int nproc, rproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);


    // set up the splitted modes in each proc
    ind[0] = rproc*(NY/nproc); ind[1] = 1; ind[2] = NY / nproc;
    sind[0] = rproc; sind[1] = nproc; sind[2] = (M+1) / nproc;
    printf("\nproc rank %d initialized...\n", rproc);
    printf("ind  = [%d %d %d],\t",  ind[0],  ind[1],  ind[2]);
    printf("sind = [%d %d %d]\n", sind[0], sind[1], sind[2]);

    MPI_Barrier(MPI_COMM_WORLD);

    printf("initializing variables...\n");

    // initialize groups
    initialize_trans_group(ind, sind, &u_gp);
    initialize_trans_group(ind, sind, &v_gp);
    initialize_trans_group(ind, sind, &curl_gp);

    // initialize spec field
    initialize_spec_field(sind, &delta_curl);
    initialize_spec_field(sind, &curl_pre);
    initialize_spec_field(sind, &curl_post);

    initialize_spec_field(sind, &psi_s);
    initialize_spec_field(sind, &pot_s);

    // set the unused potential flow to zero
    for (m = 0; m < pot_s.sind[2]; m++) {
        cur_ind = pot_s.sind[0]+m*pot_s.sind[1];
        for (n = 0; n < M+2-cur_ind; n++) {
            pot_s.spec[m][n] = 0.0;
            psi_s.spec[m][n] = 0.0;
            curl_post.spec[m][n] = 0.0;
        }
    }

    // initialize transformation operator
    initialize_transform(ind, sind, &trans);

    // initialize rotation
    initialize_grid_field(ind, &freq);
    double theta;
    for (j = 0; j < ind[2]; j++) {
        cur_ind = ind[0]+j*ind[1];
        if (cur_ind < ((NY+1)/2)) {
            theta = asin(trans.legendre_points[cur_ind]);
        } else {
            theta = asin(-trans.legendre_points[cur_ind-(NY/2)]);
        }
        for (i = 0; i < (NY*2); i++) {
            freq.grid[j*(NY*2)+i] = 2.0 * OMEGA * sin(theta);
        }
    }

    return 0;
}

int finalize_barotropic()
{
    finalize_trans_group(&u_gp);
    finalize_trans_group(&v_gp);
    finalize_trans_group(&curl_gp);

    finalize_spec_field(&delta_curl);
    finalize_spec_field(&curl_pre);
    finalize_spec_field(&curl_post);
    finalize_spec_field(&psi_s);
    finalize_spec_field(&pot_s);

    finalize_transform(&trans);
    finalize_grid_field(&freq);

    return 0;
}

// set initial values for u, v, curl, curl_pre
int set_initial_value(int num)
{
    int i, j, m, n;
    int cur_ind;

    double theta;
    int mode = 4;
    double theta0 = 45, theta_w = 15;
    double drate;

    printf("Setting up initial values!\n");

    switch(num)
    {
        case 1 :

            for (j = 0; j < u_gp.ind[2]; j++) {
                cur_ind = u_gp.ind[0]+j*u_gp.ind[1];
                if (cur_ind < ((NY+1)/2)) {
                    theta = asin(trans.legendre_points[cur_ind]);
                } else {
                    theta = asin(-trans.legendre_points[cur_ind-(NY/2)]);
                }

                for (i = 0; i < (2*NY); i++) {
                    u_gp.grid[j*(2*NY)+i] = 25.0 * cos(theta) - 30.0 * pow(cos(theta), 3) +
                                             300.0 * pow(sin(theta), 2) * pow(cos(theta), 6);
                    v_gp.grid[j*(2*NY)+i] = 0.0;

                    curl_gp.grid[j*(2*NY)+i] = 25.0 * sin(theta) - 90.0 * pow(cos(theta), 2) * sin(theta) -
                    600.0 * sin(theta) * pow(cos(theta), 7) + 1800.0 * pow(sin(theta), 3) * pow(cos(theta), 5) +
                    sin(theta) / cos(theta) * u_gp.grid[j*(2*NY)+i];
                    curl_gp.grid[j*(2*NY)+i] /= RADIUS;

                    // add perturbation
                    drate = (theta*180.0/PI - theta0) / theta_w;
                    curl_gp.grid[j*(2*NY)+i] += 8.0e-5 / 2 * cos(theta) * exp(-pow(drate, 2)) *
                                                 cos((double) mode * PI * i / (double) NY);
                }
            }

            transform_group_forward(&curl_gp, &trans);
            for (m = 0; m < curl_pre.sind[2]; m++) {
                cur_ind = curl_pre.sind[0]+m*curl_pre.sind[1];
                for (n = 0; n < M+2-cur_ind; n++) {
                    curl_pre.spec[m][n] = curl_gp.spec[m][n];
                }
            }

            break;
        default :
            printf("Invalid initialization choice!\n");
            break;
    }

    return 0;
}


// this function takes the initial values of u, v, and curl, update the tendency delta_curl
// u_gp->grid, v_gp->grid, curl_gp->grid & spec, curl_pre must be initialized;
// (u_gp, v_gp) values will be destroyed during the calculation.
int barotropic_integration()
{
    int i, j, m, n;
    int cur_ind;

    double impli;
    double sigma;

    // set fields to get divergence
    for (j = 0; j < u_gp.ind[2]; j++) {
        cur_ind = u_gp.ind[0]+j*u_gp.ind[1];

        for (i = 0; i < (2*NY); i++) {
            u_gp.grid[j*(2*NY)+i] = (-1.0) * u_gp.grid[j*(2*NY)+i] * 
                                    (curl_gp.grid[j*(2*NY)+i] + freq.grid[j*(2*NY)+i]);
            v_gp.grid[j*(2*NY)+i] = (-1.0) * v_gp.grid[j*(2*NY)+i] *
                                    (curl_gp.grid[j*(2*NY)+i] + freq.grid[j*(2*NY)+i]);
        }
    }

        // for test purpose
//        int start = 0;
//        printf("uv grid\n");
//        for (i = 0; i < NY; i+=1) {
//            printf("%12g\t%12g\n", u_gp.grid[i*2*NY+start], v_gp.grid[i*2*NY+start]);
//        }


    // get spectral representation of the divergence
    div_spec_from_uv_grid(&u_gp, &v_gp, &trans, &delta_curl);


    // add hyper-diffusivity to the tendency
    spectral_damping(&curl_pre, &delta_curl);

    for (m = 0; m < delta_curl.sind[2]; m++) {
        cur_ind = delta_curl.sind[0]+m*delta_curl.sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {

            sigma = 1;
            for (i = 0; i < damping_order; i++) {
                sigma = sigma * (double) (n+cur_ind) * (double) (n+1+cur_ind) / (RADIUS*RADIUS);
            }
//            sigma = (double) (n+cur_ind) * (double) (n+1+cur_ind) / (RADIUS*RADIUS);
//            sigma = pow(sigma, damping_order);
            
            impli = 1.0 + damping_coeff * 2 * delta_t * sigma;
            delta_curl.spec[m][n] = delta_curl.spec[m][n] / impli;

            // for test purpose
//            if (m == 0 && n == 0) {
//                printf("impli = %f\n", impli);
//            }
        }
    }

        // for test purpose
//        int rproc;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rproc);
//        if (rproc == 1) {
//        int i;
//        printf("\ntime tendency\n");
//        for (i = 0; i < M+2; i++) {
//            printf("%12g+%12gi\n", creal(delta_curl.spec[0][i]), cimag(delta_curl.spec[0][i]));
//        }
//        }

    return 0;
}


// this function sets up the values for the next updating circle
// Input: curl_post
// Output: u_gp, v_gp, curl_gp, curl_pre
int barotropic_update()
{
    int m, n;
    int cur_ind;

    for (m = 0; m < psi_s.sind[2]; m++) {
        cur_ind = psi_s.sind[0]+m*psi_s.sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {
            // set up the values for psi_s -- divergence
            psi_s.spec[m][n] = curl_post.spec[m][n];

            // initialize new vorticity from Robert filter
            curl_pre.spec[m][n] = (1.0-2.0*robert) * curl_gp.spec[m][n] + 
                             robert * (curl_post.spec[m][n] + curl_pre.spec[m][n]);

            // update the new spectral vorticity at new time
            curl_gp.spec[m][n] = curl_post.spec[m][n]; 
        }
    }

    // update vorticity at grid points at new time
    transform_group_backward(&curl_gp, &trans);

    // update velocity field from div & curl at new time
    uv_grid_from_div_curl(&pot_s, &psi_s, &trans, &u_gp, &v_gp);

    return 0;
}
