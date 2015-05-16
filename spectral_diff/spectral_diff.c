/*
 * =====================================================================================
 *
 *       Filename:  spectral_diff.c
 *
 *    Description:  This part of codes calculates the divergence and curl of vector field (u, v) and also the inverse transform for (u, v) field. It will also calculate forward & backward Laplacian operator.
 *
 *        Version:  1.0
 *        Created:  05/14/2015 13:41:01
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
#include "spectral_diff.h"

// This function calculates the divergence of vector field (u, v) in spectral domain
// input: u_gp, v_gp, trans, NOTE the values in u_gp, v_gp will be destroyed during calculation!
// output: div_s
int div_spec_from_uv_grid(trans_group *u_gp, trans_group *v_gp, legendre *trans, spec_field *div_s)
{
    int i, j, m, n;
    double theta;
    double eps_mn, eps_mn1;

    int cur_ind;
    // devide the flow field by cos(theta)
    for (i = 0; i < u_gp->ind[2]; i++) {
        // calculate theta on each layer
        cur_ind = u_gp->ind[0] + i*u_gp->ind[1];
        if (cur_ind < (NY+1)/2) {
            theta = asin(trans->legendre_points[cur_ind]);
        } else {
            theta = asin(-trans->legendre_points[cur_ind-(NY/2)]);
        }

        for (j = 0; j < (NY*2); j++) {
            u_gp->grid[i*(2*NY)+j] /= cos(theta);  
            v_gp->grid[i*(2*NY)+j] /= cos(theta);
        }
    }

    // get spectral representation of (u/cos^2(theta), v/cos^(theta))
    transform_group_forward(u_gp, trans);
    transform_group_forward(v_gp, trans);

    // calculate divergence
    for (m = 0; m < div_s->sind[2]; m++) {
        cur_ind = div_s->sind[0]+m*div_s->sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {
            eps_mn1 = sqrt((double) ((n+1+cur_ind)*(n+1+cur_ind) - cur_ind * cur_ind) /
                           (double) (4 * (n+1+cur_ind)*(n+1+cur_ind) - 1));
            div_s->spec[m][n] = (I * (double) cur_ind) * u_gp->spec[m][n] +
                                (double) (n + cur_ind) * eps_mn1 * v_gp->spec[m][n+1];
            
            if ((n+cur_ind) > 0) {
                eps_mn  = sqrt((double) ((n+cur_ind)*(n+cur_ind) - cur_ind * cur_ind) /
                           (double) (4 * (n+cur_ind)*(n+cur_ind) - 1));
                div_s->spec[m][n] -= (double) (n+1+cur_ind) * eps_mn * v_gp->spec[m][n-1];
            }

            div_s->spec[m][n] /= RADIUS;
        }
        // set the additional unused mode to be zero
        div_s->spec[m][M_rsv+1-cur_ind] = 0;
    }

    return 0;
}

// This function calculate the vector field (u, v) in physical domain from its divergence and curl
// input: div_s, curl_s, trans, NOTE: the values in div_s, curl_s will be destroyed during calculation
// output: u_gp, v_gp -- flow velocity fields
int uv_grid_from_div_curl(spec_field *div_s, spec_field *curl_s, legendre *trans, trans_group *u_gp, trans_group *v_gp)
{
    int i, j, m, n;
    double theta;
    double eps_mn, eps_mn1;

    int cur_ind;

    // transform div, curl back to psi & xi
    laplace_backward(div_s, div_s);
    laplace_backward(curl_s, curl_s);

    // calculate spectral U/a, V/a
    for (m = 0; m < div_s->sind[2]; m++) {
        cur_ind = div_s->sind[0]+m*div_s->sind[1];
        for (n = 0; n < M_rsv+2-cur_ind; n++) {
            u_gp->spec[m][n] = I * (double) cur_ind * div_s->spec[m][n];
            v_gp->spec[m][n] = I * (double) cur_ind * curl_s->spec[m][n];

            if ((n+cur_ind) > 0) {
                eps_mn  = sqrt((double) ((n+cur_ind)*(n+cur_ind) - cur_ind * cur_ind) /
                           (double) (4 * (n+cur_ind)*(n+cur_ind) - 1));

                u_gp->spec[m][n] += (double) (n-1+cur_ind) * eps_mn * curl_s->spec[m][n-1];
                v_gp->spec[m][n] -= (double) (n-1+cur_ind) * eps_mn * div_s->spec[m][n-1];
            }
            if ((n+cur_ind) < (M_rsv+1)) {
                eps_mn1 = sqrt((double) ((n+1+cur_ind)*(n+1+cur_ind) - cur_ind * cur_ind) /
                           (double) (4 * (n+1+cur_ind)*(n+1+cur_ind) - 1));

                u_gp->spec[m][n] -= (double) (n+2+cur_ind) * eps_mn1 * curl_s->spec[m][n+1]; 
                v_gp->spec[m][n] += (double) (n+2+cur_ind) * eps_mn1 * div_s->spec[m][n+1];
            }
        }
    }

    // transform fields back to physical domain
    transform_group_backward(u_gp, trans);
    transform_group_backward(v_gp, trans);

    // normalize the velocity field
    for (i = 0; i < u_gp->ind[2]; i++) {
        // calculate theta on each layer
        cur_ind = u_gp->ind[0] + i*u_gp->ind[1];
        if (cur_ind < (NY+1)/2) {
            theta = asin(trans->legendre_points[cur_ind]);
        } else {
            theta = asin(-trans->legendre_points[cur_ind-(NY/2)]);
        }

        for (j = 0; j < (2*NY); j++) {
            u_gp->grid[i*(2*NY)+j] *= RADIUS/cos(theta);
            v_gp->grid[i*(2*NY)+j] *= RADIUS/cos(theta);
        }
    }

    return 0;
}


// Note these functions can be only applied to (stream func, potential func) or (div, curl)
int laplace_forward(spec_field *psi_in, spec_field *curl_out)
{
    int m, n;
    int cur_ind;

    for (m = 0; m < psi_in->sind[2]; m++) {
        cur_ind = psi_in->sind[0]+m*psi_in->sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {
            curl_out->spec[m][n] = (-1.0) * (double) (n+cur_ind) * (double) (n+1+cur_ind) * psi_in->spec[m][n];
        }

        // set the unused mode zero
        curl_out->spec[m][M_rsv+1-cur_ind] = 0;
    }
    return 0;
}
int laplace_backward(spec_field *curl_in, spec_field *psi_out)
{
    int m, n;
    int cur_ind;

    for (m = 0; m < curl_in->sind[2]; m++) {
        cur_ind = curl_in->sind[0]+m*curl_in->sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {
            psi_out->spec[m][n] = (-1.0) * curl_in->spec[m][n] / ((double) (n+cur_ind) * (double) (n+1+cur_ind));
        }

        // set the unused mode zero
        psi_out->spec[m][M_rsv+1-cur_ind] = 0;
    }
    return 0;
}

int laplace_order(spec_field *psi_in, spec_field *psi_out, int lorder)
{
    int m, n, i;
    int cur_ind;
    double sigma;

    for (m = 0; m < psi_in->sind[2]; m++) {
        cur_ind = psi_in->sind[0]+m*psi_in->sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {
            sigma = 1;
            for (i = 0; i < lorder; i++) {
                sigma = (-1.0) * sigma * (double) (n+cur_ind) * (double) (n+1+cur_ind) / (RADIUS*RADIUS);
            }
//            sigma = (double) (n+cur_ind) * (double) (n+1+cur_ind) / (RADIUS*RADIUS);
//            sigma = pow(sigma, damping_order);

            psi_out->spec[m][n] = sigma * psi_in->spec[m][n];
        }

        // set the unused mode zero
        psi_out->spec[m][M_rsv+1-cur_ind] = 0;
    }
    
    return 0;
}

// calculate spectral damping, additional hyperdiffusion in psi_in is added to increament incr_out
int spectral_damping(spec_field *psi_in, spec_field *incr_out)
{
    int m, n, i;
    int cur_ind;
    double sigma;

    for (m = 0; m < psi_in->sind[2]; m++) {
        cur_ind = psi_in->sind[0]+m*psi_in->sind[1];
        for (n = 0; n < M_rsv+1-cur_ind; n++) {

            sigma = 1;
            for (i = 0; i < damping_order; i++) {
                sigma = sigma * (double) (n+cur_ind) * (double) (n+1+cur_ind) / (RADIUS*RADIUS);
            }
//            sigma = (double) (n+cur_ind) * (double) (n+1+cur_ind) / (RADIUS*RADIUS);
//            sigma = pow(sigma, damping_order);
            incr_out->spec[m][n] = incr_out->spec[m][n] - damping_coeff * sigma * psi_in->spec[m][n];
        }

        // for safety, set the unused mode zero
        incr_out->spec[m][M_rsv+1-cur_ind] = 0;
    }
    return 0;
}


int multi_mu_to_spec(spec_field *psi_in, spec_field *multi_psi)
{
    int m, n;
    double eps_mn, eps_mn1;

    int cur_ind;

    // calculate spectral modes
    for (m = 0; m < psi_in->sind[2]; m++) {
        cur_ind = psi_in->sind[0]+m*psi_in->sind[1];
        for (n = 0; n < M_rsv+2-cur_ind; n++) {
            multi_psi->spec[m][n] = 0;

            if ((n+cur_ind) > 0) {
                eps_mn  = sqrt((double) ((n+cur_ind)*(n+cur_ind) - cur_ind * cur_ind) /
                           (double) (4 * (n+cur_ind)*(n+cur_ind) - 1));

                multi_psi->spec[m][n] += eps_mn * psi_in->spec[m][n-1];
            }
            if ((n+cur_ind) < (M_rsv+1)) {
                eps_mn1 = sqrt((double) ((n+1+cur_ind)*(n+1+cur_ind) - cur_ind * cur_ind) /
                           (double) (4 * (n+1+cur_ind)*(n+1+cur_ind) - 1));

                multi_psi->spec[m][n] = eps_mn1 * psi_in->spec[m][n+1]; 
            }
        }
    }

    return 0;
}
