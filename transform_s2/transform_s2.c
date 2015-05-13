/*
 * =====================================================================================
 *
 *       Filename:  transform_s2.c
 *
 *    Description:  This file calculates the transform between physical domain and spectral domain under spherical harmonics
 *
 *        Version:  1.0
 *        Created:  04/21/2015 21:54:37
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
#include <fftw3.h>
#include <omp.h>
#include <mpi.h>
#include "gauss_legendre.h"
#include "legendre_polynomial.h"
#include "../field_management/field_management.h"
#include "transform_s2.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

int initialize_transform(int ind[], int sind[], legendre *trans)
{
    double *weight;
    double x_pos[NY];
    double *pmn_vec;
    int m = (NY+1)>>1;
    int i, j, k;
    int ind_count, sind_count;

    // allocate memory
    weight = (double*) malloc(sizeof(double) * NY);
    trans->legendre_points = (double*) malloc(sizeof(double) * m);
    // by default we use triangular truncation
    trans->leg_trans_forward = (double***) malloc(sizeof(double**) * (M+1));
    for (i = 0; i < M+1; i++) {
        trans->leg_trans_forward[i] = (double**) malloc(sizeof(double*) * (M+2-i));
        for (j = 0; j < M+2-i; j++)
            trans->leg_trans_forward[i][j] = (double*) malloc(sizeof(double) * ind[2]);
    }
    trans->leg_trans_backward = (double***) malloc(sizeof(double**) * NY);
    for (i = 0; i < NY; i++) {
        trans->leg_trans_backward[i] = (double**) malloc(sizeof(double*) * sind[2]);
        for (j = 0; j < sind[2]; j++)
            trans->leg_trans_backward[i][j] = (double*) malloc(sizeof(double) * (M+2-sind[0]-j*sind[1]));
    }

    // initialize index
    trans->ind[0] = ind[0];
    trans->ind[1] = ind[1];
    trans->ind[2] = ind[2];
    trans->sind[0] = sind[0];
    trans->sind[1] = sind[1];
    trans->sind[2] = sind[2];


    // calculate the abscissas and weights for Gauss-Legendre quadrature
    gauss_legendre_tbl(NY, trans->legendre_points, weight, EPS);

    // setup the Legendre transformation matrix
    /* x_pos stores the gaussian quadrature points in the order of (0,x{1}<x{2}<...<x{(NY+1)/2},-x{1}>-x{2}>...>-x{(NY+1)/2}), 0 only appears when NY is odd.
     * */
    for (i = 0; i < m; i++) {
        x_pos[i] = trans->legendre_points[i];
        if (NY%2 == 0) {
            x_pos[m+i] = -trans->legendre_points[i];
            weight[m+i] = weight[i];
        } else if (i > 0) {
            x_pos[m+i-1] = -trans->legendre_points[i];
            weight[m+i-1] = weight[i];
        }
    }
    // check for initialization
    /*
    for (i = 0; i < NY; i++) {
        printf("x[%d]=%f,\tw[%d]=%f\n", i, x_pos[i], i, weight[i]);
    }
    */

    sind_count = 0;
    for (j = 0; j < M+1; j++) {
        pmn_vec = pmn_polynomial_value(NY, M+1, j, x_pos);

        ind_count = 0;
        for (i = 0; i < NY; i++) {
            if ((ind_count < ind[2]) && (i == ind[0]+ind_count*ind[1])) {
                for (k = j; k < M+2; k++) {
                    trans->leg_trans_forward[j][k-j][ind_count] = pmn_vec[i+k*NY]*weight[i];
                }
                ind_count++; 
            }
        }

        if ((sind_count < sind[2]) && (j == sind[0]+sind_count*sind[1])) {
            for (i = 0; i < NY; i++) {
                for (k = j; k < M+2; k++) {
                    trans->leg_trans_backward[i][sind_count][k-j] = pmn_vec[i+k*NY];
                }
            }
            sind_count++;
        }

//        printf("ind_count =%d, j=%d\n", ind_count, j);
//        printf("sind_count=%d, j=%d\n", sind_count, j);
        free(pmn_vec);
    }

    free(weight);
    return 0;
}

int finalize_transform(legendre *trans)
{
    int i, j;

    free(trans->legendre_points);
    for (i = 0; i < M+1; i++) {
        for (j = 0; j < M+2-i; j++)
            free(trans->leg_trans_forward[i][j]);
        free(trans->leg_trans_forward[i]);
    }
    free(trans->leg_trans_forward);

    for (i = 0; i < NY; i++) {
        for (j = 0; j < trans->sind[2]; j++)
            free(trans->leg_trans_backward[i][j]);
        free(trans->leg_trans_backward[i]);
    }
    free(trans->leg_trans_backward);

    return 0;
}

int transform_g2s(sphere_grid *xi_g, sphere_four *xi_f, sphere_spec *xi_s, legendre *trans)
{
    int i, m, n;
    int cur_p, cur_mul, cur_ind;
    int bsize = xi_s->sind[2]*(M_rsv+2);

    int rproc, nproc; // define the mpi processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // transform from grids to Fourier domain, with normalization 2*NY
    fftw_execute(xi_g->p_forward);
    // normalization
    for (i = 0; i < (NY+1)*xi_f->ind[2]; i++) xi_f->u_four[i] /= (2*NY);

    // calculate spherical coeff. through Legendre transform
    // piece of transform inside the proc
    double *leg_sum_r, *leg_sum_i;
    leg_sum_r = (double*) calloc((M+1)*(M_rsv+2), sizeof(double*));
    leg_sum_i = (double*) calloc((M+1)*(M_rsv+2), sizeof(double*));

    for (m = 0; m < M+1; m++) {
        // reorder the zonal wavenumber in order to save them in the same order in each proc
        // further improvement for efficiency needs to store the modes better
        cur_p = m % xi_s->sind[1];
        cur_mul = (m-cur_p) / xi_s->sind[1];
        cur_ind = cur_p*xi_s->sind[2] + cur_mul;
        for (n = 0; n < M_rsv+2-m; n++) {

            for (i = 0; i < xi_f->ind[2]; i++) {
                leg_sum_r[cur_ind*(M_rsv+2)+n] += trans->leg_trans_forward[m][n][i] * creal(xi_f->u_four[i*(NY+1)+m]);
                leg_sum_i[cur_ind*(M_rsv+2)+n] += trans->leg_trans_forward[m][n][i] * cimag(xi_f->u_four[i*(NY+1)+m]);
            }
        }
    }

    // reduce with other processors.
    int root = 0;
//    MPI_Allreduce(MPI_IN_PLACE, leg_sum, (M_rsv+1)*(M_rsv+2), MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    if (rproc == root) {
        MPI_Reduce(MPI_IN_PLACE, leg_sum_r, (M+1)*(M_rsv+2), MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, leg_sum_i, (M+1)*(M_rsv+2), MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(leg_sum_r,    leg_sum_r, (M+1)*(M_rsv+2), MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(leg_sum_i,    leg_sum_i, (M+1)*(M_rsv+2), MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    }

    double *temp_order_r, *temp_order_i; // temp to store the reduced data
    temp_order_r = (double*) calloc(xi_s->sind[2]*(M_rsv+2), sizeof(double));
    temp_order_i = (double*) calloc(xi_s->sind[2]*(M_rsv+2), sizeof(double));

    // Scatter data to each proc
    MPI_Scatter(leg_sum_r, bsize, MPI_DOUBLE, temp_order_r, bsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Scatter(leg_sum_i, bsize, MPI_DOUBLE, temp_order_i, bsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    // reallocate the array
    for (m = 0; m < xi_s->sind[2]; m++) {
        cur_ind = xi_s->sind[0]+m*xi_s->sind[1];
        for (i = 0; i < M_rsv+2-cur_ind; i++) {
            xi_s->u_spec[m][i] = temp_order_r[m*(M_rsv+2)+i] + 
                                 I * temp_order_i[m*(M_rsv+2)+i];
        }
//        for (i = MAX(M_rsv+2-cur_ind, 0); i < M+2-cur_ind; i++) {
//            xi_s->u_spec[m][i] = 0;
//        }
    }

    free(temp_order_r);
    free(temp_order_i);
    free(leg_sum_r);
    free(leg_sum_i);
    return 0;
}

int transform_s2g(sphere_spec *xi_s, sphere_four *xi_f, sphere_grid *xi_g, legendre *trans)
{
    int i, mi, n;
    int cur_ind;

    int bsize = xi_f->ind[2]*xi_s->sind[2];
    int rproc, nproc; // define the mpi processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // calculate the Fourier coeff. from available spherical modes sind[0]:sind[1]:sind[0]+sind[2]*sind[1]
    double *leg_dec_r, *leg_dec_i;
    leg_dec_r = (double*) calloc(NY * (xi_s->sind[2]), sizeof(double));
    leg_dec_i = (double*) calloc(NY * (xi_s->sind[2]), sizeof(double));

    // For safety: set the useless mode to 0 to avoid error; this should be guaranteed before enter
//    for (mi = 0; mi < xi_s->sind[2]; mi++) {
//        cur_ind = xi_s->sind[0]+mi*xi_s->sind[1];
//        for (i = MAX(M_rsv+2-cur_ind, 0); i < M+2-cur_ind; i++) 
//            xi_s->u_spec[mi][i] = 0;
//    }

    // for test only
//    for (i = 0; i < M+2; i++) {
//        printf("%12g+%12gi\n", creal(xi_s->u_spec[0][i]), cimag(xi_s->u_spec[0][i]));
//    }
//    printf("u_spec test\n");
//    for (i = 0; i < M+2; i++) {
//        printf("%12g %12g\n", trans->leg_trans_backward[0][0][i], trans->leg_trans_backward[NY-1][0][i]);
//    }
//    printf("legendre test\n");

    for (i = 0; i < NY; i++) {
        for (mi = 0; mi < xi_s->sind[2]; mi++) {
            cur_ind = xi_s->sind[0]+mi*xi_s->sind[1];
            
            for (n = 0; n < M_rsv+2-cur_ind; n++) {
                leg_dec_r[i*xi_s->sind[2]+mi] += trans->leg_trans_backward[i][mi][n] * creal(xi_s->u_spec[mi][n]);
                leg_dec_i[i*xi_s->sind[2]+mi] += trans->leg_trans_backward[i][mi][n] * cimag(xi_s->u_spec[mi][n]);
            }
        }
    }

    // for test only
//    for (i = 0; i < NY; i++) {
//        printf("%12g+%12gi\n", leg_dec_r[i*(M+1)+0],    leg_dec_i[i*(M+1)+0]);
//    }
//    printf("next test\n");

    // communicate with other processors; note the constraint nproc*sind[2] = M+1
    double *temp_order_r, *temp_order_i; // temp to store the disordered collection
    temp_order_r = (double*) malloc(sizeof(double) * (xi_f->ind[2]*(M+1)));
    temp_order_i = (double*) malloc(sizeof(double) * (xi_f->ind[2]*(M+1)));
    MPI_Alltoall(leg_dec_r, bsize, MPI_DOUBLE, temp_order_r, bsize, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoall(leg_dec_i, bsize, MPI_DOUBLE, temp_order_i, bsize, MPI_DOUBLE, MPI_COMM_WORLD);
    // reorder the received data, the orders in the initial setting-up MATTERS here!
    /* write in the order of output (write) file
    for (i = 0; i < xi_f->ind[2]; i++) {
        for (mi = 0; mi < xi_s->sind[2]; mi++) {
            for (n = 0; n < nproc; n++) {
                xi_f->u_four[i*(NY+1) + mi*nproc + n] = temp_order[n*bsize + i*xi_s->sind[2] + mi];
            }
        }
    } */

    // for test only
//    for (i = 0; i < NY; i++) {
//        printf("%12g+%12gi  %12g+%12gi\n", temp_order_r[i*(M+1)+0], temp_order_i[i*(M+1)+0], 
//                                           leg_dec_r[i*(M+1)+0],    leg_dec_i[i*(M+1)+0]);
//    }
//    printf("next test\n");

    /* read in the order of input (read) file */
    for (n = 0; n < nproc; n++) {
        for (i = 0; i < xi_f->ind[2]; i++) {
            for (mi = 0; mi < xi_s->sind[2]; mi++) {
                xi_f->u_four[i*(NY+1) + mi*nproc + n] = temp_order_r[n*bsize + i*xi_s->sind[2] + mi] + 
                                                        I * temp_order_i[n*bsize + i*xi_s->sind[2] + mi];
            }
        }
    }

    // de-aliasing
    // or it may need to start at a smaller value M_res to show the real number of resolved mode
    for (i = 0; i < xi_f->ind[2]; i++) {
        for (mi = M_rsv+1; mi < NY+1; mi++) {
            xi_f->u_four[i*(NY+1) + mi] = 0;
        }
    }

    // for test only
//    for (i = 0; i < NY; i++) {
//        printf("%12g+%12gi\n", creal(xi_f->u_four[i*(NY+1)+0]), cimag(xi_f->u_four[i*(NY+1)+0]));
//    }

    // tansform from Fourier domain to grids, no normalization
    // NOTE original xi_f->u_four is destroyed after the inverse transform
    fftw_execute(xi_f->p_backward);
//    fftw_execute_dft_c2r(xi_f->p_backward, xi_f->u_four, xi_g->u_grid);

    free(temp_order_r);
    free(temp_order_i);
    free(leg_dec_r);
    free(leg_dec_i);
    return 0;
}


