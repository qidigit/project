/*
 * =====================================================================================
 *
 *       Filename:  transform_s2.c
 *
 *    Description:  This file calculates the transform between physical domain and spectral domain under spherical harmonics
 *    NOTE: THIS IS THE PREVIOUS TEMP FILE FOR TEST, NO ACTUAL USE!
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
#include "gauss_legendre.h"
#include "legendre_polynomial.h"
#include "transform_s2.h"

int initialize_transform(int ind[], int sind[], legendre *trans)
{
    double *weight;
    double x_pos[NY];
    double *pmn_vec;
    int m = (NY+1)>>1;
    int i, j, k;
    int ind_count, sind_count;

    // allocate memory
    weight = (double*) malloc(sizeof(double) * m);
    trans->legendre_points = (double*) malloc(sizeof(double) * m);
    // by default we use triangular truncation
    trans->leg_trans_forward = (double***) malloc(sizeof(double**) * ind[2]);
    for (i = 0; i < ind[2]; i++) {
        trans->leg_trans_forward[i] = (double**) malloc(sizeof(double*) * (M+1));
        for (j = 0; j < M+1; j++)
            trans->leg_trans_forward[i][j] = (double*) malloc(sizeof(double) * (M+2-j));
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
        } else if (i > 0) {
            x_pos[m+i-1] = -trans->legendre_points[i];
        }
    }
    ind_count = 0;
    sind_count = 0;
    for (j = 0; j < M+1; j++) {
        pmn_vec = pmn_polynomial_value(NY, M+1, j, x_pos);

        for (i = 0; i < NY; i++) {
            if ((ind_count < ind[2]) && (i == ind[0]+ind_count*ind[1])) {
                for (k = j; k < M+2; k++) {
                    trans->leg_trans_forward[ind_count][j][k-j] = pmn_vec[i+k*NY]*weight[i];
                }
                ind_count++; 
            }

            if ((sind_count < sind[2]) && (j == sind[0]+sind_count*sind[1])) {
                for (k = j; k < M+2; k++) {
                    trans->leg_trans_backward[i][sind_count][k-j] = pmn_vec[i+k*NY];
                }
                sind_count++;
            }
        }

        free(pmn_vec);
    }

    free(weight);
    return 0;
}

int finalize_transform(legendre *trans)
{
    int i, j;

    free(trans->legendre_points);
    for (i = 0; i < trans->ind[2]; i++) {
        for (j = 0; j < M+1; j++)
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
}


