/*
 * =====================================================================================
 *
 *       Filename:  field_management.c
 *
 *    Description:  This file contains functions that initialize and finalize the allocated fields. FFTW plan is created in the initialization step for each structure (may not be efficient for multi-field with similar fields.
 *
 *        Version:  1.0
 *        Created:  04/21/2015 12:48:40
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
#include <fftw3.h>
#include <omp.h>
#include "../model_constants/model_constants.h"
#include "field_management.h"


int initialize_fields(int ind[], int sind[], sphere_grid *xi_g, sphere_four *xi_f, sphere_spec *xi_s)
{
    int i; // iteration index

    // initialize grid field
    xi_g->u_grid = (double*) fftw_malloc(sizeof(double) * (2*NY) * ind[2]);
    xi_g->ind[0] = ind[0];
    xi_g->ind[1] = ind[1];
    xi_g->ind[2] = ind[2];

    if (xi_g->u_grid == NULL) {
        printf("grid field initilization failure!\n");
        return(EXIT_FAILURE);
    }


    // initialize fourier field
    xi_f->u_four = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (NY+1) * ind[2]);
    xi_f->ind[0] = ind[0];
    xi_f->ind[1] = ind[1];
    xi_f->ind[2] = ind[2];

    if (xi_f->u_four == NULL) {
        printf("fourier field initilization failure!\n");
        return(EXIT_FAILURE);
    }


    // initialize FFTW plan
    initialize_fft_plan(xi_g, xi_f);


    // initialize spectral field
    switch(trunc_type)
    {
    case 'r' :
        xi_s->u_spec = (double complex**) malloc(sizeof(double complex*) * sind[2]);
        if (xi_s->u_spec == NULL) {
            printf("spectral field initilization failure!\n");
            return(EXIT_FAILURE);
        }

        for (i = 0; i < sind[2]; i++) {
            xi_s->u_spec[i] = (double complex*) malloc(sizeof(double complex) * (M+2));
            if (xi_s->u_spec[i] == NULL) {
                printf("spectral field initilization failure!\n");
                return(EXIT_FAILURE);
            }
        }
        break;
    case 't' :
        xi_s->u_spec = (double complex**) malloc(sizeof(double complex*) * sind[2]);
        if (xi_s->u_spec == NULL) {
            printf("spectral field initilization failure!\n");
            return(EXIT_FAILURE);
        }

        for (i = 0; i < sind[2]; i++) {
            xi_s->u_spec[i] = (double complex*) malloc(sizeof(double complex) * (M+2-sind[0]-i*sind[1]));
            if (xi_s->u_spec[i] == NULL) {
                printf("spectral field initilization failure!\n");
                return(EXIT_FAILURE);
            }
        }
        break;
    default :
        printf("Invalid trunctation type, must be 't' or 'r'\n");
        break;
    }
    xi_s->sind[0] = sind[0];
    xi_s->sind[1] = sind[1];
    xi_s->sind[2] = sind[2];

    return 0;
}


int finalize_fields(sphere_grid *xi_g, sphere_four *xi_f, sphere_spec *xi_s)
{
    fftw_destroy_plan(xi_g->p_forward);
    fftw_destroy_plan(xi_f->p_backward);
    fftw_free(xi_f->u_four);
    fftw_free(xi_g->u_grid);

    int i;
    for (i = 0; i < xi_s->sind[2]; i++)
        free(xi_s->u_spec[i]);
    free(xi_s->u_spec);

    return 0;
}


int initialize_group_fft_plan(trans_group *xi)
{
    int size = 2*NY;
    int howmany = xi->ind[2];
    int idist = size, odist = size/2+1;
    int istride = 1, ostride = 1;
    int inembed = size, onembed = size;

    // initialize omp for parallel
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    printf("Start OpenMP for FFT, Number of Threads = %d\n\n", omp_get_max_threads());

    xi->gp_forward  = fftw_plan_many_dft_r2c(1,         &size,      howmany,
                                             xi->grid,  &inembed,   istride,   idist,
                                             xi->four,  &onembed,   ostride,   odist,
                                             FFTW_MEASURE);
    xi->gp_backward = fftw_plan_many_dft_c2r(1,         &size,      howmany,
                                             xi->four,  &onembed,   ostride,   odist,
                                             xi->grid,  &inembed,   istride,   idist,
                                             FFTW_MEASURE);


    return 0;
}

int initialize_fft_plan(sphere_grid *xi_g, sphere_four *xi_f)
{
    int size = 2*NY;
    int howmany = xi_g->ind[2];
    int idist = size, odist = size/2+1;
    int istride = 1, ostride = 1;
    int inembed = size, onembed = size;

    // initialize omp for parallel
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    printf("Start OpenMP for FFT, Number of Threads = %d\n\n", omp_get_max_threads());

    xi_g->p_forward  = fftw_plan_many_dft_r2c(1,            &size,      howmany,
                                             xi_g->u_grid,  &inembed,   istride,   idist,
                                             xi_f->u_four,  &onembed,   ostride,   odist,
                                             FFTW_MEASURE);
    xi_f->p_backward = fftw_plan_many_dft_c2r(1,            &size,      howmany,
                                             xi_f->u_four,  &onembed,   ostride,   odist,
                                             xi_g->u_grid,  &inembed,   istride,   idist,
                                             FFTW_MEASURE);


    return 0;
}




int initialize_trans_group(int ind[], int sind[], trans_group *xi)
{
    int i; // iteration index

    // initialize grid field
    xi->grid = (double*) fftw_malloc(sizeof(double) * (2*NY) * ind[2]);
    xi->ind[0] = ind[0];
    xi->ind[1] = ind[1];
    xi->ind[2] = ind[2];

    if (xi->grid == NULL) {
        printf("group grid field initilization failure!\n");
        return(EXIT_FAILURE);
    }


    // initialize fourier field
    xi->four = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (NY+1) * ind[2]);
    if (xi->four == NULL) {
        printf("group fourier field initilization failure!\n");
        return(EXIT_FAILURE);
    }


    // initialize FFTW plan
    initialize_group_fft_plan(xi);


    // initialize spectral field
    switch(trunc_type)
    {
    case 'r' :
        xi->spec = (double complex**) malloc(sizeof(double complex*) * sind[2]);
        if (xi->spec == NULL) {
            printf("group spectral field initilization failure!\n");
            return(EXIT_FAILURE);
        }

        for (i = 0; i < sind[2]; i++) {
            xi->spec[i] = (double complex*) malloc(sizeof(double complex) * (M+2));
            if (xi->spec[i] == NULL) {
                printf(" group spectral field initilization failure!\n");
                return(EXIT_FAILURE);
            }
        }
        break;
    case 't' :
        xi->spec = (double complex**) malloc(sizeof(double complex*) * sind[2]);
        if (xi->spec == NULL) {
            printf("group spectral field initilization failure!\n");
            return(EXIT_FAILURE);
        }

        for (i = 0; i < sind[2]; i++) {
            xi->spec[i] = (double complex*) malloc(sizeof(double complex) * (M+2-sind[0]-i*sind[1]));
            if (xi->spec[i] == NULL) {
                printf("group spectral field initilization failure!\n");
                return(EXIT_FAILURE);
            }
        }
        break;
    default :
        printf("Invalid trunctation type, must be 't' or 'r'\n");
        break;
    }
    xi->sind[0] = sind[0];
    xi->sind[1] = sind[1];
    xi->sind[2] = sind[2];

    return 0;
}


int finalize_trans_group(trans_group *xi)
{
    fftw_destroy_plan(xi->gp_forward);
    fftw_destroy_plan(xi->gp_backward);
    fftw_free(xi->four);
    fftw_free(xi->grid);

    int i;
    for (i = 0; i < xi->sind[2]; i++)
        free(xi->spec[i]);
    free(xi->spec);

    return 0;
}


int initialize_grid_field(int ind[], grid_field *u_g)
{
    // initialize grid field
    u_g->grid = (double*) fftw_malloc(sizeof(double) * (2*NY) * ind[2]);
    u_g->ind[0] = ind[0];
    u_g->ind[1] = ind[1];
    u_g->ind[2] = ind[2];

    if (u_g->grid == NULL) {
        printf("grid field initilization failure!\n");
        return(EXIT_FAILURE);
    }

    return 0;
}


int finalize_grid_field(grid_field *u_g)
{
    fftw_free(u_g->grid);

    return 0;
}

int initialize_spec_field(int sind[], spec_field *u_s)
{
    int i; // iteration index

    // initialize spectral field
    switch(trunc_type)
    {
    case 'r' :
        u_s->spec = (double complex**) malloc(sizeof(double complex*) * sind[2]);
        if (u_s->spec == NULL) {
            printf("spectral field initilization failure!\n");
            return(EXIT_FAILURE);
        }

        for (i = 0; i < sind[2]; i++) {
            u_s->spec[i] = (double complex*) malloc(sizeof(double complex) * (M+2));
            if (u_s->spec[i] == NULL) {
                printf("spectral field initilization failure!\n");
                return(EXIT_FAILURE);
            }
        }
        break;
    case 't' :
        u_s->spec = (double complex**) malloc(sizeof(double complex*) * sind[2]);
        if (u_s->spec == NULL) {
            printf("spectral field initilization failure!\n");
            return(EXIT_FAILURE);
        }

        for (i = 0; i < sind[2]; i++) {
            u_s->spec[i] = (double complex*) malloc(sizeof(double complex) * (M+2-sind[0]-i*sind[1]));
            if (u_s->spec[i] == NULL) {
                printf("spectral field initilization failure!\n");
                return(EXIT_FAILURE);
            }
        }
        break;
    default :
        printf("Invalid trunctation type, must be 't' or 'r'\n");
        break;
    }
    u_s->sind[0] = sind[0];
    u_s->sind[1] = sind[1];
    u_s->sind[2] = sind[2];

    return 0;
}


int finalize_spec_field(spec_field *u_s)
{
    int i;
    for (i = 0; i < u_s->sind[2]; i++)
        free(u_s->spec[i]);
    free(u_s->spec);

    return 0;
}


