/*
 * =====================================================================================
 *
 *       Filename:  field_management.h
 *
 *    Description:  This is the header file to create and destroy fields.
 *
 *        Version:  1.0
 *        Created:  04/21/2015 12:13:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _FIELD_MANAGEMENT_H
#define _FIELD_MANAGEMENT_H

// model size defined from the input
extern const int M, NY;
extern const char trunc_type;


// define the fields in physical and spectral domains
/* ind[0] = index_start, ind[1] = increment, ind[2] = size of stencil.
 * */

// field in physical domain, of size (2*NY*ind[2]) 1d.
typedef struct {
    double *u_grid;
    fftw_plan p_forward;
    int ind[3];
} sphere_grid;

// field in spectral domain under spherical harmonics, of size (sind[2]) * (M+2-m) 2d, the spherical modes is of wavenumber m:M+1 according to zonal mode m>=0.
typedef struct {
    complex **u_spec;
    int sind[3];
} sphere_spec;

// field with zonal direction in Fourier domain, of size ((NY+1)*ind[2]) 1d, note that only the first M+1 (0:M) modes are actually of use and the rest part is set to be 0 for de-aliasing.
typedef struct {
    fftw_complex *u_four;
    fftw_plan p_backward;
    int ind[3];
} sphere_four;


int initialize_fields(int ind[], int sind[], sphere_grid *xi_g, sphere_four *xi_f, sphere_spec *xi_s);

int finalize_fields(sphere_grid *xi_g, sphere_four *xi_f, sphere_spec *xi_s);

int initialize_fft_plan(sphere_grid *xi_g, sphere_four *xi_f);

#endif /* FIELD_MANAGEMENT_H */
