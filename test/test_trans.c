/*
 * =====================================================================================
 *
 *       Filename:  test_trans.c
 *
 *    Description:  test run function to check the transformation between physical domain and spectral domain.
 *
 *        Version:  1.0
 *        Created:  04/21/2015 13:10:00
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
#include <time.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <omp.h>
#include <mpi.h>
#include "../field_management/field_management.h"
#include "../transform_s2/transform_s2.h"

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif

#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif


const int NY = 256;
const int M_rsv = 170; // (NY*2-1) / 3;
const int M = 173;
const char trunc_type = 't';

int error_func() {
    int i, j, ll;
    int status_t = 1, status_f = 1;
    int ind[3], sind[3];

    int nproc, rproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // define the number of transfered modes M, required to be a multiple of nproc.
//    M = (int) ceil((double) M_rsv/nproc) * nproc - 1;
    printf("\nThe setting-ups for the computation is NY = %d, M_rsv = %d, M = %d\n", NY, M_rsv, M);

    // set up the splitted modes in each proc
    ind[0] = rproc*(NY/nproc); ind[1] = 1; ind[2] = NY / nproc;
    sind[0] = rproc; sind[1] = nproc; sind[2] = (M+1) / nproc;
    printf("proc rank %d is assigned index block:\n", rproc);
    printf("ind  = [%d %d %d]\n",  ind[0],  ind[1],  ind[2]);
    printf("sind = [%d %d %d]\n", sind[0], sind[1], sind[2]);

    MPI_Barrier(MPI_COMM_WORLD);

    sphere_grid xi_g;
    sphere_spec xi_s;
    sphere_four xi_f;
    status_f = initialize_fields(ind, sind, &xi_g, &xi_f, &xi_s);

    legendre trans;
    status_t=initialize_transform(ind, sind, &trans);
    
    if (status_t == 1 || status_f == 1) {
        printf("initilization failure!\n");
        return(EXIT_FAILURE);
    }    

//    time_t t;
//    srand((unsigned) time(&t));
    int size = (NY*2);
    int cur_ind;
    double input;
    double theta;
    for (j = 0; j < xi_g.ind[2]; j++) {
        cur_ind = xi_g.ind[0]+j*xi_g.ind[1];
        if (cur_ind < (NY/2)) {
            theta = asin(trans.legendre_points[cur_ind]);
        } else {
            // in this setting, NY must be even
            theta = asin(-trans.legendre_points[cur_ind-(NY/2)]);
        }
        for (i = 0; i < (NY*2); i++) {
            input = 25 * cos(theta) - 30 * pow(cos(theta), 3) + 300 * pow(sin(theta), 2) * pow(cos(theta), 6);
            xi_g.u_grid[j*(NY*2)+i] = input * cos(5*i*2*PI/(2*NY)) * sin(10*i*2*PI/(2*NY));
        }
//        input = rand() % (int) 1E10;
//        xi_g.u_grid[i] = (double) input / 1E10;

    }

    double tic, toc, elapsed;
    double xi_ori[(2*NY)*xi_g.ind[2]];
    double xi_tran[(2*NY)*xi_g.ind[2]];
    double xi_back[(2*NY)*xi_g.ind[2]];

    // set printing parameters
    int start = NY/20, col = 0;

    for (i = 0; i < (2*NY)*xi_g.ind[2]; i++) xi_ori[i] = xi_g.u_grid[i];
    tic = omp_get_wtime();

    transform_g2s(&xi_g, &xi_f, &xi_s, &trans);

//    for (i = 0; i < NY/nproc; i++) {
//        printf("%12g+%12gi\n", creal(xi_f.u_four[i*(NY+1)+start]), cimag(xi_f.u_four[i*(NY+1)+start]));
//    }

    transform_s2g(&xi_s, &xi_f, &xi_g, &trans);

//    printf("\n  input               trans           spec\n");
//    for (i = 0; i < M+2; i++) {
//        printf("%12g    %12g+%12gi  %12g\n", xi_g.u_grid[i*2*NY+start],creal(xi_f.u_four[i*(NY+1)+start]), cimag(xi_f.u_four[i*(NY+1)+start]), cimag(xi_s.u_spec[start][i]));
//    }

    for (i = 0; i < (2*NY)*xi_g.ind[2]; i++) xi_tran[i] = xi_g.u_grid[i];

    transform_g2s(&xi_g, &xi_f, &xi_s, &trans);
    transform_s2g(&xi_s, &xi_f, &xi_g, &trans);
    for (i = 0; i < (2*NY)*xi_g.ind[2]; i++) xi_back[i] = xi_g.u_grid[i];
    toc = omp_get_wtime();

    
    double error[(2*NY)*xi_g.ind[2]];
    double error1[(2*NY)*xi_g.ind[2]];
    double terr = 0, terr1 = 0;
    for (i = 0; i < (2*NY)*xi_g.ind[2]; i++) error[i] = xi_back[i]-xi_tran[i];
    for (i = 0; i < (2*NY)*xi_g.ind[2]; i++) error1[i] = (xi_back[i]-xi_ori[i]) / xi_ori[i];
    elapsed = toc-tic;
   
    for (j = 0; j < nproc; j++) {
    if (j == rproc) {

    for (ll = 0; ll < (2*NY)*xi_g.ind[2]; ll++) terr  += error[ll];
    for (ll = 0; ll < (2*NY)*xi_g.ind[2]; ll++) terr1 += error1[ll];

    printf("\nproc rank %d output:\n", rproc);
    printf("run_time = %f seconds.\n\n", elapsed);

    if (j == 0) {

    printf("\n  input               trans            back             error         error1\n");
    for (i = 0; i < NY/nproc; i+=1) {
    
            printf("%12g    %12g    %12g    %12g    %12g\n", xi_ori[i*2*NY+start],xi_tran[i*2*NY+start], xi_back[i*2*NY+start], error[i*2*NY+start], error1[i*2*NY+start]);
    }

    }

    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    }

    if (rproc == 0) {
        printf("\ntotal error = %f, total rela. error = %f\n", terr/nproc, terr1/nproc);
    }
    
    status_t = finalize_transform(&trans);
    status_f = finalize_fields(&xi_g, &xi_f, &xi_s);
    if (status_t == 1 || status_f == 1) {
        printf("finalize failure!\n");
        return(EXIT_FAILURE);
    }
   
    return 0; 
}

int main(int argc, char **argv) {

    int nproc, rproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);

    printf("rank %d return status %d\n", rproc, error_func());

    MPI_Finalize(); 
//    sphere_grid xi_g = error_func();
//    printf("xi_grid begin =  %d\t, end = %d\n", xi_g.ind[0], xi_g.ind[1] );
//
    return 0;

}
