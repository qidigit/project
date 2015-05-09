/*
 * =====================================================================================
 *
 *       Filename:  test1.c
 *
 *    Description:  
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
#include <complex.h>
#include <fftw3.h>
#include <omp.h>
#include <mpi.h>
#include "field_management.h"
#include "transform_s2.h"

const int M = 100, NY = 128;
const char trunc_type = 't';
const int M_res = 85;

int error_func() {
    int ind[3]={0,1,20}, sind[3]={0,2,6};
    int status = 1, status_f = 1;

    int rproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);


    sphere_grid xi_g;
    sphere_spec xi_s;
    sphere_four xi_f;
    status_f = initialize_fields(ind, sind, &xi_g, &xi_f, &xi_s);

    legendre trans;
    status=initialize_transform(ind, sind, &trans);
    
    if (status == 1 || status_f == 1) {
        printf("initilization failure!\n");
        return(EXIT_FAILURE);
    }    

    time_t t;
    srand((unsigned) time(&t));
    int i;
    int size = (NY*2) * xi_g.ind[2];
    double input;
    for (i = 0; i < (NY*2) * xi_g.ind[2]; i++) {
        input = 1 / (1 + (double) i*i/size/size);
        xi_g.u_grid[i] = input;
//        input = rand() % (int) 1E10;
//        xi_g.u_grid[i] = (double) input / 1E10;

    }

    double tic, toc, elapsed;
    tic = omp_get_wtime();
    transform_g2s(&xi_g, &xi_f, &xi_s, &trans);
    double xi_ori[(2*NY)*xi_g.ind[2]];
    for (i = 0; i < (2*NY)*xi_g.ind[2]; i++) xi_ori[i] = xi_g.u_grid[i];
    complex xi_hat[(NY+1)*xi_f.ind[2]];
    for (i = 0; i < (NY+1)*xi_f.ind[2]; i++) xi_hat[i] = xi_f.u_four[i];

    transform_s2g(&xi_s, &xi_f, &xi_g, &trans);
    toc = omp_get_wtime();

    if (rproc == 0) {
    double error[size];
    for (i = 0; i < size; i++) error[i] = xi_g.u_grid[i]-xi_ori[i];
    elapsed = toc-tic;
    printf("\nrun_time = %f seconds.\n", elapsed);

    int start = 0, col = 0;
    printf("\n  input           fourier         output          error\n");
    for (i = 0; i < 20; i++) {
        printf("%12g    %12g+%12gi  %12g    %12g\n", xi_ori[col*2*NY+start+i],creal(xi_hat[col*(NY+1)+start+i]),cimag(xi_hat[col*(NY+1)+start+i]), xi_g.u_grid[col*2*NY+start+i], error[col*(NY+1)+start+i]);
    }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    /*  
    for (i = 0; i < (NY+1)/2; i++) {
        printf("x[%d]=%10.10e\n", i, trans.legendre_points[i]);
    }
    int m_prt = 0;
    int n_prt = 0-m_prt;
    printf ( "\n     N     M        X               Pmn(N,M,X)\n" );
    for (i = 0; i < ind[2]; i++) {
        printf("  %4d  %4d  %12g  %24.16g  %24.16g\n", n_prt+m_prt, m_prt, 
                trans.legendre_points[ind[0]+i*ind[1]], 
                trans.leg_trans_forward[m_prt][n_prt][i], 
                trans.leg_trans_backward[ind[0]+i*ind[1]][(m_prt-sind[0])/sind[1]][n_prt]);
    }

    int j, k;
    int count = 0;
    for (i = 0; i < M+1; i++) {
        for (j = 0; j < M+2-i; j++) {
            for (k = 0; k < ind[2]; k++) {
                if (trans.leg_trans_forward[i][j][k] == 0) {
                    printf("non-initialized array at (i,j,k)=(%d,%d,%d)\n", i, j, k);
                } else {
                    count++;
                }
            }
        }
    }
    printf("non-zero entries = %d, size of array = %lu\n", count, sizeof(trans.leg_trans_forward)/sizeof(double));
    count = 0;
    for (i = 0; i < NY; i++) {
        for (j = 0; j < sind[2]; j++) {
            for (k = 0; k < M+2-(sind[0]+j*sind[1]); k++) {
                if (trans.leg_trans_backward[i][j][k] == 0) {
                    printf("non-initialized array at (i,j,k)=(%d,%d,%d)\n", i, j, k);
                } else {
                    count++;
                }
            }
        }
    }
    printf("non-zero entries = %d, size of array = %lu\n", count, sizeof(trans.leg_trans_backward)/sizeof(double));
    */

    status = finalize_transform(&trans);
    status_f = finalize_fields(&xi_g, &xi_f, &xi_s);
    if (status == 1 || status_f == 1) {
        printf("initilization failure!\n");
        return(EXIT_FAILURE);
    }
   
    return 0; 
}

void main(int argc, char **argv) {

    int nproc, rproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);

    printf("\n rank %d return status %d\n", rproc, error_func());

    MPI_Finalize(); 
//    sphere_grid xi_g = error_func();
//    printf("xi_grid begin =  %d\t, end = %d\n", xi_g.ind[0], xi_g.ind[1] );

}
