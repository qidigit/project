/*
 * =====================================================================================
 *
 *       Filename:  test.c
 *
 *    Description:  test constants
 *
 *        Version:  1.0
 *        Created:  05/14/2015 17:46:34
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
#include <mpi.h>
#include "model_constants.h"

int main(int argc, char **argv)
{
    int nproc, rproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);

    NY = 2048;
    initial_resolution();

    printf("proc %d initialized: NY = %d, M = %d, M_rsv = %d...\n", rproc, NY, M, M_rsv);
    printf("type %c\n", trunc_type);
    printf("damp_order = %d, damp_coeff=%.12f\n", damping_order, damping_coeff);
    printf("RADIUS = %.12f, OMEGA = %.12f, GRAV = %.12f, PI = %.12f\n\n", RADIUS, OMEGA, GRAV, PI);

    MPI_Finalize();
    return 0;
}
