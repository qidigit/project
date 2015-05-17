/*
 * =====================================================================================
 *
 *       Filename:  diagnostics.c
 *
 *    Description:  diagnositcs
 *
 *        Version:  1.0
 *        Created:  05/16/2015 13:37:20
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
#include "../barotropic/barotropic.h"
#include "diagnostics.h"

char ind_name[256];
char ug_name[256];
char vg_name[256];
char vorg_name[256];
char vors_name[256];

int diag_initial(char f_status)
{
    int rproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rproc);

    snprintf(ind_name, 256, "../model_data/ind_r%03d.bin", rproc);
    snprintf(ug_name, 256, "../model_data/u_r%03d.bin", rproc);
    snprintf(vg_name, 256, "../model_data/v_r%03d.bin", rproc);
    snprintf(vorg_name, 256, "../model_data/vorg_r%03d.bin", rproc);
    snprintf(vors_name, 256, "../model_data/vors_r%03d.bin", rproc);


    if (f_status == 'w') {
        FILE *pFile = NULL;
        pFile = fopen(ind_name, "wb");
        // write index of this process and the Gaussian quadrature points to this file
        fwrite(trans.ind, sizeof(trans.ind[0]), sizeof(trans.ind)/sizeof(trans.ind[0]), pFile);
        fwrite(trans.sind, sizeof(trans.sind[0]), sizeof(trans.sind)/sizeof(trans.sind[0]), pFile);
        fwrite(trans.legendre_points, sizeof(double), (NY+1)/2, pFile);
        fclose (pFile);

        // write initial fields
        pFile = fopen(ug_name, "wb");
        fwrite(u_gp.grid, sizeof(double), u_gp.ind[2]*(2*NY), pFile);
        fclose (pFile);

        pFile = fopen(vg_name, "wb");
        fwrite(v_gp.grid, sizeof(double), v_gp.ind[2]*(2*NY), pFile);
        fclose (pFile);

        pFile = fopen(vorg_name, "wb");
        fwrite(curl_gp.grid, sizeof(double), curl_gp.ind[2]*(2*NY), pFile);
        fclose (pFile);

    }

/*     pFile = fopen(vors_name, "wb");
 *     int m, n, cur_ind;
 *     double temp;
 *     for (m = 0; m < curl_gp.sind[2]; m++) {
 *         cur_ind = curl_gp.sind[0]+m*curl_gp.sind[1];
 *         for (n = 0; n < M+2-cur_ind; n++) {
 *             temp = creal(curl_gp.spec[m][n]);
 *             fwrite(&temp, sizeof(double), 1, pFile);
 *             temp = cimag(curl_gp.spec[m][n]);
 *             fwrite(&temp, sizeof(double), 1, pFile);
 *         }
 *     }
 *     fclose (pFile);
 */
    return 0;
}

int diag_writeData()
{
    FILE *pFile = NULL;
    // write fields
    pFile = fopen(ug_name, "ab");
    fwrite(u_gp.grid, sizeof(double), u_gp.ind[2]*(2*NY), pFile);
    fclose(pFile);

    pFile = fopen(vg_name, "ab");
    fwrite(v_gp.grid, sizeof(double), v_gp.ind[2]*(2*NY), pFile);
    fclose(pFile);

    pFile = fopen(vorg_name, "ab");
    fwrite(curl_gp.grid, sizeof(double), curl_gp.ind[2]*(2*NY), pFile);
    fclose(pFile);

    return 0;
}
int diag_readData()
{
    FILE *pFile = NULL;

    pFile = fopen(ug_name, "rb");
    fread(u_gp.grid, sizeof(double), u_gp.ind[2]*(2*NY), pFile);
    fclose (pFile);

    pFile = fopen(vg_name, "rb");
    fread(v_gp.grid, sizeof(double), v_gp.ind[2]*(2*NY), pFile);
    fclose (pFile);

    pFile = fopen(vorg_name, "rb");
    fread(curl_gp.grid, sizeof(double), curl_gp.ind[2]*(2*NY), pFile);
    fclose(pFile);

    return 0;
}

