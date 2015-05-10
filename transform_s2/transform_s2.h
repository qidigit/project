/*
 * =====================================================================================
 *
 *       Filename:  transform_s2.h
 *
 *    Description:  This is the header file for functions conducting FFTs on the 2-sphere.
 *
 *        Version:  1.0
 *        Created:  04/21/2015 11:52:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _TRANSFORM_S2_H
#define _TRANSFORM_S2_H

// model size defined from the input
extern const int M, NY;
extern const char trunc_type;
extern const int M_rsv; // the real number of resolved modes (due to the equi-partition of modes in each processor)

// setup for gaussian quadrature points precision 
#define EPS 1E-10

/* struc legendre defines the legendre transform in latitude direction. 
 * legendre_points defines the NY gauss-legendre quadrature points (with symmetry between 0);
 * leg_trans_forward defines P_mn(mu_i)*w(mu_i) the forward legendre transform matrix, of size (ind_begin:ind_end) * (M+1) * (Nm+1);
 * leg_trans_backward defines P_mn(mu_i) the backward legendre transform of size NY * (sind_begin:sind_end) * (Nm+1).
 * */

/* Convention of indexing ind[3] & sind[3]
 * ind takes values in {0,1,...,NY-1} is the index for latitude grid number corresponding to the latitude grids mu_(ind). ind[0] = index_begin, ind[1] = index_increasement>=1, ind[2] = index_length>0.
 * sind takes values in {0,1,...,M} is the positive wavenumbers in zonal direction.
 * */
typedef struct {
    int ind[3];
    int sind[3];

    double *legendre_points;          // mu_j
    double ***leg_trans_forward;      // P_mn(mu_j)*w(mu_j)
    double ***leg_trans_backward;     // P_mn(mu_j)
} legendre;

int initialize_transform(int ind[], int sind[], legendre *trans);

int transform_g2s(sphere_grid *xi_g, sphere_four *xi_f, sphere_spec *xi_s, legendre *trans);
int transform_s2g(sphere_spec *xi_s, sphere_four *xi_f, sphere_grid *xi_g, legendre *trans);

int finalize_transform(legendre *trans);

#endif /* TRANSFORM_S2_H */
