/*
 * =====================================================================================
 *
 *       Filename:  constants.h
 *
 *    Description:  This file defines useful constants for Earth
 *
 *        Version:  1.0
 *        Created:  05/14/2015 17:35:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _MODEL_CONSTANTS_H
#define _MODEL_CONSTANTS_H

// Model setups: define model resolution and model truncation method
extern int M, NY; // the number of resolved zonal waves and meridional grid points
extern int M_rsv; // actual resolved number of zonal waves
extern char trunc_type; // truncation type
extern double TT; // total time to calculate

// diffusion & damping
extern int damping_order;
extern double damping_coeff;

// integration time step
extern double delta_t;
// Robert filter weight
extern double robert;

// physical constants
#define RADIUS 6371.0e3
#define OMEGA 7.292e-5
#define GRAV 9.80
#define RDGAS 287.04
#define KAPPA 2.0/7.0
#define CP_AIR 287.04/2.0*7.0 
#define CP_OCEAN 3989.24495292815
#define RHO0 1.035e3
#define RHO0R 1.0/1.035e3
#define RHO_CP 1.035e3*3989.24495292815

// other important constants
#define PI 3.1415926535897932384626433832795028841971693993751
#define KELVIN 273.15


int initial_resolution();
#endif /* MODEL_CONSTANTS */
