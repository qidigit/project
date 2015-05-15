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

// diffusion & damping
extern int damping_order;
extern double damping_coeff;

// physical constants
extern const double RADIUS;
extern const double OMEGA;
extern const double GRAV;
extern const double RDGAS;
extern const double KAPPA;
extern const double CP_AIR; 
extern const double CP_OCEAN;
extern const double RHO0;
extern const double RHO0R;
extern const double RHO_CP;

// other important constants
extern const double PI;
extern const double DELVIN;

int initial_resolution();
#endif /* MODEL_CONSTANTS */
