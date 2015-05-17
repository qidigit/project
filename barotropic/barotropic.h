/*
 * =====================================================================================
 *
 *       Filename:  barotropic.h
 *
 *    Description:  This is the main part to calculate the barotropic flow on the sphere
 *
 *        Version:  1.0
 *        Created:  05/15/2015 16:43:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _BAROTROPIC_H
#define _BAROTROPIC_H

extern trans_group u_gp, v_gp, curl_gp;
extern spec_field delta_curl, curl_pre, curl_post;
extern spec_field psi_s, pot_s;
extern grid_field freq;
extern legendre trans;

int initialize_barotropic();
int finalize_barotropic();

// set initial value for barotropic model input
int set_initial_value(int num);

// one time integration for the barotropic equation from t-dt, t --> t+dt
int barotropic_integration();
int barotropic_update();
#endif /* BAROTROPIC_H */
