/*
 * =====================================================================================
 *
 *       Filename:  time_inte.h
 *
 *    Description:  This function runs the leapfrog time integration
 *
 *        Version:  1.0
 *        Created:  05/15/2015 16:18:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _TIME_INTE_H
#define _TIME_INTE_H

// this function calculates the leapfrog time integration.
// Input: xi_pre at time t-dt;
//        delta_xi = (xi_post - xi_pre) / (2 * dt);
//        dt.
// Output: xi_post at time t+dt.
int leapfrog_time_inte(spec_field *xi_pre, spec_field *delta_xi, double dt, spec_field *xi_post);
#endif /* TIME_INTE_H */
