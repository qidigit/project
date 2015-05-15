/*
 * =====================================================================================
 *
 *       Filename:  spectral_diff.h
 *
 *    Description:  This part of codes calculates the divergence and curl of the vector field (u, v) from the streamfunctions and also the inverse transform.
 *
 *        Version:  1.0
 *        Created:  05/13/2015 20:30:34
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _SPECTRAL_DIFF_H
#define _SPECTRAL_DIFF_H


// notion for variables
/* stream function     psi -- with normalization a^2
 * potential function  xi  -- with normalization a^2
 * flow vairables      (u, v) -- with normalization a
 * divergence and rotation (div, curl)
 */

// This function calculates the divergence of vector field (u, v) in spectral domain
int div_spec_from_uv_grid(trans_group *u_gp, trans_group *v_gp, legendre *trans, spec_field *div_s);

// This function calculate the vector field (u, v) in physical domain from its divergence and curl
int uv_grid_from_div_curl(spec_field *div_s, spec_field *curl_s, legendre *trans, trans_group *u_gp, trans_group *v_gp);

// These two functions calculate the forward & backward Laplacian operator in spectral domain
// Note these functions can be only applied to psi, xi (or div, curl)
int laplace_forward(spec_field *psi_in, spec_field *curl_out);
int laplace_backward(spec_field *curl_in, spec_field *psi_out);
#endif /* SPECTRAL_DIFF_H */
