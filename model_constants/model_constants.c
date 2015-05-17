/*
 * =====================================================================================
 *
 *       Filename:  constants.c
 *
 *    Description:  define useful constants
 *
 *        Version:  1.0
 *        Created:  05/14/2015 17:38:42
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
#include <math.h>
#include <mpi.h>
#include "model_constants.h"

int NY = 256;
int M = 171;
int M_rsv = 170;
char trunc_type = 't';
double TT = 1; // number of days to run

int damping_order = 4;
double damping_coeff = 1.0e-3;

double delta_t = 1800.0;
double delta_at = 1800.0;
double robert = 0.04;


/* ------------ physical constants ---------------
 *  <DATA NAME="RADIUS" UNITS="m" TYPE="real" DEFAULT="6371.e3">
 *    radius of the earth
 *  </DATA>
 *  <DATA NAME="OMEGA" UNITS="1/s" TYPE="real" DEFAULT="7.292e-5">
 *    rotation rate of the planet (earth)
 *  </DATA>
 *  <DATA NAME="GRAV" UNITS="m/s^2" TYPE="real" DEFAULT="9.80">
 *    acceleration due to gravity
 *  </DATA>
 *  <DATA NAME="RDGAS" UNITS="J/kg/deg" TYPE="real" DEFAULT="287.04">
 *    gas constant for dry air
 *  </DATA>
 *  <DATA NAME="KAPPA" TYPE="real" DEFAULT="2./7.">
 *    RDGAS / CP_AIR
 *  </DATA>
 *  <DATA NAME="CP_AIR" UNITS="J/kg/deg" TYPE="real" DEFAULT="RDGAS/KAPPA">
 *    specific heat capacity of dry air at constant pressure
 *  </DATA>
 *  <DATA NAME="CP_OCEAN" UNITS="J/kg/deg" TYPE="real" DEFAULT="3989.24495292815">
 *    specific heat capacity taken from McDougall (2002) "Potential Enthalpy ..."
 *  </DATA>
 *  <DATA NAME="RHO0" UNITS="kg/m^3" TYPE="real" DEFAULT="1.035e3">
 *    average density of sea water
 *  </DATA>
 *  <DATA NAME="RHO0R" UNITS="m^3/kg" TYPE="real" DEFAULT="1.0/RHO0">
 *    reciprocal of average density of sea water
 *  </DATA>
 *  <DATA NAME="RHO_CP" UNITS="J/m^3/deg" TYPE="real" DEFAULT="RHO0*CP_OCEAN">
 *    (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C)
 *  </DATA>
 */

/* const double RADIUS = 6371.0e3;
 * const double OMEGA  = 7.292e-5;
 * const double GRAV   = 9.80;
 * const double RDGAS  = 287.04;
 * const double KAPPA  = 2.0/7.0;
 * const double CP_AIR = 287.04/2.0*7.0; 
 * const double CP_OCEAN = 3989.24495292815;
 * const double RHO0    = 1.035e3;
 * const double RHO0R   = 1.0/1.035e3;
 * const double RHO_CP  = 1.035e3*3989.24495292815;
 */


/* const double PI = 3.1415926535897932384626433832795028841971693993751;
 * const double KELVIN = 273.15;
 */


// initialize distributed modes and resolved modes M, M_rsv from NY
int initial_resolution()
{
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    M_rsv = (2 * NY - 1) / 3;
    M = (int) ceil((double) (M_rsv+1)/nproc) * nproc - 1; 
    return 0;
}
