/*
 * =====================================================================================
 *
 *       Filename:  test2.c
 *
 *    Description:  test2
 *
 *        Version:  1.0
 *        Created:  04/22/2015 16:04:56
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

void main() {
    int i;
    int NY =13;
    int m = (NY+1)>>1;
    double x_pos[NY];
    double x[m];

    for (i = 0; i < m; i++)
        x[i] = i;

    for (i = 0; i < m; i++) {
        x_pos[i] = x[i];
        if (NY%2 == 0) {
            x_pos[m+i] = -x[i];
        } else if (i > 0) {
            x_pos[m+i-1] = -x[i];
        }
    }

    for (i = 0; i < NY; i++) {
        printf("x_pos[%d] = %f.\n", i+1, x_pos[i]);
    }

    return;
}
