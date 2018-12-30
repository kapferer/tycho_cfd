/*
 * para_coef.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Here the coefficients for the parabel functions are determined.
 */
int para_coef(int nmin, int nmax, int flag1, double *a_coef, double *dx,
        double *ai_coef, double *b_coef, double *bi_coef, double *c_coef,
        double *d_x, double **para, double *ci_coef) {
    int i;

    for (i = nmin - 2; i <= nmax; i++) {
        a_coef[i] = dx[i] + dx[i + 1];
        //printf("%i %i %e %e\n", flag1, i, dx[i], dx[i+1]);
        ai_coef[i] = 1.0 / a_coef[i];
        b_coef[i] = a_coef[i] + dx[i];
        bi_coef[i] = 1.0 / b_coef[i];
        c_coef[i] = a_coef[i] + dx[i + 1];
        ci_coef[i] = 1.0 / c_coef[i];
    }
    //a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * da(j+1) + C3 * da(j)
    for (i = nmin - 1; i <= nmax; i++) {
        d_x[i] = 1.0 / (a_coef[i - 1] + a_coef[i + 1]);
        para[i][0] = dx[i] * ai_coef[i] + 2.0 * dx[i + 1] * dx[i] * d_x[i] * ai_coef[i] * (a_coef[i - 1] * bi_coef[i] - a_coef[i + 1] * ci_coef[i]);

        // to debug
        //printf("%i %i %f %f %f\n", i, flag1, d_x[i], a_coef[i-1], a_coef[i+1]);
        //printf("in para_coef: %f %f %f %f\n %f %f %f %f %f %i flag:%i\n", para[i][0], dx[i],ai_coef[i],dx[i+1], d_x[i], a_coef[i-1], bi_coef[i], a_coef[i+1], ci_coef[i], i, flag1);

        para[i][1] = -1 * d_x[i] * dx[i] * a_coef[i - 1] * bi_coef[i];
        para[i][2] = d_x[i] * dx[i + 1] * a_coef[i + 1] * ci_coef[i];
    }

    //da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))
    for (i = nmin - 1; i <= nmax; i++) {
        d_x[i] = dx[i] / (a_coef[i - 1] + dx[i + 1]);
        para[i][3] = d_x[i] * b_coef[i - 1] * ai_coef[i];
        para[i][4] = d_x[i] * c_coef[i] * ai_coef[i - 1];
    }

    // to debug
    //for (n=0; n<max_array_length+12; n++) printf("in para_coef: %g %g %g %g %g %i where:%i\n",  para[n][0], para[n][1], para[n][2], para[n][3], para[n][4], n, flag1);

    return 0;
}
