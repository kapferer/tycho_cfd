/*
 * flatten.c
 *
 *      Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Inline function to calculate the maximum of two doubles
 */
inline double max_flatten(double a, double b) {
    double tmp = 0.0;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 Inline function to calculate the minimum of two doubles
 */
inline double min_flatten(double a, double b) {
    double tmp = 0.0;

    if (a < b) tmp = a;
    if (a == b)tmp = a;
    if (a > b) tmp = b;

    return tmp;
}

/*!
 Look for presence of a shock using pressure gradient and sign of
 velocity jump:  shock = 1 if there is a shock in the zone, else shock = 0
 Compute steepness parameter based on steepness of pressure jump IF
 there is a shock. The flatting is needed to eliminated post-shock
 oscillations.
 */
int flatten(int nmin, int nmax, double *pre_1D, double *vx_1D, double *steep, double *flat) {

    int i;
    double temp1, temp2, temp3;
    double delp1, delp2, shock;
    double epsilon, omega1, omega2;

    epsilon = 0.3;
    omega1 = 0.75;
    omega2 = 5.0;

    /*
    Look for presence of a shock using pressure gradient and sign of
    velocity jump:  shock = 1 if there is a shock in the zone, else shock = 0
    Compute steepness parameter based on steepness of pressure jump IF
    there is a shock.
     */

    for (i = nmin - 4; i <= nmax + 4; i++) {
        // delp1 und delp2 are the pressure difference over the next and next after next cells
        delp1 = pre_1D[i + 1] - pre_1D[i - 1];
        delp2 = pre_1D[i + 2] - pre_1D[i - 2];
        
        // the actual shock detection
        if (fabs(delp2) < small) delp2 = small;
        shock = fabs(delp1) / min_flatten(pre_1D[i + 1], pre_1D[i - 1]) - epsilon;
        shock = max_flatten(0.0, shock);
        if (shock > 0.0) shock = 1.0;
        if (vx_1D[i - 1] < vx_1D[i + 1]) shock = 0.0;
        temp1 = (delp1 / delp2 - omega1) * omega2;

        // the steepness parameter
        steep[i] = shock * max_flatten(0.0, temp1);
    }

    //Set phony boundary conditions for the steepness parameter

    steep[nmin - 5] = steep[nmin - 4];
    steep[nmax + 5] = steep[nmax + 4];

    //Set flattening coefficient based on the steepness in neighboring zones
    for (i = nmin - 4; i <= nmax + 4; i++) {
        temp2 = max_flatten(steep[i - 1], steep[i]);
        temp3 = max_flatten(steep[i + 1], temp2);
        flat[i] = max_flatten(0.0, min_flatten(0.5, temp3));
        if (flat[i] > 0) shock_or_not = 1;
    }


    return 0;
}