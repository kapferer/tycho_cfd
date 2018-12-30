/*
 * advection.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Advect the marker field with the velocities of the hydro-solver
 with a Second-order upwind scheme
 */
int advect(int nmin, int nmax, int flag, double *dx, double *vx_1D, double *marker_1D) {

    int n;
    double spacing;

    if (flag == 0) spacing = 2.0 * xmax / (double) x;
    if (flag == 1) spacing = 2.0 * ymax / (double) y;
    if (flag == 2) spacing = 2.0 * zmax / (double) z;

    //Advect the marker field
    for (n = nmin; n <= nmax; n++) {
        if (vx_1D[n] > 0.0)
            marker_1D[n] = marker_1D[n]-(dt * vx_1D[n]*(3 * marker_1D[n] - 4 * marker_1D[n - 1] + marker_1D[n - 2])) / spacing;
        if (vx_1D[n] < 0.0)
            marker_1D[n] = marker_1D[n]-(dt * vx_1D[n]*(-1 * marker_1D[n + 2] + 4 * marker_1D[n + 1] - 3 * marker_1D[n])) / spacing;
        //to prevent negative marker densities
        if (marker_1D[n] < 0.0) marker_1D[n] = 0.0;
    }

    //Now the ghost cell resions are cleared
    for (n = 1; n <= 6; n++) marker_1D[nmin - n] = 0.0;
    for (n = 1; n <= 6; n++) marker_1D[nmax + n] = 0.0;

    return 0;
}