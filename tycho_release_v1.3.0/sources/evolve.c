/*
 * evolve.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Inline function to calculate the maximum of two doubles
 */
inline double max_evolve(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 Evolution of velocities and energies due to 
 pressure acceleration and external
 forces.
 */
int evolve(int nmin, int nmax, int flag, double *rho_1D, double *dvol,
        double *dm, double *dtbdm, double *xa1, double *xa,
        double *dvol1, double *umid, double *upmid,
        double *pmid, double *xa2, double *dx, double *xa3,
        double *vx_1D_old, double *vx_1D, double *vy_1D,
        double *vz_1D, double *eng_1D,
        double *e_int_1D, double *pre_1D) {

    int i;

    // Calculate the mass and upmid in a cell
    for (i = nmin - 3; i <= nmax + 4; i++) {
        dm[i] = rho_1D[i] * dvol[i];
        dtbdm[i] = dt / dm[i];
        xa1[i] = xa[i];
        dvol1[i] = dvol[i];
        xa[i] = xa[i] + dt * umid[i];
        upmid[i] = umid[i] * pmid[i];
    }

    xa1[nmin - 4] = xa[nmin - 4];
    xa1[nmax + 5] = xa[nmax + 5];

    // Calculate the cell centered coordinates
    for (i = nmin - 4; i <= nmax + 5; i++) {
        xa2[i] = xa1[i] + 0.5 * dx[i];
        dx[i] = xa[i + 1] - xa[i];
        xa3[i] = xa[i] + 0.5 * dx[i];
    }

    for (i = nmin - 3; i <= nmax + 4; i++) {
        dvol[i] = dx[i];
    }

    // new densities
    for (i = nmin - 3; i <= nmax + 3; i++) {
        rho_1D[i] = rho_1D[i]*(dvol1[i] / dvol[i]);
        rho_1D[i] = max_evolve(rho_1D[i], smallr);

        // velocity evolution due to pressure acceleration and external forces.
        vx_1D_old[i] = vx_1D[i];
        if (flag == 1) {
            if (gravity_on_off == 1) {
                vx_1D[i] = vx_1D[i] - dtbdm[i]*(pmid[i + 1] - pmid[i]) + dt * grav_acc;
                eng_1D[i] = eng_1D[i] - dtbdm[i]*(upmid[i + 1] - upmid[i]) + 0.5 * dt *
                        grav_acc * (vx_1D_old[i] + vx_1D[i]);
            }
            if (gravity_on_off == 0) {
                vx_1D[i] = vx_1D[i] - dtbdm[i]*(pmid[i + 1] - pmid[i]);
                eng_1D[i] = eng_1D[i] - dtbdm[i]*(upmid[i + 1] - upmid[i]);
            }
        } else {
            vx_1D[i] = vx_1D[i] - dtbdm[i]*(pmid[i + 1] - pmid[i]);
            eng_1D[i] = eng_1D[i] - dtbdm[i]*(upmid[i + 1] - upmid[i]);
        }

        // total energy evolution
        e_int_1D[i] = eng_1D[i] - 0.5 * (pow(vx_1D[i], 2) + pow(vy_1D[i], 2) + pow(vz_1D[i], 2));
        pre_1D[i] = max_evolve(rho_1D[i] * e_int_1D[i] * Gamma, smallp);
    }

    return 0;
}