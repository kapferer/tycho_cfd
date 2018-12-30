/*
 * remap.c
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
inline double max_remap(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 Remap the new hydro-quantities on the computational domain grid.
 */
int remap(int nmin, int nmax, int flag, double *a_coef, double *dx,
        double *ai_coef, double *b_coef, double *bi_coef, double *c_coef,
        double *d_x, double **para, double *ci_coef, double *dr, double *r6,
        double *rl, double *diffa, double *da, double *ar, double *flat,
        double *scratch1, double *scratch2, double *scratch3, double *du,
        double *u6, double *ul, double *dv, double *v6, double *vl, double *w6,
        double *wl, double *dq, double *q6, double *ql, double *de, double *e6,
        double *el, double *xa, double *xa0, double *delta, double *fluxr,
        double *fluxu, double *fluxv, double *fluxw, double *dw, double *fluxe,
        double *fluxq, double *dm, double *rho_1D, double *dvol, double *dvol0,
        double *dm0, double * vx_1D, double *vy_1D, double *vz_1D, double *eng_1D,
        double *e_int_1D, double *pre_1D) {

    int n, nn;

    double fractn, fractn2, ekin;
    double deltx;
    double fourthd = 4.0 / 3.0;

    //Generate interpolation functions, saving da, al for
    //constructing left and right total energy states.
    para_coef(nmin - 1, nmax + 1, 1, a_coef, dx,
            ai_coef, b_coef, bi_coef, c_coef,
            d_x, para, ci_coef);
    parabola(nmin - 1, nmax + 1, rho_1D, dr, r6, rl, 3, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 1, nmax + 1, vx_1D, du, u6, ul, 4, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 1, nmax + 1, vy_1D, dv, v6, vl, 5, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 1, nmax + 1, vz_1D, dw, w6, wl, 6, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 1, nmax + 1, e_int_1D, dq, q6, ql, 7, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 1, nmax + 1, eng_1D, de, e6, el, 8, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);

    for (n = nmin; n <= nmax + 1; n++) {
        delta[n] = xa[n] - xa0[n];
    }


    //Calculate the total mass (fluxr), momentum (fluxu), and energy (fluxe)
    //in the subshell created by the overlap of the Lagrangian and Eulerian grids.
    //If the zone face has moved to the left (deltx > 0), use the integral from the
    //left side of zone n (fluxrr).  If the zone face has moved to the right
    //(deltx < 0), use the integral from the right side of zone nn=n-1 (fluxrl).
    for (n = nmin; n <= nmax + 1; n++) {
        deltx = xa[n] - xa0[n];
        if (deltx >= 0.0) {
            nn = n - 1;
            fractn = 0.5 * deltx / dx[nn];
            fractn2 = 1.0 - fourthd*fractn;
            fluxr[n] = (rl[nn] + dr[nn] - fractn * (dr[nn] - fractn2 * r6[nn])) * delta[n];
            fluxu[n] = (ul[nn] + du[nn] - fractn * (du[nn] - fractn2 * u6[nn])) * fluxr[n];
            fluxv[n] = (vl[nn] + dv[nn] - fractn * (dv[nn] - fractn2 * v6[nn])) * fluxr[n];
            fluxw[n] = (wl[nn] + dw[nn] - fractn * (dw[nn] - fractn2 * w6[nn])) * fluxr[n];
            fluxe[n] = (el[nn] + de[nn] - fractn * (de[nn] - fractn2 * e6[nn])) * fluxr[n];
            fluxq[n] = (ql[nn] + dq[nn] - fractn * (dq[nn] - fractn2 * q6[nn])) * fluxr[n];
        } else {
            fractn = 0.5 * deltx / dx[n];
            fractn2 = 1.0 + fourthd*fractn;
            fluxr[n] = (rl[n] - fractn * (dr[n] + fractn2 * r6[n])) * delta[n];
            fluxu[n] = (ul[n] - fractn * (du[n] + fractn2 * u6[n])) * fluxr[n];
            fluxv[n] = (vl[n] - fractn * (dv[n] + fractn2 * v6[n])) * fluxr[n];
            fluxw[n] = (wl[n] - fractn * (dw[n] + fractn2 * w6[n])) * fluxr[n];
            fluxe[n] = (el[n] - fractn * (de[n] + fractn2 * e6[n])) * fluxr[n];
            fluxq[n] = (ql[n] - fractn * (dq[n] + fractn2 * q6[n])) * fluxr[n];
        }
    }

    //Advect mass, momentum, and energy by moving the sub-shell quantities
    //into the appropriate Eulerian zone.
    for (n = nmin; n <= nmax; n++) {
        dm[n] = rho_1D[n] * dvol[n];
        dm0[n] = dm[n] + fluxr[n] - fluxr[n + 1];
        rho_1D[n] = dm0[n] / dvol0[n];
        rho_1D[n] = max_remap(smallr, rho_1D[n]);
        dm0[n] = 1.0 / (rho_1D[n] * dvol0[n]);
        vx_1D[n] = (vx_1D[n] * dm[n] + fluxu[n] - fluxu[n + 1]) * dm0[n];
        vy_1D[n] = (vy_1D[n] * dm[n] + fluxv[n] - fluxv[n + 1]) * dm0[n];
        vz_1D[n] = (vz_1D[n] * dm[n] + fluxw[n] - fluxw[n + 1]) * dm0[n];
        eng_1D[n] = (eng_1D[n] * dm[n] + fluxe[n] - fluxe[n + 1]) * dm0[n];
        e_int_1D[n] = (e_int_1D[n] * dm[n] + fluxq[n] - fluxq[n + 1]) * dm0[n];
        //If flow is highly supersonic remap on internal energy, else on total energy
        ekin = 0.5 * (pow(vx_1D[n], 2) + pow(vy_1D[n], 2) + pow(vz_1D[n], 2));
        if (ekin / e_int_1D[n] < 100.0) e_int_1D[n] = eng_1D[n] - ekin;
        pre_1D[n] = Gamma * rho_1D[n] * e_int_1D[n];
        pre_1D[n] = max_remap(smallp, pre_1D[n]);
    }

    return 0;
}