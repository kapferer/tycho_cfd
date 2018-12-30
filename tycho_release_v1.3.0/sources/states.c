/*
 * states.c
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
inline double max_states(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 Get hydroquantities at cell centres and faces.
 */
int states(int nmin, int nmax, int flag, double *pre_1D, double *rho_1D,
        double *dx, double *Cdtdx, double *fCdtdx, double *plft,
        double *pl, double *dp, double *p6, double *ulft, double *ul,
        double *du, double *u6, double *rlft, double *rl, double *dr,
        double *r6, double *prgh, double *urgh, double *rrgh) {
    int i, ii;
    double hdt, svel;
    double fourthd;

    fourthd = 4.0 / 3.0;
    svel = 0.0;
    hdt = 0.5 * dt;

    for (i = nmin - 4; i <= nmax + 4; i++) {
        Cdtdx[i] = sqrt(Gamma1 * pre_1D[i] / rho_1D[i]) / (dx[i]);
        svel = max_states(svel, Cdtdx[i]);
        Cdtdx[i] = Cdtdx[i] * hdt;
        fCdtdx[i] = 1.0 - fourthd * Cdtdx[i];
    }

    for (i = nmin - 4; i <= nmax + 4; i++) {
        ii = i + 1;
        plft[ii] = pl[i] + dp[i] - Cdtdx[i]*(dp[i] - fCdtdx[i] * p6[i]);
        ulft[ii] = ul[i] + du[i] - Cdtdx[i]*(du[i] - fCdtdx[i] * u6[i]);
        rlft[ii] = rl[i] + dr[i] - Cdtdx[i]*(dr[i] - fCdtdx[i] * r6[i]);
        plft[ii] = max_states(smallp, plft[ii]);
        rlft[ii] = max_states(smallr, rlft[ii]);
        if (flag == 0) {
            ulft[ii] = ulft[ii] + hdt * grav_acc;
        } else {
            ulft[ii] = ulft[ii];
        }
        prgh[i] = pl[i] + Cdtdx[i]*(dp[i] + fCdtdx[i] * p6[i]);
        urgh[i] = ul[i] + Cdtdx[i]*(du[i] + fCdtdx[i] * u6[i]);
        rrgh[i] = rl[i] + Cdtdx[i]*(dr[i] + fCdtdx[i] * r6[i]);
        prgh[i] = max_states(smallp, prgh[i]);
        rrgh[i] = max_states(smallr, rrgh[i]);
        if (flag == 1) {
            if (gravity_on_off == 1) urgh[i] = urgh[i] + hdt * grav_acc;
            if (gravity_on_off == 0) urgh[i] = urgh[i];
        } else {
            urgh[i] = urgh[i];
        }
    }
    
  
    return 0;
}