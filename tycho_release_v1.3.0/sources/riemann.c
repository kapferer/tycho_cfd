/*
 * riemann.c
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
inline double max_riemann(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 The Riemann solver for the shock tube problem.
  
 Input variables are:
 lmin   = zone number of first physical zone
 lmax   = zone number of first ghost zone on right (lmax=nmax+1)
 gamma    = equation of state gamma
 prgh   = pressure state on the right side of the boundary
 plft   = pressure state on the left side of the boundary
 urgh   = velocity state on the right side of the boundary
 ulft   = velocity state on the left side of the boundary
 vrgh   = density state on the right side of the boundary
 rlft   = density state on the left side of the boundary
 (rlft and vrgh are inverted to get the specific volume)

 Output variables are:
 umid   = time_sim averaged velocity at the zone interface
 pmid   = time_sim averaged pressure at the zone interface
 */
int riemann(int nmin, int nmax, double *clft, double *crgh,
        double *rlft, double *rrgh, double *plfti,
        double *prghi, double *pmid, double *pmold,
        double *plft, double *wrgh, double *prgh, double *wlft,
        double *zlft, double *zrgh, double *umidl, double *umidr,
        double *umid, double *urgh, double *ulft) {
    int n, l;

    double tol = 1.0E-3;
    double gamfac2, gamfac1;

    gamfac2 = Gamma1 + 1.0;
    gamfac1 = 0.5 * (gamfac2) / Gamma1;


    //Obtain first guess for Pmid by assuming Wlft, Wrgh = Clft, Crgh
    for (n = nmin; n <= nmax; n++) {
        clft[n] = sqrt(Gamma1 * plft[n] * rlft[n]);
        crgh[n] = sqrt(Gamma1 * prgh[n] * rrgh[n]);
        rlft[n] = 1.0 / rlft[n];
        rrgh[n] = 1.0 / rrgh[n];
        plfti[n] = 1.0 / plft[n];
        prghi[n] = 1.0 / prgh[n];
        pmid[n] = prgh[n] - plft[n] - crgh[n]*(urgh[n] - ulft[n]);
        pmid[n] = plft[n] + pmid[n] * clft[n] / (clft[n] + crgh[n]);
        pmid[n] = max_riemann(smallp, pmid[n]);
    }

    /*
    Iterate up to 8 time_sims using Newton's method to converge on correct Pmid
    -use previous guess for pmid to get wavespeeds: wlft, wrgh
    -find the slope in the u-P plane for each state: zlft, zrgh
    -use the wavespeeds and pmid to guess umid on each side: umidl, umidr
    -project tangents from (pmid,umidl) and (pmid,umidr) to get new pmid
    -make sure pmid does not fall below floor value for pressure
     */
    for (n = nmin; n <= nmax; n++) {
        for (l = 0; l < 12; l++) {
            pmold[n] = pmid[n];
            wlft[n] = 1.0 + gamfac1 * (pmid[n] - plft[n]) * plfti[n];
            wrgh[n] = 1.0 + gamfac1 * (pmid[n] - prgh[n]) * prghi[n];
            wlft[n] = clft[n] * sqrt(wlft[n]);
            wrgh[n] = crgh[n] * sqrt(wrgh[n]);
            zlft[n] = 4.0 * rlft[n] * wlft[n] * wlft[n];
            zrgh[n] = 4.0 * rrgh[n] * wrgh[n] * wrgh[n];
            zlft[n] = -zlft[n] * wlft[n] / (zlft[n] - gamfac2 * (pmid[n] - plft[n]));
            zrgh[n] = zrgh[n] * wrgh[n] / (zrgh[n] - gamfac2 * (pmid[n] - prgh[n]));
            umidl[n] = ulft[n] - (pmid[n] - plft[n]) / wlft[n];
            umidr[n] = urgh[n] + (pmid[n] - prgh[n]) / wrgh[n];
            pmid[n] = pmid[n] + (umidr[n] - umidl[n])*(zlft[n] * zrgh[n]) / (zrgh[n] - zlft[n]);
            pmid[n] = max_riemann(smallp, pmid[n]);
            if (fabs(pmid[n] - pmold[n]) / pmid[n] < tol) break;
        }
    }


    //Calculate umid by averaging umidl, umidr based on new pmid
    for (n = nmin; n <= nmax; n++) {
        umidl[n] = ulft[n] - (pmid[n] - plft[n]) / wlft[n];
        umidr[n] = urgh[n] + (pmid[n] - prgh[n]) / wrgh[n];
        umid [n] = 0.5 * (umidl[n] + umidr[n]);
    }

    return 0;
}