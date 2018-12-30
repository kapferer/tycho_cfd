/*
 * parabola.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//special Headerfiles
#include "variables_global.h"
#include "prototypes.h"

/*!
 Inline function to calculate the signum of a given double
 */
inline double signum_parabola(double a) {
    double tmp = 0.0;

    if (a < 0) tmp = -1.0;
    if (a == 0) tmp = 0.0;
    if (a > 0) tmp = +1.0;

    return tmp;
}

/*!
 Inline function to calculate the minimum of two doubles
 */
inline double min_par(double a, double b) {
    double tmp = 0.0;

    if (a < b) tmp = a;
    if (a == b)tmp = a;
    if (a > b) tmp = b;

    return tmp;

}

/*! 
 Here the picewise parabolic interpolated quantities 
 for the PPM scheme are calculated.
 */
int parabola(int nmin, int nmax, double *a, double *deltaa,
        double *a6, double *al, int flag1, double *diffa,
        double *da, double **para, double *ar, double *flat,
        double *scratch1, double *scratch2, double *scratch3) {
    int i;
    double tmp1, tmp2, tmp3, onemfl;


    for (i = nmin - 2; i <= nmax + 1; i++) {
        diffa[i] = a[i + 1] - a[i];
    }

    //Equation 1.7 of Colella & Woodward, Journal of Computational Physics 53, 174-201. 1984
    //da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))
    for (i = nmin - 1; i <= nmax + 1; i++) {
        da[i] = para[i][3] * diffa[i] + para[i][4] * diffa[i - 1];
        tmp1 = min_par(fabs(da[i]), 2.0 * fabs(diffa[i - 1]));
        tmp2 = min_par(2.0 * fabs(diffa[i]), tmp1);
        tmp3 = signum_parabola(da[i]);
        da[i] = tmp2*tmp3;
    }


    //zero da(n) if a(n) is a local max/min
    for (i = nmin - 1; i <= nmax + 1; i++) {
        if (diffa[i - 1] * diffa[i] < 0.0) da[i] = 0.0;
    }


    //Equation 1.6 of of Colella & Woodward, Journal of Computational Physics 53, 174-201. 1984
    //a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * dma(j+1) + C3 * dma(j)
    //MONOT: Limit ar(n) to the range defined by a(n) and a(n+1)
    for (i = nmin - 1; i <= nmax; i++) {
        ar[i] = a[i] + para[i][0] * diffa[i] + para[i][1] * da[i + 1] + para[i][2] * da[i];
        al[i + 1] = ar[i];
    }

    //eqn. 4.1 - flatten interpolation in zones with a shock ( flat(n)->1. )
    for (i = nmin; i <= nmax; i++) {
        onemfl = 1.0 - flat[i];
        ar[i] = flat[i] * a[i] + onemfl * ar[i];
        al[i] = flat[i] * a[i] + onemfl * al[i];
    }

    /*
     MONOTONICITY constraints:
     compute delta_a, a_6
     MONOT: if a is a local max/min, flatten zone structure ar,al -> a.
     MONOT: compute monotonzsed values using eq. 1.10 of of Colella & Woodward, Journal of Computational Physics 53, 174-201. 1984
     if parabola exceeds al/ar, reset ar/al so that slope -> 0.
     Recalculate delta_a and a_6
     */
    for (i = nmin; i <= nmax; i++) {

        onemfl = 1.0 - flat[i];
        ar[i] = flat[i] * a[i] + onemfl * ar[i];
        al[i] = flat[i] * a[i] + onemfl * al[i];

        deltaa[i] = ar[i] - al[i];
        a6[i] = 6.0 * (a[i] - 0.5 * (al[i] + ar[i]));
        scratch1[i] = (ar[i] - a[i])*(a[i] - al[i]);
        scratch2[i] = deltaa[i] * deltaa[i];
        scratch3[i] = deltaa[i] * a6[i];

        if (scratch1[i] <= 0.0) {
            ar[i] = a[i];
            al[i] = a[i];
        }
        if (scratch2[i] < +scratch3[i]) al[i] = 3.0 * a[i] - 2.0 * ar[i];
        if (scratch2[i] < -scratch3[i]) ar[i] = 3.0 * a[i] - 2.0 * al[i];


        deltaa[i] = ar[i] - al[i];
        a6[i] = 6.0 * (a[i] - 0.5 * (al[i] + ar[i]));

    }

    for (i = nmin; i <= nmax; i++) {
        if (scratch1[i] <= 0.0) {
            ar[i] = a[i];
            al[i] = a[i];
        }
        if (scratch2[i] < +scratch3[i]) al[i] = 3.0 * a[i] - 2.0 * ar[i];
        if (scratch2[i] < -scratch3[i]) ar[i] = 3.0 * a[i] - 2.0 * al[i];
    }



    for (i = nmin; i <= nmax; i++) {
        deltaa[i] = ar[i] - al[i];
        a6[i] = 6.0 * (a[i] - 0.5 * (al[i] + ar[i]));
    }

    return 0;
}
