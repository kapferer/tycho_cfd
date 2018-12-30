/*
 * volume.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>

#include "prototypes.h"
#include "variables_global.h"

/*
 Fills the dvol, dvol0 arrays
 */
int volume(int nmin, int nmax, double *dvol, double *dx, double *dvol0, double *dx0) {
    int i;

    for (i = nmin - 3; i <= nmax + 4; i++) {
        dvol[i] = dx[i];
        dvol0[i] = dx0[i];
    }

    return 0;
}