/*
 * boundary_velo_corrector.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 example function to suppress boundary streams;
 !!highly dangerous, only for testing!!
 */
int boundary_velo_corrector(int x, int y, int z) {

    int i, j, k;

    for (j = 0; j < y; j++) {
        for (i = 0; i < x; i++) {
            if (wind_direction == 2) {
                vx[0][j][k] = small;
                vx[x - 1][j][k] = small;
            }
        }
    }

    return 0;
}