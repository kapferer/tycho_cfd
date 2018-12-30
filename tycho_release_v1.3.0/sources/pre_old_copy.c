/*
 * pressure_on_solid_reset.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "variables_global.h"
#include "prototypes.h"

/*!
 Here the present pressure is copied into the pre_old
 */
int pre_old_copy(int x, int y, int z) {

    int i, j, k;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, pre_old, pre)
    {
#endif

        for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif    

            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    pre_old[i][j][k] = pre[i][j][k];
                }
            }
        }

#ifdef _OPENMP
    }
#endif

    return 0;
}