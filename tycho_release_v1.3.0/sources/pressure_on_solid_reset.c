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
 Here the pressure on solid arrays
 pressure_on_solid,  , 
   are set to zero.
 */
int pressure_on_solid_reset(int x, int y, int z) {

    int i, j, k;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, pressure_on_solid)
    {
#endif

        for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif    

            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    pressure_on_solid[i][j][k] = 0.0;
                }
            }
        }

#ifdef _OPENMP
    }
#endif

    return 0;
}