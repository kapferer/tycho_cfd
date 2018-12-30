/*
 * copy arrays for viscosity.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 Copies the density and velocity arrays before each hydro sweep for the 
 viscosity calculation.
 */


int copy_arrays_for_viscosity(int x, int y, int z) {

    int i, j, k;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x,y,z, rho_visc, \
                vx_visc, vy_visc, vz_visc, \
                rho, vx, vy, vz)
    {
#endif

        for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif  
            for (j = 0; j < y; j++) {
                
                for (k = 0; k < z; k++) {
                    rho_visc[i][j][k] = rho[i][j][k];
                    vx_visc[i][j][k] = vx[i][j][k];
                    vy_visc[i][j][k] = vy[i][j][k];
                    vz_visc[i][j][k] = vz[i][j][k];
                }
            }
        }

#ifdef _OPENMP        
    }
#endif 

    return 0;
}