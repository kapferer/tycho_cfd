/*
 * hydro_sweeps.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 * The forward and backward hydro sweeps the direction
 * gives the forward [0] or backward [1] sweeps
 */
int hydro_sweeps(int x, int y, int z, int direction) {
    // Here we switch the wind on
    if (wind_on_off == 1) wind(x, y, z);

    if (direction == 0) {
        if (viscosity_on_off == 1) copy_arrays_for_viscosity(x, y, z);
        sweep_x(x, y, z, 0);
        if (viscosity_on_off == 1) copy_arrays_for_viscosity(x, y, z);
        if (dimension > 1) sweep_y(x, y, z, 1);
        if (viscosity_on_off == 1) copy_arrays_for_viscosity(x, y, z);
        if (dimension == 3) sweep_z(x, y, z, 2);
    }
    
    if (direction == 1) {
        if (viscosity_on_off == 1) copy_arrays_for_viscosity(x, y, z);
        if (dimension == 3) sweep_z(x, y, z, 2);
        if (viscosity_on_off == 1) copy_arrays_for_viscosity(x, y, z);
        if (dimension > 1) sweep_y(x, y, z, 1);
        if (viscosity_on_off == 1) copy_arrays_for_viscosity(x, y, z);
        sweep_x(x, y, z, 0);
    }
    
    return 0;
}
