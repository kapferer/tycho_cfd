/*
 * kinematic_viscosity.c
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
 In order to enable viscosity in an explicit way, this
 routine calculates the kinematic_viscosity for a given temperature
 according to Sutherland's law.
 */
double kinematic_viscosity(double temperature, double density) {

    double kinematic_viscosity;

    if ((temperature >= 200) && (temperature <= 1400)) {
        kinematic_viscosity = (C1_visc * pow(temperature, 1.5)) / (temperature + S_visc);
        kinematic_viscosity = kinematic_viscosity / density;
    } else {
        kinematic_viscosity = 0.0;
    }

    return kinematic_viscosity;
}

/*!
 In order to enable viscosity in an explicit way, this
 routine calculates the dynamic_viscosity for a given temperature
 according to Sutherland's law.
 */
double dynamic_viscosity(double temperature) {

    double dynamic_viscosity;

    if ((temperature >= 200) && (temperature <= 1400)) {
        dynamic_viscosity = (C1_visc * pow(temperature, 1.5)) / (temperature + S_visc);
    } else {
        dynamic_viscosity = 0.0;
    }

    return dynamic_viscosity;
}