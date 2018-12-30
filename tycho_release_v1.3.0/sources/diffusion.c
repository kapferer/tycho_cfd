/*
 * diffusion.c
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
 In order to enable diffusion in an explicit way, this
 routine solves the diffusion equation in a very
 basic way using operator splitting. Use it with caution!
 This is only done within the computational domain of
 the obstacles.
 */
int diffusion(int x, int y, int z) {

    int i, j, k, n;
    double dxfac, dyfac, dzfac;
    double coefficient_solid;
    double alpha_solid_x, alpha_solid_y, alpha_solid_z;
    double *temperature_1D, *temperature_1D_future;
    double *mass;
    double energy_left, energy_right;
    double temperature_right_gas, temerpature_left_gas;
    double A;
    double heat_transfer_coefficient;
    double volume;
    int *dom_1D;
    int dom_counter;
    int i_local_min, i_local_max;
    int condition;

    dom_counter = 0;
    i_local_min = 0;
    i_local_max = 0;
    temperature_right_gas = 0.0;
    temerpature_left_gas = 0.0;
    A = 0.0;
    heat_transfer_coefficient = 0.0;
    condition = 0;

    // spacing in all three directions
    dxfac = (xmax - xmin) / (double) x;
    dyfac = (ymax - ymin) / (double) y;
    dzfac = (zmax - zmin) / (double) z;
    if (dimension == 1) volume = dxfac;
    if (dimension == 2) volume = dxfac * dyfac;
    if (dimension == 3) volume = dxfac * dyfac * dzfac;


    dxfac = dxfac*dxfac;
    dyfac = dyfac*dyfac;
    dzfac = dzfac*dzfac;

    coefficient_solid = obstacle_heat_conductivity;

    if (dimension > 0) alpha_solid_x = coefficient_solid * dt / (dxfac);
    if (dimension > 1) alpha_solid_y = coefficient_solid * dt / (dyfac);
    if (dimension > 2) alpha_solid_z = coefficient_solid * dt / (dzfac);


    if ((alpha_solid_x >= 0.3) && (dimension > 0)) {
        dt = (0.3 * dxfac) / coefficient_solid;
    }
    if ((alpha_solid_y >= 0.3) && (dimension > 1)) {
        dt = (0.3 * dyfac) / coefficient_solid;
    }
    if ((alpha_solid_z >= 0.3) && (dimension > 2)) {
        dt = (0.3 * dzfac) / coefficient_solid;
    }


    //THE X-SWEEP=========================================================================================
#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, n,  temperature_1D, dom_1D, \
                        temperature_1D_future, mass, dom_counter, \
                        i_local_min, i_local_max, energy_left, energy_right, \
                        condition, heat_transfer_coefficient) \
                shared(x, y, z, dom, rho, pre, dt, dxfac, \
                       dyfac, dzfac, alpha_solid_x, gasconstant, A, \
                       dimension, volume, \
                       specific_heat_capacity_gas, \
                       specific_heat_capacity_obstacle)
    {
#endif

        init_diffusion(&temperature_1D, &temperature_1D_future, &mass, &dom_1D);

        if (dimension == 1) A = 1;
        if (dimension == 2) A = dyfac;
        if (dimension == 3) A = dyfac * dzfac;


        //the x-sweep
        for (k = 0; k < z; k++) {


#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif      

            for (j = 0; j < y; j++) {

                condition = 0;
                i_local_min = 0;
                i_local_max = 0;

                while (condition == 0) {

                    dom_counter = 0;
                    energy_left = 0.0;
                    energy_right = 0.0;

                    for (i = i_local_min; i < x; i++) {

                        temperature_1D[i] = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                        mass[i] = rho[i][j][k] * volume;
                        dom_1D[i] = dom[i][j][k];

                        if (dom_1D[i] == 1) {
                            //Now we have the beginning- and the end-index of the solid rod
                            if (dom_counter == 0) i_local_min = i;
                            dom_counter++;
                        } else {
                            if (dom_counter > 0) {
                                i_local_max = i - 1;
                                break;
                            }
                        }
                        if (i == (x - 1)) {
                            i_local_max = i;
                            if (dom_counter == 0) i_local_min = i;
                        }
                    }

                    if (i == x) {
                        condition = 1;
                    }

                    //in the boundary of the obstacle
                    //case if i_local_min > 0
                    if (i_local_min != 0) {
                        heat_transfer_coefficient = alpha_heat_transfer(i_local_min - 1, j, k, 0);
                        energy_left = heat_transfer_coefficient * A * (temperature_1D[i_local_min - 1] - temperature_1D[i_local_min]) * dt;
                    }
                    if (i_local_max != (x - 1)) {
                        heat_transfer_coefficient = alpha_heat_transfer(i_local_max + 1, j, k, 0);
                        energy_right = heat_transfer_coefficient * A * (temperature_1D[i_local_max + 1] - temperature_1D[i_local_max]) * dt;
                    }

                    //temperature in the obstacle higher than in the gas left of it
                    if ((energy_left <= 0) && (i_local_min != 0)) {
                        temperature_1D[i_local_min - 1] = temperature_1D[i_local_min - 1] + fabs(energy_left / (specific_heat_capacity_gas * mass[i_local_min - 1]));
                        temperature_1D[i_local_min] = temperature_1D[i_local_min] - fabs(energy_left / (specific_heat_capacity_obstacle * mass[i_local_min]));

                        pre[i_local_min - 1][j][k] = (gasconstant * rho[i_local_min - 1][j][k]) * temperature_1D[i_local_min - 1];
                        pre[i_local_min][j][k] = (gasconstant * rho[i_local_min][j][k]) * temperature_1D[i_local_min];

                    }

                    //temperature in the obstacle lower than in the gas left of it
                    if ((energy_left > 0) && (i_local_min != 0)) {
                        temperature_1D[i_local_min - 1] = temperature_1D[i_local_min - 1] - fabs(energy_left / (specific_heat_capacity_gas * mass[i_local_min - 1]));
                        temperature_1D[i_local_min] = temperature_1D[i_local_min] + fabs(energy_left / (specific_heat_capacity_obstacle * mass[i_local_min]));

                        pre[i_local_min - 1][j][k] = (gasconstant * rho[i_local_min - 1][j][k]) * temperature_1D[i_local_min - 1];
                        pre[i_local_min][j][k] = (gasconstant * rho[i_local_min][j][k]) * temperature_1D[i_local_min];

                    }

                    //temperature in the gas is higher than in the gas right of it
                    if ((energy_right > 0) && (i_local_max != (x - 1))) {
                        temperature_1D[i_local_max] = temperature_1D[i_local_max] + fabs(energy_right / (specific_heat_capacity_obstacle * mass[i_local_max]));
                        temperature_1D[i_local_max + 1] = temperature_1D[i_local_max + 1] - fabs(energy_right / (specific_heat_capacity_gas * mass[i_local_max + 1]));

                        pre[i_local_max][j][k] = (gasconstant * rho[i_local_max][j][k]) * temperature_1D[i_local_max];
                        pre[i_local_max + 1][j][k] = (gasconstant * rho[i_local_max + 1][j][k]) * temperature_1D[i_local_max + 1];

                    }

                    //temperature in the gas lower than in the obstacle left of it
                    if ((energy_right <= 0) && (i_local_max != (x - 1))) {
                        temperature_1D[i_local_max] = temperature_1D[i_local_max] - fabs(energy_right / (specific_heat_capacity_obstacle * mass[i_local_max]));
                        temperature_1D[i_local_max + 1] = temperature_1D[i_local_max + 1] + fabs(energy_right / (specific_heat_capacity_gas * mass[i_local_max + 1]));

                        pre[i_local_max][j][k] = (gasconstant * rho[i_local_max][j][k]) * temperature_1D[i_local_max];
                        pre[i_local_max + 1][j][k] = (gasconstant * rho[i_local_max + 1][j][k]) * temperature_1D[i_local_max + 1];

                    }

                    if (dom_counter > 2) {
                        for (n = i_local_min + 1; n < i_local_max; n++) {
                            temperature_1D_future[n] = alpha_solid_x * temperature_1D[n + 1] +
                                    (1 - 2 * alpha_solid_x) * temperature_1D[n] +
                                    alpha_solid_x * temperature_1D[n - 1];
                        }

                        for (i = i_local_min + 1; i < i_local_max; i++) {
                            pre[i][j][k] = (gasconstant * rho[i][j][k]) * temperature_1D_future[i];
                        }
                    }

                    for (i = i_local_min; i < i_local_max; i++) {
                        temperature_1D[i] = 0.0;
                        temperature_1D_future[i] = 0.0;
                        mass[i] = 0.0;
                        dom_1D[i] = 0;
                    }

                    i_local_min = i_local_max + 1;
                }
            }
        }
        free_diffusion(&temperature_1D, &temperature_1D_future, &mass, &dom_1D);

#ifdef _OPENMP        
    }
#endif    


    //THE Y-SWEEP=========================================================================================
#ifdef _OPENMP      
#pragma omp parallel default(none) \
                        private(j, i, k, n,  temperature_1D, dom_1D, \
                                temperature_1D_future, mass, dom_counter, \
                                i_local_min, i_local_max, energy_left, \
                                condition, energy_right, \
                                heat_transfer_coefficient) \
                        shared(x, y, z, dom, rho, pre, dt, dxfac, \
                               dyfac, dzfac, alpha_solid_y, \
                               gasconstant, A, dimension, volume, \
                               specific_heat_capacity_gas, \
                               specific_heat_capacity_obstacle)
    {
#endif

        init_diffusion(&temperature_1D, &temperature_1D_future, &mass, &dom_1D);

        if (dimension == 1) A = 1;
        if (dimension == 2) A = dxfac;
        if (dimension == 3) A = dxfac * dzfac;

        //the y-sweep
        for (k = 0; k < z; k++) {

#ifdef _OPENMP  
#pragma omp for schedule(static)
#endif

            for (i = 0; i < x; i++) {

                condition = 0;
                i_local_min = 0;
                i_local_max = 0;

                while (condition == 0) {

                    dom_counter = 0;
                    energy_left = 0.0;
                    energy_right = 0.0;


                    for (j = i_local_min; j < y; j++) {

                        temperature_1D[j] = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                        mass[j] = rho[i][j][k] * volume;
                        dom_1D[j] = dom[i][j][k];

                        if (dom_1D[j] == 1) {
                            //Now we have the beginning- and the end-index of the solid rod
                            if (dom_counter == 0) i_local_min = j;
                            dom_counter++;
                        } else {
                            if (dom_counter > 0) {
                                i_local_max = j - 1;
                                break;
                            }
                        }
                        if (j == (y - 1)) {
                            i_local_max = j;
                            if (dom_counter == 0) i_local_min = j;
                        }
                    }

                    if (j == y) {
                        condition = 1;
                    }

                    //Now the heat transfer is calcuated to set the the right temerpature
                    //in the boundary of the obstacle
                    //case if i_local_min > 0
                    if (i_local_min != 0) {
                        heat_transfer_coefficient = alpha_heat_transfer(i, i_local_min - 1, k, 1);
                        energy_left = heat_transfer_coefficient * A * (temperature_1D[i_local_min - 1] - temperature_1D[i_local_min]) * dt;
                    }
                    if (i_local_max != (y - 1)) {
                        heat_transfer_coefficient = alpha_heat_transfer(i, i_local_max + 1, k, 1);
                        energy_right = heat_transfer_coefficient * A * (temperature_1D[i_local_max + 1] - temperature_1D[i_local_max]) * dt;
                    }

                    //temperature in the obstacle higher than in the gas left of it
                    if ((energy_left <= 0) && (i_local_min != 0)) {
                        temperature_1D[i_local_min - 1] = temperature_1D[i_local_min - 1] + fabs(energy_left / (specific_heat_capacity_gas * mass[i_local_min - 1]));
                        temperature_1D[i_local_min] = temperature_1D[i_local_min] - fabs(energy_left / (specific_heat_capacity_obstacle * mass[i_local_min]));

                        pre[i][i_local_min - 1][k] = (gasconstant * rho[i][i_local_min - 1][k]) * temperature_1D[i_local_min - 1];
                        pre[i][i_local_min][k] = (gasconstant * rho[i][i_local_min][k]) * temperature_1D[i_local_min];
                    }

                    //temperature in the obstacle lower than in the gas left of it
                    if ((energy_left > 0) && (i_local_min != 0)) {
                        temperature_1D[i_local_min - 1] = temperature_1D[i_local_min - 1] - fabs(energy_left / (specific_heat_capacity_gas * mass[i_local_min - 1]));
                        temperature_1D[i_local_min] = temperature_1D[i_local_min] + fabs(energy_left / (specific_heat_capacity_obstacle * mass[i_local_min]));

                        pre[i][i_local_min - 1][k] = (gasconstant * rho[i][i_local_min - 1][k]) * temperature_1D[i_local_min - 1];
                        pre[i][i_local_min][k] = (gasconstant * rho[i][i_local_min][k]) * temperature_1D[i_local_min];
                    }

                    //temperature in the obstacle higher than in the gas right of it
                    if ((energy_right > 0) && (i_local_max != (y - 1))) {
                        temperature_1D[i_local_max] = temperature_1D[i_local_max] + fabs(energy_right / (specific_heat_capacity_obstacle * mass[i_local_max]));
                        temperature_1D[i_local_max + 1] = temperature_1D[i_local_max + 1] - fabs(energy_right / (specific_heat_capacity_gas * mass[i_local_max + 1]));

                        pre[i][i_local_max][k] = (gasconstant * rho[i][i_local_max][k]) * temperature_1D[i_local_max];
                        pre[i][i_local_max + 1][k] = (gasconstant * rho[i][i_local_max + 1][k]) * temperature_1D[i_local_max + 1];
                    }

                    //temperature in the obstacle lower than in the gas right of it
                    if ((energy_right <= 0) && (i_local_max != (y - 1))) {
                        temperature_1D[i_local_max] = temperature_1D[i_local_max] - fabs(energy_right / (specific_heat_capacity_obstacle * mass[i_local_max]));
                        temperature_1D[i_local_max + 1] = temperature_1D[i_local_max + 1] + fabs(energy_right / (specific_heat_capacity_gas * mass[i_local_max + 1]));

                        pre[i][i_local_max][k] = (gasconstant * rho[i][i_local_max][k]) * temperature_1D[i_local_max];
                        pre[i][i_local_max + 1][k] = (gasconstant * rho[i][i_local_max + 1][k]) * temperature_1D[i_local_max + 1];
                    }

                    //printf("energy_left: %g energy_right: %g dt: %g\n", energy_left, energy_right, dt);

                    if (dom_counter > 2) {

                        for (n = i_local_min + 1; n < i_local_max; n++) {
                            //here the heat conductivity constant are determined
                            temperature_1D_future[n] = alpha_solid_y * temperature_1D[n + 1] +
                                    (1 - 2 * alpha_solid_y) * temperature_1D[n] +
                                    alpha_solid_y * temperature_1D[n - 1];
                        }

                        for (j = i_local_min + 1; j < i_local_max; j++) {
                            pre[i][j][k] = (gasconstant * rho[i][j][k]) * temperature_1D_future[j];
                        }


                        for (j = i_local_min; j < i_local_max; j++) {
                            temperature_1D[j] = 0.0;
                            temperature_1D_future[j] = 0.0;
                            mass[j] = 0.0;
                            dom_1D[j] = 0;
                        }
                    }

                    i_local_min = i_local_max + 1;

                }
            }
        }
        free_diffusion(&temperature_1D, &temperature_1D_future, &mass, &dom_1D);

#ifdef _OPENMP         
    }
#endif    

    //THE Z-SWEEP=========================================================================================
#ifdef _OPENMP  
#pragma omp parallel default(none) \
                        private(j, i, k, n,  temperature_1D, dom_1D, \
                                temperature_1D_future, mass, dom_counter, \
                                i_local_min, i_local_max, energy_left, \
                                condition, energy_right, \
                                heat_transfer_coefficient) \
                        shared(x, y, z, dom, rho, pre, dt, dxfac, \
                               dyfac, dzfac, alpha_solid_z, gasconstant, A, \
                               dimension, volume, \
                               specific_heat_capacity_gas, \
                               specific_heat_capacity_obstacle)
    {
#endif  
        init_diffusion(&temperature_1D, &temperature_1D_future, &mass, &dom_1D);

        if (dimension == 1) A = 1;
        if (dimension == 2) A = dxfac;
        if (dimension == 3) A = dxfac * dyfac;

        //the z-sweep
        for (j = 0; j < y; j++) {

#ifdef _OPENMP  
#pragma omp for schedule(static)
#endif 

            for (i = 0; i < x; i++) {

                condition = 0;
                i_local_min = 0;
                i_local_max = 0;

                while (condition == 0) {

                    dom_counter = 0;
                    energy_left = 0.0;
                    energy_right = 0.0;

                    for (k = i_local_min; k < z; k++) {

                        temperature_1D[k] = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                        mass[k] = rho[i][j][k] * volume;
                        dom_1D[k] = dom[i][j][k];

                        if (dom_1D[k] == 1) {
                            //Now we have the beginning- and the end-index of the solid rod
                            if (dom_counter == 0) i_local_min = k;
                            dom_counter++;
                        } else {
                            if (dom_counter > 0) {
                                i_local_max = k - 1;
                                break;
                            }
                        }
                        if (k == (z - 1)) {
                            i_local_max = k;
                            if (dom_counter == 0) i_local_min = k;
                        }
                    }

                    if (k == z) {
                        condition = 1;
                    }

                    //Now the heat transfer is calculated to set the the right temperature
                    //in the boundary of the obstacle
                    //case if i_local_min > 0
                    if (i_local_min != 0) {
                        heat_transfer_coefficient = alpha_heat_transfer(i, j, i_local_min - 1, 2);
                        energy_left = heat_transfer_coefficient * A * (temperature_1D[i_local_min - 1] - temperature_1D[i_local_min]) * dt;
                    }
                    if (i_local_max != (z - 1)) {
                        heat_transfer_coefficient = alpha_heat_transfer(i, j, i_local_max + 1, 2);
                        energy_right = heat_transfer_coefficient * A * (temperature_1D[i_local_max + 1] - temperature_1D[i_local_max]) * dt;
                    }

                    //temperature in the obstacle higher than in the gas left of it
                    if ((energy_left <= 0) && (i_local_min != 0)) {
                        temperature_1D[i_local_min - 1] = temperature_1D[i_local_min - 1] + fabs(energy_left / (specific_heat_capacity_gas * mass[i_local_min - 1]));
                        temperature_1D[i_local_min] = temperature_1D[i_local_min] - fabs(energy_left / (specific_heat_capacity_obstacle * mass[i_local_min]));

                        pre[i][j][i_local_min - 1] = (gasconstant * rho[i][j][i_local_min - 1]) * temperature_1D[i_local_min - 1];
                        pre[i][j][i_local_min] = (gasconstant * rho[i][j][i_local_min]) * temperature_1D[i_local_min ];
                    }

                    //temperature in the obstacle lower than in the gas left of it
                    if ((energy_left > 0) && (i_local_min != 0)) {
                        temperature_1D[i_local_min - 1] = temperature_1D[i_local_min - 1] - fabs(energy_left / (specific_heat_capacity_gas * mass[i_local_min - 1]));
                        temperature_1D[i_local_min] = temperature_1D[i_local_min] + fabs(energy_left / (specific_heat_capacity_obstacle * mass[i_local_min]));

                        pre[i][j][i_local_min - 1] = (gasconstant * rho[i][j][i_local_min - 1]) * temperature_1D[i_local_min - 1];
                        pre[i][j][i_local_min] = (gasconstant * rho[i][j][i_local_min]) * temperature_1D[i_local_min ];
                    }

                    //temperature in the obstacle higher than in the gas right of it
                    if ((energy_right > 0) && (i_local_max != (z - 1))) {
                        temperature_1D[i_local_max] = temperature_1D[i_local_max] + fabs(energy_right / (specific_heat_capacity_obstacle * mass[i_local_max]));
                        temperature_1D[i_local_max + 1] = temperature_1D[i_local_max + 1] - fabs(energy_right / (specific_heat_capacity_gas * mass[i_local_max + 1]));

                        pre[i][j][i_local_max] = (gasconstant * rho[i][j][i_local_max]) * temperature_1D[i_local_max ];
                        pre[i][j][i_local_max + 1] = (gasconstant * rho[i][j][i_local_max + 1]) * temperature_1D[i_local_max + 1];
                    }

                    //temperature in the obstacle lower than in the gas right of it
                    if ((energy_right <= 0) && (i_local_max != (z - 1))) {
                        temperature_1D[i_local_max] = temperature_1D[i_local_max] - fabs(energy_right / (specific_heat_capacity_obstacle * mass[i_local_max]));
                        temperature_1D[i_local_max + 1] = temperature_1D[i_local_max + 1] + fabs(energy_right / (specific_heat_capacity_gas * mass[i_local_max + 1]));

                        pre[i][j][i_local_max] = (gasconstant * rho[i][j][i_local_max]) * temperature_1D[i_local_max ];
                        pre[i][j][i_local_max + 1] = (gasconstant * rho[i][j][i_local_max + 1]) * temperature_1D[i_local_max + 1];
                    }

                    //to debug
                    //printf("energy_left: %g energy_right: %g dt: %g\n", energy_left, energy_right, dt);

                    if (dom_counter > 2) {

                        for (n = i_local_min + 1; n < i_local_max; n++) {
                            //here the heat conductivity constants are determined
                            temperature_1D_future[n] = alpha_solid_z * temperature_1D[n + 1] +
                                    (1 - 2 * alpha_solid_z) * temperature_1D[n] +
                                    alpha_solid_z * temperature_1D[n - 1];
                        }

                        for (k = i_local_min + 1; k < i_local_max; k++) {
                            pre[i][j][k] = (gasconstant * rho[i][j][k]) * temperature_1D_future[k];
                        }

                        for (k = i_local_min; k < i_local_max; k++) {
                            temperature_1D[k] = 0.0;
                            temperature_1D_future[k] = 0.0;
                            mass[k] = 0.0;
                            dom_1D[k] = 0;
                        }
                    }

                    i_local_min = i_local_max + 1;

                }
            }
        }

        free_diffusion(&temperature_1D, &temperature_1D_future, &mass, &dom_1D);

#ifdef _OPENMP 
    }
#endif 


    return 0;
}