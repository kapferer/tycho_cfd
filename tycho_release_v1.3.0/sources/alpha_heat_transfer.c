/*
 * alpha_heat_transfer.c
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
 Calculates the heat-transfer coefficient for the 
 heat transfer between obstacles and gas. 
 Here you can put your special physics for
 the heat transfer coefficient in.
 Example for flat wall.
 */
double alpha_heat_transfer(int x, int y, int z, int direction) {

    double heat_transfer_coefficient;
    double velocity_component_at_obstacle;

    //the x-y-z sweeps
    if (direction == 0) velocity_component_at_obstacle = sqrt(pow(vy[x][y][z],2)+pow(vz[x][y][z],2));
    if (direction == 1) velocity_component_at_obstacle = sqrt(pow(vx[x][y][z],2)+pow(vz[x][y][z],2));
    if (direction == 2) velocity_component_at_obstacle = sqrt(pow(vx[x][y][z],2)+pow(vy[x][y][z],2));


    if (velocity_component_at_obstacle <= 5) heat_transfer_coefficient = 5.6 + 4.0 * velocity_component_at_obstacle;
    if (velocity_component_at_obstacle > 5) heat_transfer_coefficient = 7.2 + pow(velocity_component_at_obstacle, 0.78);

    return heat_transfer_coefficient;
}
