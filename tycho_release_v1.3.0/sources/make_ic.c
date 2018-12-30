/*
 * make_ic.c
 *
 *      Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 An initial model is build in this function. If the parameter file defines strat_const_atmos
 a stratified atmosphere is build. If the parameter is set to zero
 a constant pressure/density distribution is realized.
 */
int make_ic(int x, int y, int z) {
    int i, j, k;
    double velo_x, velo_y, velo_z;

    double T0, T, A, rho0, pre0;

    //T(z)=T_0 + A*(z-z0)
    T0 = 288.15; //[K]
    A = -6.5E-3; //[K/m]

    rho0 = 1.229; // Density at sea level [kg/m^3]
    pre0 = 1.013E5; //pressure at sea level [N/m^2]

    velo_x = small;
    velo_y = small;
    velo_z = small;

    if (starting_flow == 1) {
        if (bound.down == 4) velo_y = inflow_velocity;
        if (bound.up == 4) velo_y = inflow_velocity;
        if (bound.left == 4) velo_x = inflow_velocity;
        if (bound.right == 4) velo_x = inflow_velocity;
        if (bound.front == 4) velo_z = inflow_velocity;
        if (bound.back == 4) velo_z = inflow_velocity;
    }

    if (strat_const_atmos == 1) {
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {

                    T = T0 + A * ((ymax / y) * (y - j));

                    rho[i][j][k] = rho0 * pow((T / T0), (-grav_acc / (gasconstant * A) + 1));
                    pre[i][j][k] = pre0 * pow((T / T0), (-grav_acc / (gasconstant * A)));

                    vx[i][j][k] = velo_x;
                    vy[i][j][k] = velo_y;
                    vz[i][j][k] = velo_z;

                }
            }
        }
    }

    if (strat_const_atmos == 0) {
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {

                    rho[i][j][k] = 1.229;
                    pre[i][j][k] = 300 * gasconstant * rho[i][j][k];

                    vx[i][j][k] = velo_x;
                    vy[i][j][k] = velo_y;
                    vz[i][j][k] = velo_z;

                }
            }
        }
    }

    // If you need a conversion from temperature to pressure
    // calculate_pressure(x, y, z, temperature);

    //just to add explosions
    //insert_pressure(x, y, z);

    return 0;
}

/*!
 If one needs the pressure from temperature.
 */
int calculate_pressure(int x, int y, int z, double temperature) {
    int i, j, k;


    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                pre[i][j][k] = temperature * gasconstant * rho[i][j][k];
            }
        }
    }

    return 0;
}

/*!
 The form of the obstacle is defined in this function
 example of just a round thing.
 */
int initiate_domain(int x, int y, int z) {
    int i, j, k;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                int i2 = i - x / 2.0;
                int j2 = j - y / 2.0;
                int xy = sqrt(i2 * i2 + j2 * j2);

                if (xy < x / 20) {
                    dom[i][j][k] = 1;
                    rho[i][j][k] = obstacle_density;
                    pre[i][j][k] = obstacle_temperature;

                    vx[i][j][k] = 0.0;
                    vy[i][j][k] = 0.0;
                    vz[i][j][k] = 0.0;
                } else {
                    dom[i][j][k] = 0;
                }

            }
        }
    }

    printf("Initiate domain done\n");

    return 0;
}

/*!
 If the dom file is given the density and temperature for the solid are given here.
 Note that the pre array holds the temperature for the solid
 */
int initiate_domain_dom(int x, int y, int z) {
    int i, j, k;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {

                if (dom[i][j][k] == 1) {
                    rho[i][j][k] = obstacle_density;
                    pre[i][j][k] = obstacle_temperature;

                    vx[i][j][k] = 0.0;
                    vy[i][j][k] = 0.0;
                    vz[i][j][k] = 0.0;
                }
            }
        }
    }

    printf("Initiate domain done\n");

    return 0;
}

/*!
 If the marker file is given the density for the marker is defined here.
 */
int initiate_domain_marker(int x, int y, int z) {
    int i, j, k;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {

                if (marker[i][j][k] == 1) {
                    marker[i][j][k] = marker_density;
                }
            }
        }
    }

    printf("Initiate domain done\n");

    return 0;
}

/*!
 A Sod Shock Tube is generated here
 */
int make_sod_ic(int x, int y, int z) {
    int i, j, k;
    float radius, radius1;

    if (dimension == 2) radius = sqrt(pow(x, 2) + pow(y, 2));
    if (dimension == 3) radius = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {

                if (dimension == 1) {
                    if (i < x / 2) {
                        rho[i][j][k] = 1.0;
                        pre[i][j][k] = 1.0 * 5E4;
                        vx[i][j][k] = 0.0;
                        vy[i][j][k] = 0.0;
                        vz[i][j][k] = 0.0;
                    } else {
                        rho[i][j][k] = 0.15;
                        pre[i][j][k] = 0.10 * 5E4;
                        vx[i][j][k] = 0.0;
                        vy[i][j][k] = 0.0;
                        vz[i][j][k] = 0.0;
                    }
                }

                if (dimension == 2) {

                    radius1 = sqrt(pow(i, 2) + pow(j, 2));

                    if (radius1 < radius / 2) {
                        rho[i][j][k] = 1.0;
                        pre[i][j][k] = 1.0 * 1E4;
                        vx[i][j][k] = 0.0;
                        vy[i][j][k] = 0.0;
                        vz[i][j][k] = 0.0;
                    } else {
                        rho[i][j][k] = 0.15;
                        pre[i][j][k] = 0.10 * 1E2;
                        vx[i][j][k] = 0.0;
                        vy[i][j][k] = 0.0;
                        vz[i][j][k] = 0.0;
                    }
                }

                if (dimension == 3) {

                    radius1 = sqrt(pow(i, 2) + pow(j, 2) + pow(k, 2));

                    if (radius1 <= radius / 2) {
                        rho[i][j][k] = 1.0;
                        pre[i][j][k] = 1.0;
                        vx[i][j][k] = 0.0;
                        vy[i][j][k] = 0.0;
                        vz[i][j][k] = 0.0;
                    } else {
                        rho[i][j][k] = 0.15;
                        pre[i][j][k] = 0.10;
                        vx[i][j][k] = 0.0;
                        vy[i][j][k] = 0.0;
                        vz[i][j][k] = 0.0;
                    }
                }

            }
        }
    }

    return 0;
}

/*!
 A simple Kelvin Helmholtz Instability Setup
 */
int make_kh_instabilities(int x, int y, int z) {

    int i, j, k;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {

                if ((j < y / 4) || (j > 3 * y / 4)) {
                    rho[i][j][k] = 1;
                    pre[i][j][k] = 1E5;
                    vx[i][j][k] = -50;
                    vy[i][j][k] = 0.0;
                    vz[i][j][k] = 0.0;
                } else {
                    rho[i][j][k] = 2;
                    pre[i][j][k] = 1E5;
                    vx[i][j][k] = 0.0;
                    vy[i][j][k] = 0.0;
                    vz[i][j][k] = 0.0;
                }

                if (dimension == 3) {
                    if (j == y / 4) vy[i][j][k] = sin(xmax / x * i) * cos(zmax / z * k);
                    if (j == 3 * y / 4) vy[i][j][k] = -1 * sin(xmax / x * i) * cos(zmax / z * k);
                }

                if (dimension == 2) {
                    if (j == y / 4) vy[i][j][k] = sin(xmax / x * i);
                    if (j == 3 * y / 4) vy[i][j][k] = -1 * sin(xmax / x * i);
                }
            }
        }
    }



    return 0;
}

/*!
 To insert pressure seeds in some regions
 */
int insert_pressure(int x, int y, int z) {

    int i, j, k;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                if (dom[i][j][k] != 1) {

                    if (j < y / 4) {
                        if (round((double) rand() / (double) RAND_MAX) == 1) {
                            pre[i][j][k] = 2000 * gasconstant * rho[i][j][k];
                        }
                    }
                }
            }
        }
    }


    return 0;
}