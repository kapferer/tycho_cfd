/**
 * wind.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 A function to increase the density/pressure by a fixed factor
 */
int density_increaser(int x, int y, int z) {

    int i, j, k;
    double max;

    max = 0.0;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                rho[i][j][k] = 1.1 * rho[i][j][k];
                pre[i][j][k] = 1.1 * pre[i][j][k];
            }
        }
    }

    return 0;
}

/*!
 The Wind routine, which gives a constant flow with wind-speed velocity
 */
int wind(int x, int y, int z) {
    int i, j, k;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {

                if (wind_marker[i][j][k] == 2) {
                    vx[i][j][k] = wind_speed;
                    vy[i][j][k] = vy[i][j][k];
                    vz[i][j][k] = vz[i][j][k];
                }
                if (wind_marker[i][j][k] == 3) {
                    vx[i][j][k] = -1 * wind_speed;
                    vy[i][j][k] = vy[i][j][k];
                    vz[i][j][k] = vz[i][j][k];
                }
                if (wind_marker[i][j][k] == 4) {
                    vx[i][j][k] = vx[i][j][k];
                    vy[i][j][k] = wind_speed;
                    vz[i][j][k] = vz[i][j][k];
                }
                if (wind_marker[i][j][k] == 5) {
                    vx[i][j][k] = vx[i][j][k];
                    vy[i][j][k] = -1 * wind_speed;
                    vz[i][j][k] = vz[i][j][k];
                }
                if (wind_marker[i][j][k] == 6) {
                    vx[i][j][k] = vx[i][j][k];
                    vy[i][j][k] = vy[i][j][k];
                    vz[i][j][k] = wind_speed;
                }
                if (wind_marker[i][j][k] == 7) {
                    vx[i][j][k] = vx[i][j][k];
                    vy[i][j][k] = vy[i][j][k];
                    vz[i][j][k] = -1 * wind_speed;
                }
            }
        }
    }


    return 0;
}
