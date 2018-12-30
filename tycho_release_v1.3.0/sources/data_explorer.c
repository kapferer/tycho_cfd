/*
 * data_explorer.c
 * 
 * Author: Wolfgang Kapferer
 * 
 * Pressure on solid 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "variables_global.h"
#include "prototypes.h"

/*!
 Quantities of the fluid are integrated here. Here you can start
 to implement you owns investigations of integrated quantities of
 the fluid.
 */
int velocity_field_analyser(int x, int y, int z) {

    FILE *fd;
    char filename[600];
    int i, j, k;
    double velocity_x, velocity_y, velocity_z;
    double velo_ratio;

    sprintf(filename, "%stotal_quantities_of_fluid.txt", output_dir);
    fd = fopen(filename, "aw");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }

    velocity_x = velocity_y = velocity_z = velo_ratio = 0.0;

    if (dimension == 3) {
#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, velocity_x, velocity_y, velocity_z, \
                       rho, dom, vx, vy, vz)
        {
#endif

            //Here we calculate a density * velocity quantity
            for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif

                for (j = 0; j < y; j++) {
                    for (k = 0; k < z; k++) {
                        if (dom[i][j][k] == 0) {
                            velocity_x += rho[i][j][k] * vx[i][j][k];
                            velocity_y += rho[i][j][k] * vy[i][j][k];
                            velocity_z += rho[i][j][k] * vz[i][j][k];
                        }
                    }
                }
            }

#ifdef _OPENMP
        }
#endif
    }

    if (dimension == 2) {
#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, velocity_x, velocity_y, velocity_z, \
                       rho, dom, vx, vy, vz)
        {
#endif

            k = 0;
            //Here we calculate a density * velocity quantity
            for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif

                for (j = 0; j < y; j++) {
                    if (dom[i][j][k] == 0) {
                        velocity_x += rho[i][j][k] * vx[i][j][k];
                        velocity_y += rho[i][j][k] * vy[i][j][k];
                        velocity_z += rho[i][j][k] * vz[i][j][k];
                    }
                }
            }

#ifdef _OPENMP
        }
#endif
    }

    if (dimension == 2) {
        velo_ratio = velocity_x / velocity_y;
        fprintf(fd, "%e %e %e %e\n", time_sim, velo_ratio, velocity_x, velocity_y);
    }

    if (dimension == 3) {
        fprintf(fd, "%e %e %e %e\n", time_sim, velocity_x, velocity_y, velocity_z);
    }

    fclose(fd);

    return 0;
}

/*!
 Here the total force due to the pressure on the solid is determined.
 */
int pressure_on_solid_calc(int x, int y, int z) {

    FILE *fd;
    char filename[600];
    int i, j, k, l;
    double total_pressure_force;

    total_pressure_force = 0.0;

    sprintf(filename, "%stotal_pressure_on_solid.txt", output_dir);
    fd = fopen(filename, "aw");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }

    if (dimension == 3) {
#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, l) \
                shared(x, y, z, total_pressure_force, \
                pressure_on_solid, dimension, dom, xmax, \
                ymax)
        {
#endif

            for (i = 1; i < x - 1; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif    

                for (j = 1; j < y - 1; j++) {
                    for (k = 1; k < z - 1; k++) {
                        if (dom[i][j][k] == 1) {
                            for (l = -1; l < 2; l++) {
                                if (l == 0) continue;
                                if (dom[i + l][j][k] == 0) {
                                    total_pressure_force += pressure_on_solid[i + l][j][k]*(xmax / x)*(ymax / y);
                                }
                                if (dom[i][j + l][k] == 0) {
                                    total_pressure_force += pressure_on_solid[i][j + l][k]*(xmax / x)*(ymax / y);
                                }
                                if (dom[i][j][k + l] == 0) {
                                    total_pressure_force += pressure_on_solid[i][j][k + l]*(xmax / x)*(ymax / y);
                                }
                            }
                        }
                    }
                }
            }

#ifdef _OPENMP
        }
#endif
    }


    if (dimension == 2) {
#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, l, \
                pressure_on_solid, \
                dimension, dom, \
                total_pressure_force, \
                xmax, ymax)
        {
#endif


            k = 0;
            for (i = 1; i < x - 1; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif    

                for (j = 1; j < y - 1; j++) {
                    if (dom[i][j][k] == 1) {
                        for (l = -1; l < 2; l++) {
                            if (l == 0) continue;
                            if (dom[i + l][j][k] == 0) {
                                total_pressure_force += pressure_on_solid[i + l][j][k]*(xmax / x);
                            }
                            if (dom[i][j + l][k] == 0) {
                                total_pressure_force += pressure_on_solid[i][j + l][k]*(xmax / x);
                            }
                        }
                    }
                }
            }

#ifdef _OPENMP
        }
#endif
    }

    fprintf(fd, "%e %e\n", time_sim, total_pressure_force);
    fclose(fd);

    return 0;
}


