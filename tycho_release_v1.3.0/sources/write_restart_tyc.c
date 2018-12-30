/*
 * write_restart_tyc.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Writes the restart files
 */
int write_restart_tyc(int x, int y, int z) {
    FILE *fd;
    char filename[600];
    int i, j, k, tmp1;
    double tmp;

    sprintf(filename, "%srho_restart.tyc", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "TYCHO Density File\n%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                tmp = rho[i][j][k];
                if (isnan(tmp)) {
                    printf("NaN in density array %i %i %i\n", i, j, k);
                    exit(0);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
            }
        }
    }

    fclose(fd);

    sprintf(filename, "%stemp_restart.tyc", output_dir);
    fd = fopen(filename, "wb");
    fprintf(fd, "TYCHO Temperature File\n%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                tmp = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                if (isnan(tmp)) {
                    printf("NaN in temperature array %i %i %i\n", i, j, k);
                    exit(0);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    sprintf(filename, "%svel_restart.tyc", output_dir);
    fd = fopen(filename, "wb");
    fprintf(fd, "TYCHO Velocity File%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                tmp = vx[i][j][k];
                fwrite(&tmp, 1, sizeof (double), fd);
                tmp = vy[i][j][k];
                fwrite(&tmp, 1, sizeof (double), fd);
                tmp = vz[i][j][k];
                fwrite(&tmp, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    if (with_obstacles == 1) {
        sprintf(filename, "%sdom_restart.tyc", output_dir);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO Dom File%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    tmp1 = dom[i][j][k];
                    if (isnan(tmp1)) {
                        printf("NaN in dom array %i %i %i\n", i, j, k);
                        exit(1);
                    }
                    fwrite(&tmp1, 1, sizeof (int), fd);
                }
            }
        }
        fclose(fd);
    }

    if (with_sound_emitter == 1) {
        sprintf(filename, "%ssound_emitter_restart.tyc", output_dir);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO Dom File%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    tmp1 = soundemitter[i][j][k];
                    if (isnan(tmp1)) {
                        printf("NaN in sound_emitter array %i %i %i\n", i, j, k);
                        exit(1);
                    }
                    fwrite(&tmp1, 1, sizeof (int), fd);
                }
            }
        }
        fclose(fd);
        
        sprintf(filename, "%spressure_integrated_restart.tyc", output_dir);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO Pressure Integrated File%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    tmp = pressure_integrated[i][j][k];
                    if (isnan(tmp)) {
                        printf("NaN in spressure_integrated_restart array %i %i %i\n", i, j, k);
                        exit(1);
                    }
                    fwrite(&tmp, 1, sizeof (double), fd);
                }
            }
        }
        fclose(fd);
    }

    if (advection == 1) {

        sprintf(filename, "%smarker_restart.tyc", output_dir);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO Marker File%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    tmp = marker[i][j][k];
                    if (isnan(tmp)) {
                        printf("NaN in marker array %i %i %i\n", i, j, k);
                        exit(1);
                    }
                    fwrite(&tmp, 1, sizeof (double), fd);
                }
            }
        }
        fclose(fd);

    }

    if (wind_on_off == 1) {

        sprintf(filename, "%swind_restart.tyc", output_dir);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO Wind File%g\n%i\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter, counter_restart);
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    tmp = wind_marker[i][j][k];
                    if (isnan(tmp)) {
                        printf("NaN in wind marker array %i %i %i\n", i, j, k);
                        exit(1);
                    }
                    fwrite(&tmp, 1, sizeof (double), fd);
                }
            }
        }
        fclose(fd);

    }

    return 0;

}