/**
 * write_raw.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 In this routine the output-files are written. 
 */
int write_tyc(int x, int y, int z, int counter) {
    FILE *fd;
    char filename[600];
    int i, j, k, tmp1;
    double tmp;

    sprintf(filename, "%srho_%i.tyc", output_dir, counter);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "TYCHO Density File\n%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
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


    sprintf(filename, "%stemp_%i.tyc", output_dir, counter);
    fd = fopen(filename, "wb");
    fprintf(fd, "TYCHO Temperature File\n%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                tmp = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                if (isnan(tmp)) {
                    printf("NaN in temperature array %i %i %i\n", i, j, k);
                    printf("%g %g %i\n", pre[i][j][k], rho[i][j][k], dom[i][j][k]);
                    exit(1);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    sprintf(filename, "%spressure_%i.tyc", output_dir, counter);
    fd = fopen(filename, "wb");
    fprintf(fd, "TYCHO Pressure File\n%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                if (dom[i][j][k] == 0) tmp = pre[i][j][k];
                if (dom[i][j][k] == 1) tmp = 0.0;
                if (isnan(tmp)) {
                    printf("NaN in pressure array %i %i %i\n", i, j, k);
                    exit(1);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);


    sprintf(filename, "%svel_%i.tyc", output_dir, counter);
    fd = fopen(filename, "wb");
    fprintf(fd, "TYCHO Velocity File\n%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                tmp = vx[i][j][k];
                if (isnan(tmp)) {
                    printf("NaN in vx array %i %i %i\n", i, j, k);
                    exit(1);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
                tmp = vy[i][j][k];
                if (isnan(tmp)) {
                    printf("NaN in vy array %i %i %i\n", i, j, k);
                    exit(1);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
                tmp = vz[i][j][k];
                if (isnan(tmp)) {
                    printf("NaN in vz array %i %i %i\n", i, j, k);
                    exit(1);
                }
                fwrite(&tmp, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);


    if (with_obstacles == 1) {
        sprintf(filename, "%spressure_on_solid_%i.tyc", output_dir, counter);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO pressure on solid File\n%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    tmp = pressure_on_solid[i][j][k];
                    fwrite(&tmp, 1, sizeof (double), fd);
                }
            }
        }
        fclose(fd);

        if ((with_sound_emitter == 1) && (dB_Map_ready == 1)) {
            sprintf(filename, "%sdB_map_%i.tyc", output_dir, counter);
            fd = fopen(filename, "wb");
            fprintf(fd, "TYCHO dba_map\n%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
            fseek(fd, 200, SEEK_SET);
            for (i = 0; i < x; i++) {
                for (j = 0; j < y; j++) {
                    for (k = 0; k < z; k++) {
                        tmp = dB_map[i][j][k];
                        if (isnan(tmp)) {
                            printf("NaN in dB_map array %i %i %i\n", i, j, k);
                            exit(1);
                        }

                        fwrite(&tmp, 1, sizeof (double), fd);
                    }
                }
            }
            fclose(fd);
        }
    }

    if (advection == 1) {
        sprintf(filename, "%smarker_%i.tyc", output_dir, counter);
        fd = fopen(filename, "wb");
        fprintf(fd, "TYCHO Marker File%g\n%i\n%i\n%i\n%i\n", time_sim, x, y, z, counter);
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

    return 0;
}
