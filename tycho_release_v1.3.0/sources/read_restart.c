/*
 * read_restart.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
Read restart from an native TYCHO file.
 */
int read_restart(int x, int y, int z) {
    int i, j, k, tmp1;
    double tmp;
    char filename[200];
    char filename1[200];
    char buffer[10];
    FILE *fd;
    FILE *fd1;
    int x1, y1, z1;

    //The density distribution.
    sprintf(filename, "%srho_restart.tyc", output_dir);
    fd = fopen(filename, "rb");
    if (fd == NULL) {
        printf("-----------------------------\n");
        printf("No restart density file found\n");
        printf("-----------------------------\n");
        exit(12);
    }

    fgets(buffer, 100, fd);
    fgets(buffer, 100, fd);
    time_sim = atof(buffer);
    fgets(buffer, 100, fd);
    x1 = atoi(buffer);
    fgets(buffer, 100, fd);
    y1 = atoi(buffer);
    fgets(buffer, 100, fd);
    z1 = atoi(buffer);
    fgets(buffer, 100, fd);
    counter = atoi(buffer);
    fgets(buffer, 100, fd);
    counter_restart = atoi(buffer);
    fseek(fd, 200, SEEK_SET);
    fseek(fd1, 200, SEEK_SET);

    if ((x != x1) || (y != y1) || (z != z1)) {
        printf("Restart File dimension mismatch.");
        exit(21);
    }

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp, 1, sizeof (double), fd);
                rho[i][j][k] = tmp;
            }
        }
    }

    fclose(fd);

    //The temperature distribution
    sprintf(filename, "%stemp_restart.tyc", output_dir);
    fd = fopen(filename, "rb");
    if (fd == NULL) {
        printf("---------------------------------\n");
        printf("No restart pressure file found\n");
        printf("---------------------------------\n");
        exit(13);
    }
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp, 1, sizeof (double), fd);
                pre[i][j][k] = tmp * (gasconstant * rho[i][j][k]);
            }
        }
    }
    fclose(fd);


    sprintf(filename, "%svel_restart.tyc", output_dir);
    fd = fopen(filename, "rb");
    if (fd == NULL) {
        printf("------------------------------\n");
        printf("No restart velocity file found\n");
        printf("------------------------------\n");
        exit(14);
    }
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp, 1, sizeof (double), fd);
                vx[i][j][k] = tmp;
                fread(&tmp, 1, sizeof (double), fd);
                vy[i][j][k] = tmp;
                fread(&tmp, 1, sizeof (double), fd);
                vz[i][j][k] = tmp;
            }
        }
    }
    fclose(fd);


    if (with_obstacles == 1) {
        //The obstacles distribution
        sprintf(filename, "%sdom_restart.tyc", output_dir);
        fd = fopen(filename, "rb");
        if (fd == NULL) {
            printf("---------------------------------\n");
            printf("No restart obstacles file found\n");
            printf("---------------------------------\n");
            exit(13);
        }
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    fread(&tmp1, 1, sizeof (int), fd);
                    dom[i][j][k] = tmp1;
                }
            }
        }
        fclose(fd);
    }

    if (with_sound_emitter == 1) {
        //The obstacles distribution
        sprintf(filename, "%ssound_emitter_restart.tyc", output_dir);
        fd = fopen(filename, "rb");
        if (fd == NULL) {
            printf("---------------------------------\n");
            printf("No restart sound-emitter file found\n");
            printf("---------------------------------\n");
            exit(13);
        }
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    fread(&tmp1, 1, sizeof (int), fd);
                    soundemitter[i][j][k] = tmp1;
                }
            }
        }
        fclose(fd);

        //The obstacles distribution
        sprintf(filename, "%spressure_integrated_restart.tyc", output_dir);
        fd = fopen(filename, "rb");
        if (fd == NULL) {
            printf("---------------------------------\n");
            printf("No restart pressure-integrated restart file found\n");
            printf("---------------------------------\n");
            exit(13);
        }
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    fread(&tmp, 1, sizeof (double), fd);
                    pressure_integrated[i][j][k] = tmp;
                }
            }
        }
        fclose(fd);
    }


    if (advection == 1) {
        //The marker distribution
        sprintf(filename, "%smarker_restart.tyc", output_dir);
        fd = fopen(filename, "rb");
        if (fd == NULL) {
            printf("---------------------------------\n");
            printf("No restart marker file found\n");
            printf("---------------------------------\n");
            exit(13);
        }
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    fread(&tmp, 1, sizeof (double), fd);
                    marker[i][j][k] = tmp;
                }
            }
        }
        fclose(fd);
    }

    if (wind_on_off == 1) {
        //The wind marker distribution
        sprintf(filename, "%swind_restart.tyc", output_dir);
        fd = fopen(filename, "rb");
        if (fd == NULL) {
            printf("---------------------------------\n");
            printf("No restart wind marker file found\n");
            printf("---------------------------------\n");
            exit(13);
        }
        fseek(fd, 200, SEEK_SET);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {
                    fread(&tmp, 1, sizeof (double), fd);
                    wind_marker[i][j][k] = tmp;
                }
            }
        }
        fclose(fd);
    }

    return 0;
}