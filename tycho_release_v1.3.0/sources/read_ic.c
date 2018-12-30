/*
 * read_ic.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 If you want to read in initial conditions from an native TYCHO file.
 */
int read_ic(int x, int y, int z) {
    int i, j, k, tmp1;
    double tmp;
    char filename[200];
    char buffer[10];
    FILE *fd;

    //The initial density distribution.
    sprintf(filename, "%s", dens_ic);
    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("-----------------------------\n");
        printf("No initial density file found\n");
        printf("-----------------------------\n");
        exit(12);
    }
    fgets(buffer, 100, fd);
    fgets(buffer, 100, fd);
    time_sim = atof(buffer);
    fgets(buffer, 100, fd);
    x = atoi(buffer);
    fgets(buffer, 100, fd);
    y = atoi(buffer);
    fgets(buffer, 100, fd);
    z = atoi(buffer);
    fgets(buffer, 100, fd);
    counter = atoi(buffer);
    fgets(buffer, 100, fd);
    counter_restart = atoi(buffer);
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp, 1, sizeof (double), fd);
                rho[i][j][k] = tmp;
            }
        }
    }
    fclose(fd);

    //The initial temperature distribution
    sprintf(filename, "%s", temp_ic);
    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("---------------------------------\n");
        printf("No initial temperature file found\n");
        printf("---------------------------------\n");
        exit(13);
    }
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp, 1, sizeof (double), fd);
                pre[i][j][k] = (tmp) * gasconstant * rho[i][j][k];
                ;
            }
        }
    }
    fclose(fd);

    //The initial velocity distribution
    if (intial_velocity_file == 1) {
        sprintf(filename, "%s", vel_ic);
        fd = fopen(filename, "r");
        if (fd == NULL) {
            printf("------------------------------\n");
            printf("No initial velocity file found\n");
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
    }

    if (with_obstacles == 1) {
        //The initial obstacles distribution
        sprintf(filename, "%s", obst_ic);
        fd = fopen(filename, "r");
        if (fd == NULL) {
            printf("---------------------------------\n");
            printf("No initial obstacles file found\n");
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

    return 0;
}

int read_dom(int x, int y, int z) {
    int i, j, k, tmp1;
    char filename[200];
    FILE *fd;

    //The initial obstacles distribution
    sprintf(filename, "%s", obst_ic);
    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("---------------------------------\n");
        printf("No initial obstacles file found\n");
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


    return 0;
}

int read_soundemitter(int x, int y, int z) {
    int i, j, k, tmp1;
    char filename[200];
    FILE *fd;

    //The initial obstacles distribution
    sprintf(filename, "%s", sound_ic);
    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("---------------------------------\n");
        printf("No initial sound-emitter file found\n");
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


    return 0;
}
int read_marker_file(int x, int y, int z) {
    int i, j, k, tmp1;
    char filename[200];
    FILE *fd;

    //The initial marker distribution
    sprintf(filename, "%s", marker_ic);
    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("---------------------------------\n");
        printf("No initial marker file found\n");
        printf("---------------------------------\n");
        exit(13);
    }
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp1, 1, sizeof (int), fd);
                marker[i][j][k] = tmp1;
            }
        }
    }
    fclose(fd);


    return 0;
}

int read_wind_file(int x, int y, int z) {
    int i, j, k, tmp1;
    char filename[200];
    FILE *fd;

    //The initial marker distribution
    sprintf(filename, "%s", wind_ic);
    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("---------------------------------\n");
        printf("No wind file found\n");
        printf("---------------------------------\n");
        exit(13);
    }
    fseek(fd, 200, SEEK_SET);
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            for (k = 0; k < z; k++) {
                fread(&tmp1, 1, sizeof (int), fd);
                wind_marker[i][j][k] = tmp1;
            }
        }
    }
    fclose(fd);


    return 0;
}