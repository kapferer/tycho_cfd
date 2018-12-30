/*
 * write_ic_ifrit.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Writes the initial model in a IFRIT
 uniform scale binary data file
 The slow order of the nested loops is mandatory here
 */
int write_ic_ifrit(int x, int y, int z) {
    FILE *fd;
    char filename[600];
    int i, j, k;
    int ntemp;
    double tmp, tmp_vx, tmp_vy, tmp_vz;
    float data1;

    sprintf(filename, "%srho_ic.bin", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    ntemp = 12;
    fwrite(&ntemp, 1, sizeof (int), fd);
    fwrite(&x, 1, sizeof (int), fd);
    fwrite(&y, 1, sizeof (int), fd);
    fwrite(&z, 1, sizeof (int), fd);
    fwrite(&ntemp, 1, sizeof (int), fd);
    ntemp = 4 * x * y * z;
    fwrite(&ntemp, 1, sizeof (int), fd);
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp = rho[i][j][k];
                data1 = (float) tmp;
                fwrite(&data1, 1, sizeof (data1), fd);
            }
        }
    }
    fwrite(&ntemp, 4, 1, fd);
    fclose(fd);

    sprintf(filename, "%stemp_ic.bin", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    ntemp = 12;
    fwrite(&ntemp, 1, sizeof (int), fd);
    fwrite(&x, 1, sizeof (int), fd);
    fwrite(&y, 1, sizeof (int), fd);
    fwrite(&z, 1, sizeof (int), fd);
    fwrite(&ntemp, 1, sizeof (int), fd);
    ntemp = 4 * x * y * z;
    fwrite(&ntemp, 1, sizeof (int), fd);
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                data1 = (float) tmp;
                fwrite(&data1, 1, sizeof (data1), fd);
            }
        }
    }
    fwrite(&ntemp, 1, sizeof (int), fd);
    fclose(fd);

    sprintf(filename, "%spressure_ic.bin", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    ntemp = 12;
    fwrite(&ntemp, 1, sizeof (int), fd);
    fwrite(&x, 1, sizeof (int), fd);
    fwrite(&y, 1, sizeof (int), fd);
    fwrite(&z, 1, sizeof (int), fd);
    fwrite(&ntemp, 1, sizeof (int), fd);
    ntemp = 4 * x * y * z;
    fwrite(&ntemp, 1, sizeof (int), fd);
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                if (dom[i][j][k] == 0) tmp = pre[i][j][k];
                if (dom[i][j][k] == 1) tmp = 0.0;
                data1 = (float) tmp;
                fwrite(&data1, 1, sizeof (data1), fd);
            }
        }
    }
    fwrite(&ntemp, 1, sizeof (int), fd);
    fclose(fd);

    sprintf(filename, "%svel_ic.bin", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    ntemp = 12;
    fwrite(&ntemp, 1, sizeof (int), fd);
    fwrite(&x, 1, sizeof (int), fd);
    fwrite(&y, 1, sizeof (int), fd);
    fwrite(&z, 1, sizeof (int), fd);
    fwrite(&ntemp, 1, sizeof (int), fd);
    ntemp = 4 * x * y*z;
    fwrite(&ntemp, 1, sizeof (int), fd);
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp_vx = vy[i][j][k];
                data1 = (float) tmp_vx;
                fwrite(&data1, 1, sizeof (data1), fd);
            }
        }
    }
    fwrite(&ntemp, 1, sizeof (int), fd);
    fwrite(&ntemp, 1, sizeof (int), fd);
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp_vy = vx[i][j][k];
                data1 = (float) tmp_vy;
                fwrite(&data1, 1, sizeof (data1), fd);
            }
        }
    }
    fwrite(&ntemp, 1, sizeof (int), fd);
    fwrite(&ntemp, 1, sizeof (int), fd);
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp_vz = vz[i][j][k];
                data1 = (float) tmp_vz;
                fwrite(&data1, 1, sizeof (data1), fd);
            }
        }
    }
    fwrite(&ntemp, 1, sizeof (int), fd);
    fclose(fd);

    if (with_obstacles == 1) {
        sprintf(filename, "%spressure_on_solid_ic.bin", output_dir);
        fd = fopen(filename, "wb");
        if (fd == NULL) {
            printf("-----------------------------------\n");
            printf("The output directory does not exist\n");
            printf("-----------------------------------\n");
            exit(13);
        }
        ntemp = 12;
        fwrite(&ntemp, 1, sizeof (int), fd);
        fwrite(&x, 1, sizeof (int), fd);
        fwrite(&y, 1, sizeof (int), fd);
        fwrite(&z, 1, sizeof (int), fd);
        fwrite(&ntemp, 1, sizeof (int), fd);
        ntemp = 4 * x * y * z;
        fwrite(&ntemp, 1, sizeof (int), fd);
        for (k = 0; k < z; k++) {
            for (j = 0; j < y; j++) {
                for (i = 0; i < x; i++) {
                    tmp = pressure_on_solid[i][j][k];
                    data1 = (float) tmp;
                    fwrite(&data1, 1, sizeof (data1), fd);
                }
            }
        }
        fwrite(&ntemp, 1, sizeof (int), fd);

        fclose(fd);

        if ((with_sound_emitter == 1) && (dB_Map_ready == 1)) {
            sprintf(filename, "%sdB_map_ic_ic.bin", output_dir);
            fd = fopen(filename, "wb");
            if (fd == NULL) {
                printf("-----------------------------------\n");
                printf("The output directory does not exist\n");
                printf("-----------------------------------\n");
                exit(13);
            }
            ntemp = 12;
            fwrite(&ntemp, 1, sizeof (int), fd);
            fwrite(&x, 1, sizeof (int), fd);
            fwrite(&y, 1, sizeof (int), fd);
            fwrite(&z, 1, sizeof (int), fd);
            fwrite(&ntemp, 1, sizeof (int), fd);
            ntemp = 4 * x * y * z;
            fwrite(&ntemp, 1, sizeof (int), fd);
            for (k = 0; k < z; k++) {
                for (j = 0; j < y; j++) {
                    for (i = 0; i < x; i++) {
                        tmp = dB_map[i][j][k];
                        data1 = (float) tmp;
                        fwrite(&data1, 1, sizeof (data1), fd);
                    }
                }
            }
            fwrite(&ntemp, 1, sizeof (int), fd);

            fclose(fd);
        }
    }

    if (advection == 1) {

        sprintf(filename, "%smarker_ic.bin", output_dir);
        fd = fopen(filename, "wb");
        if (fd == NULL) {
            printf("-----------------------------------\n");
            printf("The output directory does not exist\n");
            printf("-----------------------------------\n");
            exit(13);
        }
        ntemp = 12;
        fwrite(&ntemp, 1, sizeof (int), fd);
        fwrite(&x, 1, sizeof (int), fd);
        fwrite(&y, 1, sizeof (int), fd);
        fwrite(&z, 1, sizeof (int), fd);
        fwrite(&ntemp, 1, sizeof (int), fd);
        ntemp = 4 * x * y * z;
        fwrite(&ntemp, 1, sizeof (int), fd);
        for (k = 0; k < z; k++) {
            for (j = 0; j < y; j++) {
                for (i = 0; i < x; i++) {
                    tmp = marker[i][j][k];
                    data1 = (float) tmp;
                    fwrite(&data1, 1, sizeof (data1), fd);
                }
            }
        }
        fwrite(&ntemp, 1, sizeof (int), fd);
        fclose(fd);

    }

    return 0;
}
