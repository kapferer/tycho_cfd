/*
 * write_ic_vtk.c
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
 Writes the initial model in a RECTILINEAR_GRID VTK-File
 VTK binaries are in BigEndian format
 The slow order of the nested loops is mandatory here
 */
int write_ic_vtk(int x, int y, int z) {
    FILE *fd;
    char filename[600];
    int i, j, k;
    long long int tmp1, tmp1_vx, tmp1_vy, tmp1_vz;
    double tmp, tmp_vx, tmp_vy, tmp_vz;

    sprintf(filename, "%srho_ic.vtk", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# vtk DataFile Version 3.0\n");
    fprintf(fd, "TYCHO Density File - Time:%g\n", time_sim / minute_in_sec);
    fprintf(fd, "BINARY\n");
    fprintf(fd, "DATASET RECTILINEAR_GRID\n");
    fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
    fprintf(fd, "X_COORDINATES %i double\n", x);
    for (i = 0; i < x; i++) {
        tmp = ((double) xmax / (double) x) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nY_COORDINATES %i double\n", y);
    for (i = 0; i < y; i++) {
        tmp = ((double) ymax / (double) y) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nZ_COORDINATES %i double\n", z);
    for (i = 0; i < z; i++) {
        tmp = ((double) zmax / (double) z) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
    fprintf(fd, "SCALARS density double 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp = rho[i][j][k];
                tmp1 = ntohll(*((long long int*) & tmp));
                fwrite(&tmp1, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    sprintf(filename, "%spressure_ic.vtk", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# vtk DataFile Version 3.0\n");
    fprintf(fd, "TYCHO Pressure File - Time:%g\n", time_sim / minute_in_sec);
    fprintf(fd, "BINARY\n");
    fprintf(fd, "DATASET RECTILINEAR_GRID\n");
    fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
    fprintf(fd, "X_COORDINATES %i double\n", x);
    for (i = 0; i < x; i++) {
        tmp = ((double) xmax / (double) x) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nY_COORDINATES %i double\n", y);
    for (i = 0; i < y; i++) {
        tmp = ((double) ymax / (double) y) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nZ_COORDINATES %i double\n", z);
    for (i = 0; i < z; i++) {
        tmp = ((double) zmax / (double) z) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
    fprintf(fd, "SCALARS pressure double 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                if (dom[i][j][k] == 0) tmp = pre[i][j][k];
                if (dom[i][j][k] == 1) tmp = 0.0;
                tmp1 = ntohll(*((long long int*) & tmp));
                fwrite(&tmp1, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    sprintf(filename, "%stemp_ic.vtk", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# vtk DataFile Version 3.0\n");
    fprintf(fd, "TYCHO Temperature File - Time:%g\n", time_sim / minute_in_sec);
    fprintf(fd, "BINARY\n");
    fprintf(fd, "DATASET RECTILINEAR_GRID\n");
    fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
    fprintf(fd, "X_COORDINATES %i double\n", x);
    for (i = 0; i < x; i++) {
        tmp = ((double) xmax / (double) x) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nY_COORDINATES %i double\n", y);
    for (i = 0; i < y; i++) {
        tmp = ((double) ymax / (double) y) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nZ_COORDINATES %i double\n", z);
    for (i = 0; i < z; i++) {
        tmp = ((double) zmax / (double) z) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
    fprintf(fd, "SCALARS temperature double 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                tmp1 = ntohll(*((long long int*) & tmp));
                fwrite(&tmp1, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    sprintf(filename, "%svel_ic.vtk", output_dir);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# vtk DataFile Version 3.0\n");
    fprintf(fd, "TYCHO Velocity File - Time:%g\n", time_sim / minute_in_sec);
    fprintf(fd, "BINARY\n");
    fprintf(fd, "DATASET RECTILINEAR_GRID\n");
    fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
    fprintf(fd, "X_COORDINATES %i double\n", x);
    for (i = 0; i < x; i++) {
        tmp = ((double) xmax / (double) x) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nY_COORDINATES %i double\n", y);
    for (i = 0; i < y; i++) {
        tmp = ((double) ymax / (double) y) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nZ_COORDINATES %i double\n", z);
    for (i = 0; i < z; i++) {
        tmp = ((double) zmax / (double) z) * (double) i;
        tmp1 = ntohll(*((long long int*) & tmp));
        fwrite(&tmp1, 1, sizeof (double), fd);
    }
    fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
    fprintf(fd, "VECTORS velocity double 1\n");
    for (k = 0; k < z; k++) {
        for (j = 0; j < y; j++) {
            for (i = 0; i < x; i++) {
                tmp_vx = vx[i][j][k];
                tmp1_vx = ntohll(*((long long int*) & tmp_vx));
                tmp_vy = vy[i][j][k];
                tmp1_vy = ntohll(*((long long int*) & tmp_vy));
                tmp_vz = vz[i][j][k];
                tmp1_vz = ntohll(*((long long int*) & tmp_vz));

                fwrite(&tmp1_vx, 1, sizeof (double), fd);
                fwrite(&tmp1_vy, 1, sizeof (double), fd);
                fwrite(&tmp1_vz, 1, sizeof (double), fd);
            }
        }
    }
    fclose(fd);

    if (with_obstacles == 1) {
        sprintf(filename, "%spressure_on_solid_ic.vtk", output_dir);
        fd = fopen(filename, "wb");
        if (fd == NULL) {
            printf("-----------------------------------\n");
            printf("The output directory does not exist\n");
            printf("-----------------------------------\n");
            exit(13);
        }
        fprintf(fd, "# vtk DataFile Version 3.0\n");
        fprintf(fd, "TYCHO pressure on solid File - Time:%g\n", time_sim / minute_in_sec);
        fprintf(fd, "BINARY\n");
        fprintf(fd, "DATASET RECTILINEAR_GRID\n");
        fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
        fprintf(fd, "X_COORDINATES %i double\n", x);
        for (i = 0; i < x; i++) {
            tmp = ((double) xmax / (double) x) * (double) i;
            tmp1 = ntohll(*((long long int*) & tmp));
            fwrite(&tmp1, 1, sizeof (double), fd);
        }
        fprintf(fd, "\nY_COORDINATES %i double\n", y);
        for (i = 0; i < y; i++) {
            tmp = ((double) ymax / (double) y) * (double) i;
            tmp1 = ntohll(*((long long int*) & tmp));
            fwrite(&tmp1, 1, sizeof (double), fd);
        }
        fprintf(fd, "\nZ_COORDINATES %i double\n", z);
        for (i = 0; i < z; i++) {
            tmp = ((double) zmax / (double) z) * (double) i;
            tmp1 = ntohll(*((long long int*) & tmp));
            fwrite(&tmp1, 1, sizeof (double), fd);
        }
        fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
        fprintf(fd, "SCALARS pressure_on_solid double 1\n");
        fprintf(fd, "LOOKUP_TABLE default\n");
        for (k = 0; k < z; k++) {
            for (j = 0; j < y; j++) {
                for (i = 0; i < x; i++) {
                    tmp = pressure_on_solid[i][j][k];
                    tmp1 = ntohll(*((long long int*) & tmp));
                    fwrite(&tmp1, 1, sizeof (double), fd);
                }
            }
        }
        fclose(fd);
        if ((with_sound_emitter == 1) && (dB_Map_ready == 1)) {
            sprintf(filename, "%TYCHO dB_map_ic.vtk", output_dir);
            fd = fopen(filename, "wb");
            if (fd == NULL) {
                printf("-----------------------------------\n");
                printf("The output directory does not exist\n");
                printf("-----------------------------------\n");
                exit(13);
            }
            fprintf(fd, "# vtk DataFile Version 3.0\n");
            fprintf(fd, "TYCHO dBA map file - Time:%g\n", time_sim / minute_in_sec);
            fprintf(fd, "BINARY\n");
            fprintf(fd, "DATASET RECTILINEAR_GRID\n");
            fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
            fprintf(fd, "X_COORDINATES %i double\n", x);
            for (i = 0; i < x; i++) {
                tmp = ((double) xmax / (double) x) * (double) i;
                tmp1 = ntohll(*((long long int*) & tmp));
                fwrite(&tmp1, 1, sizeof (double), fd);
            }
            fprintf(fd, "\nY_COORDINATES %i double\n", y);
            for (i = 0; i < y; i++) {
                tmp = ((double) ymax / (double) y) * (double) i;
                tmp1 = ntohll(*((long long int*) & tmp));
                fwrite(&tmp1, 1, sizeof (double), fd);
            }
            fprintf(fd, "\nZ_COORDINATES %i double\n", z);
            for (i = 0; i < z; i++) {
                tmp = ((double) zmax / (double) z) * (double) i;
                tmp1 = ntohll(*((long long int*) & tmp));
                fwrite(&tmp1, 1, sizeof (double), fd);
            }
            fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
            fprintf(fd, "SCALARS pressure_on_solid double 1\n");
            fprintf(fd, "LOOKUP_TABLE default\n");
            for (k = 0; k < z; k++) {
                for (j = 0; j < y; j++) {
                    for (i = 0; i < x; i++) {
                        tmp = dB_map[i][j][k];
                        tmp1 = ntohll(*((long long int*) & tmp));
                        fwrite(&tmp1, 1, sizeof (double), fd);
                    }
                }
            }
            fclose(fd);
        }
    }

    if (advection == 1) {

        sprintf(filename, "%smarker_ic.vtk", output_dir);
        fd = fopen(filename, "wb");
        if (fd == NULL) {
            printf("-----------------------------------\n");
            printf("The output directory does not exist\n");
            printf("-----------------------------------\n");
            exit(13);
        }
        fprintf(fd, "# vtk DataFile Version 3.0\n");
        fprintf(fd, "TYCHO Marker File - Time:%g\n", time_sim / minute_in_sec);
        fprintf(fd, "BINARY\n");
        fprintf(fd, "DATASET RECTILINEAR_GRID\n");
        fprintf(fd, "DIMENSIONS %i %i %i\n", x, y, z);
        fprintf(fd, "X_COORDINATES %i double\n", x);
        for (i = 0; i < x; i++) {
            tmp = ((double) xmax / (double) x) * (double) i;
            tmp1 = ntohll(*((long long int*) & tmp));
            fwrite(&tmp1, 1, sizeof (double), fd);
        }
        fprintf(fd, "\nY_COORDINATES %i double\n", y);
        for (i = 0; i < y; i++) {
            tmp = ((double) ymax / (double) y) * (double) i;
            tmp1 = ntohll(*((long long int*) & tmp));
            fwrite(&tmp1, 1, sizeof (double), fd);
        }
        fprintf(fd, "\nZ_COORDINATES %i double\n", z);
        for (i = 0; i < z; i++) {
            tmp = ((double) zmax / (double) z) * (double) i;
            tmp1 = ntohll(*((long long int*) & tmp));
            fwrite(&tmp1, 1, sizeof (double), fd);
        }
        fprintf(fd, "\nPOINT_DATA %i\n", x * y * z);
        fprintf(fd, "SCALARS marker_density double 1\n");
        fprintf(fd, "LOOKUP_TABLE default\n");
        for (k = 0; k < z; k++) {
            for (j = 0; j < y; j++) {
                for (i = 0; i < x; i++) {
                    tmp = marker[i][j][k];
                    tmp1 = ntohll(*((long long int*) & tmp));
                    fwrite(&tmp1, 1, sizeof (double), fd);
                }
            }
        }
        fclose(fd);

    }

    return 0;
}
