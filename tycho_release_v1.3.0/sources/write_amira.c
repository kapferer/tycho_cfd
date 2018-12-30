/*
 * write_amira.c
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
 Writes the initial model in a Amira Mesh format
 those binaries are in BigEndian format
 The slow order of the nested loops is mandatory here
 */
int write_amira(int x, int y, int z, int counter) {
    FILE *fd;
    char filename[600];
    int i, j, k;

    long long int tmp1, tmp1_vx, tmp1_vy, tmp1_vz;
    double tmp, tmp_vx, tmp_vy, tmp_vz;

    sprintf(filename, "%srho_%i.am", output_dir, counter);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# AmiraMesh BINARY 1.0\n");
    fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
    fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
    fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
    fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
    fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
    fprintf(fd, "}\n\nLattice { double ScalarField } = @1\n\n");
    fprintf(fd, "@1\n");
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

    sprintf(filename, "%spressure_%i.am", output_dir, counter);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# AmiraMesh BINARY 1.0\n");
    fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
    fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
    fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
    fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
    fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
    fprintf(fd, "}\n\nLattice { double ScalarField } = @1\n\n");
    fprintf(fd, "@1\n");
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

    sprintf(filename, "%stemp_%i.am", output_dir, counter);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# AmiraMesh BINARY 1.0\n");
    fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
    fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
    fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
    fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
    fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
    fprintf(fd, "}\n\nLattice { double ScalarField } = @1\n\n");
    fprintf(fd, "@1\n");
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

    sprintf(filename, "%svel_%i.am", output_dir, counter);
    fd = fopen(filename, "wb");
    if (fd == NULL) {
        printf("-----------------------------------\n");
        printf("The output directory does not exist\n");
        printf("-----------------------------------\n");
        exit(13);
    }
    fprintf(fd, "# AmiraMesh BINARY 1.0\n");
    fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
    fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
    fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
    fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
    fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
    fprintf(fd, "}\n\nLattice { double[3] ScalarField } = @1\n\n");
    fprintf(fd, "@1\n");
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
        sprintf(filename, "%spressure_on_solid_%i.am", output_dir, counter);
        fd = fopen(filename, "wb");
        if (fd == NULL) {
            printf("-----------------------------------\n");
            printf("The output directory does not exist\n");
            printf("-----------------------------------\n");
            exit(13);
        }
        fprintf(fd, "# AmiraMesh BINARY 1.0\n");
        fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
        fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
        fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
        fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
        fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
        fprintf(fd, "}\n\nLattice { double ScalarField } = @1\n\n");
        fprintf(fd, "@1\n");
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
            sprintf(filename, "%sdB_map_%i.am", output_dir, counter);
            fd = fopen(filename, "wb");
            if (fd == NULL) {
                printf("-----------------------------------\n");
                printf("The output directory does not exist\n");
                printf("-----------------------------------\n");
                exit(13);
            }
            fprintf(fd, "# AmiraMesh BINARY 1.0\n");
            fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
            fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
            fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
            fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
            fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
            fprintf(fd, "}\n\nLattice { double ScalarField } = @1\n\n");
            fprintf(fd, "@1\n");
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
        sprintf(filename, "%smarker_%i.am", output_dir, counter);
        fd = fopen(filename, "wb");
        if (fd == NULL) {
            printf("-----------------------------------\n");
            printf("The output directory does not exist\n");
            printf("-----------------------------------\n");
            exit(13);
        }
        fprintf(fd, "# AmiraMesh BINARY 1.0\n");
        fprintf(fd, "\n# Dimensions in x-, y-, and z-direction\n");
        fprintf(fd, "define Lattice %i %i %i\n\n", x, y, z);
        fprintf(fd, "Parameters {\n  CoordType \"uniform\"\n");
        fprintf(fd, "  # BoundingBox is xmin xmax ymin ymax zmin zmax\n");
        fprintf(fd, "  BoundingBox %g %g %g %g %g %g\n", xmin, xmax, ymin, ymax, zmin, zmax);
        fprintf(fd, "}\n\nLattice { double ScalarField } = @1\n\n");
        fprintf(fd, "@1\n");
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
