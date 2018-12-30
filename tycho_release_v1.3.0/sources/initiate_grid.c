/*
 * initiate_grid.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 In this function the grid-parameters, such as spacing and physical coordinates for each cell are calculated
 */
int initiate_grid(int x, double xmin, double xmax, int y, double ymin, double ymax, int z, double zmin, double zmax) {
    int i, j, k;
    double dxfac, dyfac, dzfac;

    // spacing in all three directions
    dxfac = (xmax - xmin) / (double) x;
    dyfac = (ymax - ymin) / (double) y;
    dzfac = (zmax - zmin) / (double) z;

    // the coordinate at the cell boundary and the cell centre in x-direction
    for (i = 0; i < x; i++) {
        zxa[i] = xmin + (double) (i) * dxfac;
        zdx[i] = dxfac;
        zxc[i] = zxa[i] + 0.5 * zdx[i];
    }

    // the coordinate at the cell boundary and the cell centre in x-direction
    for (j = 0; j < y; j++) {
        zya[j] = ymin + (double) (j) * dyfac;
        zdy[j] = dyfac;
        zyc[j] = zya[j] + 0.5 * zdy[j];
    }

    // the coordiante at the cell boundary and the cell centre in x-direction
    for (k = 0; k < z; k++) {
        zza[k] = zmin + (double) (k) * dzfac;
        zdz[k] = dzfac;
        zzc[k] = zza[k] + 0.5 * zdz[k];
    }

    printf("Initiate grid done.\n");
    
    
    //A spacing in all directions
    spacing = dxfac;
    return 0;
}