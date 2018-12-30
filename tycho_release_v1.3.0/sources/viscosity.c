/*
 * viscosity.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 Inline function to calculate the maximum of two doubles
 */
inline double max_viscosity(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 The velocity fluxes are changed due to viscosity
 Velocity-Flux = -nu * grad(v_x)
 This routine computes the velocity fluxes from
 viscosity and alters the velocity fluxes from the 
 hydro step.
 */
int viscosity(int i, int j, int k, int flag, int nmin, int nmax,
        double *rho_1D, double * vx_1D, double *vy_1D, double *vz_1D,
        double *pre_1D, double *e_int_1D, double *eng_1D, int lefter,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback,
        double *vxdown, double *vxup, double *vxfront, double *vxback,
        double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback,
        int dimension) {

    int n;
    double viscosity_term;
    double temperature_i, temperature_ii;
    double term_x, term_y, term_z;
    double term_x_1, term_y_1, term_z_1;
    double term_x_2, term_y_2, term_z_2;
    double dx, dm, dm0, dtbdx;
    double *dxflux, *dyflux, *dzflux;

    dm = 0.0;
    dm0 = 0.0;
    dx = 0.0;
    dtbdx = 0.0;

    dxflux = calloc((max_array_length + 12), sizeof dxflux);
    dyflux = calloc((max_array_length + 12), sizeof dyflux);
    dzflux = calloc((max_array_length + 12), sizeof dzflux);


    // spacing in all three directions
    if (flag == 0) dx = (xmax - xmin) / (double) x;
    if (flag == 1) dx = (ymax - ymin) / (double) y;
    if (flag == 2) dx = (zmax - zmin) / (double) z;

    for (n = nmin - 1; n <= nmax + 1; n++) {

        term_x = 0.0;
        term_y = 0.0;
        term_z = 0.0;

        term_x_1 = 0.0;
        term_y_1 = 0.0;
        term_z_1 = 0.0;

        term_x_2 = 0.0;
        term_y_2 = 0.0;
        term_z_2 = 0.0;

        //the core computation of viscosity
        temperature_i = pre_1D[n] / (gasconstant * rho_1D[n]);
        temperature_ii = pre_1D[n - 1] / (gasconstant * rho_1D[n - 1]);

        viscosity_term = -0.5 * (dynamic_viscosity(temperature_i) + dynamic_viscosity(temperature_ii));

        //the dxflux term----------------------------------------------------------------------------------
        term_x = -viscosity_term / dx * (rho_1D[n]*(4.0 / 3.0)*(vx_1D[n] - vx_1D[n - 1]));
        if (dimension > 1) {
            term_x_1 = (rhoup[n] + rhodown[n])*(2.0 / 3.0)*(vyup[n] - vydown[n])*0.5;
            term_x_1 += (rhoup[n - 1] + rhodown[n - 1])*(2.0 / 3.0)*(vyup[n - 1] - vydown[n - 1])*0.5;
        }
        if (dimension > 2) {
            term_x_2 = (rhoback[n] + rhofront[n])*(2.0 / 3.0)*(vzfront[n] - vzback[n])*0.5;
            term_x_2 += (rhoback[n - 1] + rhofront[n - 1])*(2.0 / 3.0)*(vzfront[n - 1] - vzback[n - 1])*0.5;
        }

        dxflux[n] = term_x - 0.25 * (term_x_1 + term_x_2);
        //the dxflux term----------------------------------------------------------------------------------

        //the dyflux term----------------------------------------------------------------------------------
        if (dimension > 1) {
            term_y = -viscosity_term / dx * (rho_1D[n]*(4.0 / 3.0)*(vy_1D[n] - vy_1D[n - 1]));
            term_y_1 = (rhoup[n] + rhodown[n])*(vxup[n] - vxdown[n])*0.5;
            term_y_1 += (rhoup[n - 1] + rhodown[n - 1])*(vxup[n - 1] - vxdown[n - 1])*0.5;

            dyflux[n] = term_y + 0.25 * term_y_1;
        }
        //the dyflux term----------------------------------------------------------------------------------

        //the dzflux term----------------------------------------------------------------------------------
        if (dimension > 2) {
            term_z = -viscosity_term / dx * (rho_1D[n]*(4.0 / 3.0)*(vz_1D[n] - vz_1D[n - 1]));
            term_z_1 = (rhoback[n] - rhofront[n])*(vxback[n] - vxfront[n])*0.5;
            term_z_1 += (rhoback[n - 1] + rhofront[n - 1])*(vxback[n - 1] - vxfront[n - 1])*0.5;

            dzflux[n] = term_z + 0.25 * term_z_1;
        }
        //the dzflux term----------------------------------------------------------------------------------
    }

    //Now the viscosity effect on the hydro quantities will be calculated
    for (n = nmin; n <= nmax; n++) {

        dm0 = 1 / rho_1D[n];
        dtbdx = dt / dx;

        vx_1D[n] = (vx_1D[n] * rho_1D[n] - dtbdx * (dxflux[n + 1] - dxflux[n])) * dm0;
        vy_1D[n] = (vy_1D[n] * rho_1D[n] - dtbdx * (dyflux[n + 1] - dyflux[n])) * dm0;
        vz_1D[n] = (vz_1D[n] * rho_1D[n] - dtbdx * (dzflux[n + 1] - dzflux[n])) * dm0;

    }

    free(dxflux);
    free(dyflux);
    free(dzflux);

    return 0;
}