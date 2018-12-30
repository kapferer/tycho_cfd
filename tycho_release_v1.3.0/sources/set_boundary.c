/*
 * set_boundary.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 The filling routine for the ghost cells according to the given boundary conditions.
 */
int set_boundary(int nmin, int nmax, int flag, int bound_checker, int lefter,
        double *rho_1D, double *eng_1D, double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, double *rhodown,
        double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback, double *vzdown,
        double *vzup, double *vzfront, double *vzback, int viscosity_on_off, int dimension) {

    //* The case with no obstacles in the computational domain
    if (bound_checker == 0) {
        if ((bound.down == 5) && (flag == 1))
            left_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.down == 4) && (flag == 1))
            left_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.down == 3) && (flag == 1))
            left_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.down == 2) && (flag == 1))
            left_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.down == 1) && (flag == 1))
            left_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.down == 0) && (flag == 1))
            left_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);

        if ((bound.up == 5) && (flag == 1))
            right_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.up == 4) && (flag == 1))
            right_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.up == 3) && (flag == 1))
            right_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.up == 2) && (flag == 1))
            right_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.up == 1) && (flag == 1))
            right_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.up == 0) && (flag == 1))
            right_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);

        if ((bound.left == 5) && (flag == 0))
            left_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.left == 4) && (flag == 0))
            left_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.left == 3) && (flag == 0))
            left_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.left == 2) && (flag == 0))
            left_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.left == 1) && (flag == 0))
            left_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.left == 0) && (flag == 0))
            left_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);

        if ((bound.right == 5) && (flag == 0))
            right_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.right == 4) && (flag == 0))
            right_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.right == 3) && (flag == 0))
            right_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.right == 2) && (flag == 0))
            right_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.right == 1) && (flag == 0))
            right_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.right == 0) && (flag == 0))
            right_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);

        if ((bound.front == 5) && (flag == 2))
            left_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.front == 4) && (flag == 2))
            left_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.front == 3) && (flag == 2))
            left_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.front == 2) && (flag == 2))
            left_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.front == 1) && (flag == 2))
            left_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.front == 0) && (flag == 2))
            left_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);

        if ((bound.back == 5) && (flag == 2))
            right_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.back == 4) && (flag == 2))
            right_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.back == 3) && (flag == 2))
            right_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.back == 2) && (flag == 2))
            right_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.back == 1) && (flag == 2))
            right_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
        if ((bound.back == 0) && (flag == 2))
            right_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                viscosity_on_off, dimension);
    } else {
        // domain with obstacles
        //the case at one side an obstacle at the other side a boundary of the computational domain is lefter 1 and 2
        //at both sides of the computational domain an obstacle is lefter 3
        if ((bound.down == 5) && (flag == 1) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
        }
        if ((bound.down == 4) && (flag == 1) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
        }
        if ((bound.down == 3) && (flag == 1) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
        }
        if ((bound.down == 2) && (flag == 1) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
        }
        if ((bound.down == 1) && (flag == 1) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
        }
        if ((bound.down == 0) && (flag == 1) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
        }

        if ((bound.up == 5) && (flag == 1) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup,
                    vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.up == 4) && (flag == 1) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.up == 3) && (flag == 1) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup,
                    vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.up == 2) && (flag == 1) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.up == 1) && (flag == 1) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.up == 0) && (flag == 1) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup,
                    vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }

        if ((bound.left == 5) && (flag == 0) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.left == 4) && (flag == 0) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup,
                    vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.left == 3) && (flag == 0) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag,
                    rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback,
                    viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.left == 2) && (flag == 0) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.left == 1) && (flag == 0) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.left == 0) && (flag == 0) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront,
                    rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off,
                    dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup,
                    vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }

        if ((bound.right == 5) && (flag == 0) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.right == 4) && (flag == 0) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.right == 3) && (flag == 0) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.right == 2) && (flag == 0) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.right == 1) && (flag == 0) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.right == 0) && (flag == 0) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }

        if ((bound.front == 5) && (flag == 2) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.front == 4) && (flag == 2) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.front == 3) && (flag == 2) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.front == 2) && (flag == 2) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.front == 1) && (flag == 2) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.front == 0) && (flag == 2) && (bound_checker == 1) && (lefter == 1)) {
            left_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }

        if ((bound.back == 5) && (flag == 2) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_periodic(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.back == 4) && (flag == 2) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_inflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup,
                    vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.back == 3) && (flag == 2) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_outflow(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.back == 2) && (flag == 2) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_small_padding(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback,
                    vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.back == 1) && (flag == 2) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_reflecting(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        if ((bound.back == 0) && (flag == 2) && (bound_checker == 1) && (lefter == 2)) {
            right_boundary_zero_gradient(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown,
                    vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
        //if both sides are belonging to an obstacle
        if ((bound_checker == 1) && (lefter == 3)) {
            right_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
            left_boundary_reflecting_on_obstacle(nmin, nmax, rho_1D, eng_1D, pre_1D, vx_1D, vy_1D, vz_1D,
                    pressure_solid_1D, xa0, dx0, xa, dx, flag, rhodown, rhoup, rhofront, rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback,
                    vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);
        }
    }


    return 0;
}

/*!
 Zero gradient boundary condition
 */
int left_boundary_zero_gradient(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup,
        double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback, int viscosity_on_off,
        int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmin - n] = rho_1D[nmin];
        pre_1D[nmin - n] = pre_1D[nmin];
        eng_1D[nmin - n] = eng_1D[nmin];

        vx_1D[nmin - n] = vx_1D[nmin];
        vy_1D[nmin - n] = vy_1D[nmin];
        vz_1D[nmin - n] = vz_1D[nmin];

        dx[nmin - n] = dx[nmin];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmin];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = rhodown[nmin];
                rhoup[nmin - n] = rhoup[nmin];
                vxdown[nmin - n] = vxdown[nmin];
                vxup[nmin - n] = vxup[nmin];
                vydown[nmin - n] = vydown[nmin];
                vyup[nmin - n] = vyup[nmin];
                vzdown[nmin - n] = vzdown[nmin];
                vzup[nmin - n] = vzup[nmin];
            }
            if (dimension > 2) {
                rhofront[nmin - n] = rhofront[nmin];
                rhoback[nmin - n] = rhoback[nmin];
                vxfront[nmin - n] = vxfront[nmin];
                vxback[nmin - n] = vxback[nmin];
                vyfront[nmin - n] = vyfront[nmin];
                vyback[nmin - n] = vyback[nmin];
                vzfront[nmin - n] = vzfront[nmin];
                vzback[nmin - n] = vzback[nmin];
            }
        }

    }

    return 0;
}

/*!
 Zero gradient boundary condition or if stratisfied atmosphere needed the
 atmosphere is extened into the ghost cells.
 */
int right_boundary_zero_gradient(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup,
        double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback, int viscosity_on_off,
        int dimension) {

    double T0, T, A, rho0, pre0;

    //T(z)=T_0 + A*(z-z0)
    T0 = 288.15; //[K]
    A = -6.5E-3; //[K/m]

    rho0 = 1.229; // Density at sea level [kg/m^3]
    pre0 = 1.013E5; //pressure at sea level [N/m^2]

    int n;

    if ((strat_const_atmos == 1) && (flag == 1)) {
        for (n = 1; n < 7; n++) {

            T = T0 + A * ((ymax / y) * (y + n - 1));

            rho_1D[nmax + n] = rho0 * pow((T / T0), (-grav_acc / (gasconstant * A) + 1));
            pre_1D[nmax + n] = pre0 * pow((T / T0), (-grav_acc / (gasconstant * A)));
            eng_1D[nmax + n] = eng_1D[nmax];

            vx_1D[nmax + n] = vx_1D[nmax];
            vy_1D[nmax + n] = vy_1D[nmax];
            vz_1D[nmax + n] = vz_1D[nmax];

            dx[nmax + n] = dx[nmax];
            xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
            dx0[nmax + n] = dx0[nmax];
            xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

            if (viscosity_on_off == 1) {
                if (dimension > 1) {
                    rhodown[nmax + n] = rhodown[nmax];
                    rhoup[nmax + n] = rhoup[nmax];
                    vxdown[nmax + n] = vxdown[nmax];
                    vxup[nmax + n] = vxup[nmax];
                    vydown[nmax + n] = vydown[nmax];
                    vyup[nmax + n] = vyup[nmax];
                    vzdown[nmax + n] = vzdown[nmax];
                    vzup[nmax + n] = vzup[nmax];
                }
                if (dimension > 2) {
                    rhofront[nmax + n] = rhofront[nmax];
                    rhoback[nmax + n] = rhoback[nmax];
                    vxfront[nmax + n] = vxfront[nmax];
                    vxback[nmax + n] = vxback[nmax];
                    vyfront[nmax + n] = vyfront[nmax];
                    vyback[nmax + n] = vyback[nmax];
                    vzfront[nmax + n] = vzfront[nmax];
                    vzback[nmax + n] = vzback[nmax];
                }
            }


        }
    } else {
        for (n = 1; n < 7; n++) {
            rho_1D[nmax + n] = rho_1D[nmax];
            pre_1D[nmax + n] = pre_1D[nmax];
            eng_1D[nmax + n] = eng_1D[nmax];

            vx_1D[nmax + n] = vx_1D[nmax];
            vy_1D[nmax + n] = vy_1D[nmax];
            vz_1D[nmax + n] = vz_1D[nmax];

            dx[nmax + n] = dx[nmax];
            xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
            dx0[nmax + n] = dx0[nmax];
            xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

            if (viscosity_on_off == 1) {
                if (dimension > 1) {
                    rhodown[nmax + n] = rhodown[nmax];
                    rhoup[nmax + n] = rhoup[nmax];
                    vxdown[nmax + n] = vxdown[nmax];
                    vxup[nmax + n] = vxup[nmax];
                    vydown[nmax + n] = vydown[nmax];
                    vyup[nmax + n] = vyup[nmax];
                    vzdown[nmax + n] = vzdown[nmax];
                    vzup[nmax + n] = vzup[nmax];
                }
                if (dimension > 2) {
                    rhofront[nmax + n] = rhofront[nmax];
                    rhoback[nmax + n] = rhoback[nmax];
                    vxfront[nmax + n] = vxfront[nmax];
                    vxback[nmax + n] = vxback[nmax];
                    vyfront[nmax + n] = vyfront[nmax];
                    vyback[nmax + n] = vyback[nmax];
                    vzfront[nmax + n] = vzfront[nmax];
                    vzback[nmax + n] = vzback[nmax];
                }
            }
        }
    }


    return 0;
}

/*!
 Reflecting boundary condition
 */
int left_boundary_reflecting(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup,
        double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback, int viscosity_on_off,
        int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmin - n] = rho_1D[nmin + n - 1];
        pre_1D[nmin - n] = pre_1D[nmin + n - 1];
        eng_1D[nmin - n] = eng_1D[nmin + n - 1];

        vx_1D[nmin - n] = -1 * vx_1D[nmin + n - 1];
        vy_1D[nmin - n] = -1 * vy_1D[nmin + n - 1];
        vz_1D[nmin - n] = -1 * vz_1D[nmin + n - 1];

        dx[nmin - n] = dx[nmin + n - 1];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmin + n - 1];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = rhodown[nmin + n - 1];
                rhoup[nmin - n] = rhoup[nmin + n - 1];
                vxdown[nmin - n] = -vxdown[nmin + n - 1];
                vxup[nmin - n] = -vxup[nmin + n - 1];
                vydown[nmin - n] = -vydown[nmin + n - 1];
                vyup[nmin - n] = -vyup[nmin + n - 1];
                vzdown[nmin - n] = -vzdown[nmin + n - 1];
                vzup[nmin - n] = -vzup[nmin + n - 1];
            }
            if (dimension > 2) {
                rhofront[nmin - n] = rhofront[nmin + n - 1];
                rhoback[nmin - n] = rhoback[nmin + n - 1];
                vxfront[nmin - n] = -vxfront[nmin + n - 1];
                vxback[nmin - n] = -vxback[nmin + n - 1];
                vyfront[nmin - n] = -vyfront[nmin + n - 1];
                vyback[nmin - n] = -vyback[nmin + n - 1];
                vzfront[nmin - n] = -vzfront[nmin + n - 1];
                vzback[nmin - n] = -vzback[nmin + n - 1];
            }
        }
    }

    return 0;
}

/*!
 Reflecting boundary condition on obstacle
 */
int left_boundary_reflecting_on_obstacle(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n, k;
    k = 1;

    pressure_solid_1D[nmin] = pre_1D[nmin];

    for (n = 1; n < 7; n++) {

        rho_1D[nmin - n] = rho_1D[nmin + n - 1];
        pre_1D[nmin - n] = pre_1D[nmin + n - 1];
        eng_1D[nmin - n] = eng_1D[nmin + n - 1];

        vx_1D[nmin - n] = -1 * vx_1D[nmin + n - 1];
        vy_1D[nmin - n] = -1 * vy_1D[nmin + n - 1];
        vz_1D[nmin - n] = -1 * vz_1D[nmin + n - 1];

        dx[nmin - n] = dx[nmin + n - 1];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmin + n - 1];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = rhodown[nmin + n - 1];
                rhoup[nmin - n] = rhoup[nmin + n - 1];
                vxdown[nmin - n] = -vxdown[nmin + n - 1];
                vxup[nmin - n] = -vxup[nmin + n - 1];
                vydown[nmin - n] = -vydown[nmin + n - 1];
                vyup[nmin - n] = -vyup[nmin + n - 1];
                vzdown[nmin - n] = -vzdown[nmin + n - 1];
                vzup[nmin - n] = -vzup[nmin + n - 1];
            }
            if (dimension > 2) {
                rhofront[nmin - n] = rhofront[nmin + n - 1];
                rhoback[nmin - n] = rhoback[nmin + n - 1];
                vxfront[nmin - n] = -vxfront[nmin + n - 1];
                vxback[nmin - n] = -vxback[nmin + n - 1];
                vyfront[nmin - n] = -vyfront[nmin + n - 1];
                vyback[nmin - n] = -vyback[nmin + n - 1];
                vzfront[nmin - n] = -vzfront[nmin + n - 1];
                vzback[nmin - n] = -vzback[nmin + n - 1];
            }
        }
    }

    return 0;
}

/*!
 Reflecting boundary condition
 */
int right_boundary_reflecting(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup,
        double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback, int viscosity_on_off,
        int dimension) {

    int n;

    for (n = 1; n < 7; n++) {

        rho_1D[n + nmax] = rho_1D[nmax + 1 - n];
        pre_1D[n + nmax] = pre_1D[nmax + 1 - n];
        eng_1D[n + nmax] = eng_1D[nmax + 1 - n];

        vx_1D[n + nmax] = -vx_1D[nmax + 1 - n];
        vy_1D[n + nmax] = -vy_1D[nmax + 1 - n];
        vz_1D[n + nmax] = -vz_1D[nmax + 1 - n];

        dx[n + nmax] = dx[nmax + 1 - n];
        xa[n + nmax] = xa[nmax + n - 1] + dx[nmax + n - 1];
        dx0[n + nmax] = dx0[nmax + 1 - n];
        xa0[n + nmax] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[n + nmax] = rhodown[nmax + 1 - n];
                rhoup[n + nmax] = rhoup[nmax + 1 - n];
                vxdown[n + nmax] = -vxdown[nmax + 1 - n];
                vxup[n + nmax] = -vxup[nmax + 1 - n];
                vydown[n + nmax] = -vydown[nmax + 1 - n];
                vyup[n + nmax] = -vyup[nmax + 1 - n];
                vzdown[n + nmax] = -vzdown[nmax + 1 - n];
                vzup[n + nmax] = -vzup[nmax + 1 - n];
            }
            if (dimension > 2) {
                rhofront[n + nmax] = rhofront[nmax + 1 - n];
                rhoback[n + nmax] = rhoback[nmax + 1 - n];
                vxfront[n + nmax] = -vxfront[nmax + 1 - n];
                vxback[n + nmax] = -vxback[nmax + 1 - n];
                vyfront[n + nmax] = -vyfront[nmax + 1 - n];
                vyback[n + nmax] = -vyback[nmax + 1 - n];
                vzfront[n + nmax] = -vzfront[nmax + 1 - n];
                vzback[n + nmax] = -vzback[nmax + 1 - n];
            }
        }
    }

    return 0;
}

/*!
 Reflecting boundary condition on obstacle
 */
int right_boundary_reflecting_on_obstacle(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    pressure_solid_1D[nmax] = pre_1D[nmax];


    for (n = 1; n < 7; n++) {

        rho_1D[n + nmax] = rho_1D[nmax + 1 - n];
        pre_1D[n + nmax] = pre_1D[nmax + 1 - n];
        eng_1D[n + nmax] = eng_1D[nmax + 1 - n];

        vx_1D[n + nmax] = -vx_1D[nmax + 1 - n];
        vy_1D[n + nmax] = -vy_1D[nmax + 1 - n];
        vz_1D[n + nmax] = -vz_1D[nmax + 1 - n];

        dx[n + nmax] = dx[nmax + 1 - n];
        xa[n + nmax] = xa[nmax + n - 1] + dx[nmax + n - 1];
        dx0[n + nmax] = dx0[nmax + 1 - n];
        xa0[n + nmax] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmax + n] = rhodown[nmax + 1 - n];
                rhoup[nmax + n] = rhoup[nmax + 1 - n];
                vxdown[nmax + n] = -vxdown[nmax + 1 - n];
                vxup[nmax + n] = -vxup[nmax + 1 - n];
                vydown[nmax + n] = -vydown[nmax + 1 - n];
                vyup[nmax + n] = -vyup[nmax + 1 - n];
                vzdown[nmax + n] = -vzdown[nmax + 1 - n];
                vzup[nmax + n] = -vzup[nmax + 1 - n];
            }
            if (dimension > 2) {
                rhofront[nmax + n] = rhofront[nmax + 1 - n];
                rhoback[nmax + n] = rhoback[nmax + 1 - n];
                vxfront[nmax + n] = -vxfront[nmax + 1 - n];
                vxback[nmax + n] = -vxback[nmax + 1 - n];
                vyfront[nmax + n] = -vyfront[nmax + 1 - n];
                vyback[nmax + n] = -vyback[nmax + 1 - n];
                vzfront[nmax + n] = -vzfront[nmax + 1 - n];
                vzback[nmax + n] = -vzback[nmax + 1 - n];
            }
        }
    }

    return 0;
}

/*!
 Small number (e.g. 1E-50) ghost cell filling
 */
int left_boundary_small_padding(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown,
        double *vxup, double *vxfront, double *vxback, double *vydown, double *vyup,
        double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmin - n] = small;
        pre_1D[nmin - n] = small;
        eng_1D[nmin - n] = small;

        vx_1D[nmin - n] = small;
        vy_1D[nmin - n] = small;
        vz_1D[nmin - n] = small;

        dx[nmin - n] = dx[nmin];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmin];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = small;
                rhoup[nmin - n] = small;
                vxdown[nmin - n] = small;
                vxup[nmin - n] = small;
                vydown[nmin - n] = small;
                vyup[nmin - n] = small;
                vzdown[nmin - n] = small;
                vzup[nmin - n] = small;
            }
            if (dimension > 2) {
                rhofront[nmin - n] = small;
                rhoback[nmin - n] = small;
                vxfront[nmin - n] = small;
                vxback[nmin - n] = small;
                vyfront[nmin - n] = small;
                vyback[nmin - n] = small;
                vzfront[nmin - n] = small;
                vzback[nmin - n] = small;
            }
        }
    }

    return 0;
}

/*!
 Small number (e.g. 1E-50) ghost cell filling
 */
int right_boundary_small_padding(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[n + nmax] = small;
        pre_1D[n + nmax] = small;
        eng_1D[n + nmax] = small;

        vx_1D[n + nmax] = small;
        vy_1D[n + nmax] = small;
        vz_1D[n + nmax] = small;

        dx[nmax + n] = dx[nmax];
        xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
        dx0[nmax + n] = dx0[nmax];
        xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmax + n] = small;
                rhoup[nmax + n] = small;
                vxdown[nmax + n] = small;
                vxup[nmax + n] = small;
                vydown[nmax + n] = small;
                vyup[nmax + n] = small;
                vzdown[nmax + n] = small;
                vzup[nmax + n] = small;
            }
            if (dimension > 2) {
                rhofront[nmax + n] = small;
                rhoback[nmax + n] = small;
                vxfront[nmax + n] = small;
                vxback[nmax + n] = small;
                vyfront[nmax + n] = small;
                vyback[nmax + n] = small;
                vzfront[nmax + n] = small;
                vzback[nmax + n] = small;
            }
        }
    }

    return 0;
}

/*!
 outflow boundary condition
 */
int left_boundary_outflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmin - n] = rho_1D[nmin];
        pre_1D[nmin - n] = pre_1D[nmin];


        if (vx_1D[nmin] > 0) {
            vx_1D[nmin - n] = small;
            eng_1D[nmin - n] = pre_1D[nmin] / (rho_1D[nmin] * Gamma) + 0.5 *
                    (pow(vx_1D[nmin], 2) + pow(vy_1D[nmin], 2) + pow(vz_1D[nmin], 2));
        } else {
            vx_1D[nmin - n] = vx_1D[nmin];
            eng_1D[nmin - n] = eng_1D[nmin];
        }

        dx[nmin - n] = dx[nmin];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmin];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = rhodown[nmin];
                rhoup[nmin - n] = rhoup[nmin];
                vxdown[nmin - n] = vxdown[nmin];
                vxup[nmin - n] = small;
                vydown[nmin - n] = small;
                vyup[nmin - n] = vyup[nmin];
                vzdown[nmin - n] = vzdown[nmin];
                vzup[nmin - n] = vzup[nmin];
            }
            if (dimension > 2) {
                rhofront[nmin - n] = rhofront[nmin];
                rhoback[nmin - n] = rhoback[nmin];
                vxfront[nmin - n] = small;
                vxback[nmin - n] = small;
                vyfront[nmin - n] = vyfront[nmin];
                vyback[nmin - n] = vyback[nmin];
                vzfront[nmin - n] = vzfront[nmin];
                vzback[nmin - n] = vzback[nmin];
            }
        }
    }

    return 0;
}

/*!
 outflow boundary condition
 */
int right_boundary_outflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmax + n] = rho_1D[nmax];
        pre_1D[nmax + n] = pre_1D[nmax];

        if (vx_1D[nmax] < 0) {
            vx_1D[nmax + n] = small;
            eng_1D[nmax + n] = pre_1D[nmax] / (rho_1D[nmax] * Gamma) + 0.5 *
                    (pow(vx_1D[nmax], 2) + pow(vy_1D[nmax], 2) + pow(vz_1D[nmax], 2));

        } else {
            vx_1D[nmax + n] = vx_1D[nmax];
            eng_1D[nmax + n] = eng_1D[nmax];

        }

        dx[nmax + n] = dx[nmax];
        xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
        dx0[nmax + n] = dx0[nmax];
        xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmax + n] = rhodown[nmax];
                rhoup[nmax + n] = rhoup[nmax];
                vxdown[nmax + n] = small;
                vxup[nmax + n] = small;
                vyup[nmax + n] = vyup[nmax];
                vzdown[nmax + n] = vzdown[nmax];
                vzup[nmax + n] = vzup[nmax];
            }
            if (dimension > 2) {
                rhofront[nmax + n] = rhofront[nmax];
                rhoback[nmax + n] = rhoback[nmax];
                vxfront[nmax + n] = small;
                vxback[nmax + n] = small;
                vyfront[nmax + n] = vyfront[nmax];
                vyback[nmax + n] = vyback[nmax];
                vzfront[nmax + n] = vzfront[nmax];
                vzback[nmax + n] = vzback[nmax];
            }
        }
    }


    return 0;
}

/*!
 outflow boundary condition
 */
int left_boundary_inflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmin - n] = inflow_density;
        pre_1D[nmin - n] = inflow_temperature * gasconstant*inflow_density;
        eng_1D[nmin - n] = (inflow_temperature * gasconstant * inflow_density)
                / (inflow_density * Gamma) + 0.5 * (pow(inflow_velocity, 2));

        vx_1D[nmin - n] = inflow_velocity;
        vy_1D[nmin - n] = 0.0;
        vz_1D[nmin - n] = 0.0;

        dx[nmin - n] = dx[nmin];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmin];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = inflow_density;
                rhoup[nmin - n] = inflow_density;
                vxdown[nmin - n] = inflow_velocity;
                vxup[nmin - n] = inflow_velocity;
                vydown[nmin - n] = 0.0;
                vyup[nmin - n] = 0.0;
                vzdown[nmin - n] = 0.0;
                vzup[nmin - n] = 0.0;
            }
            if (dimension > 2) {
                rhofront[nmin - n] = inflow_density;
                rhoback[nmin - n] = inflow_density;
                vxfront[nmin - n] = inflow_velocity;
                vxback[nmin - n] = inflow_velocity;
                vyfront[nmin - n] = 0.0;
                vyback[nmin - n] = 0.0;
                vzfront[nmin - n] = 0.0;
                vzback[nmin - n] = 0.0;
            }
        }
    }

    return 0;
}

/*!
 outflow boundary condition
 */
int right_boundary_inflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmax + n] = inflow_density;
        pre_1D[nmax + n] = inflow_temperature * gasconstant*inflow_density;
        eng_1D[nmax + n] = (inflow_temperature * gasconstant * inflow_density)
                / (inflow_density * Gamma) + 0.5 * (pow(inflow_velocity, 2));

        vx_1D[nmax + n] = inflow_velocity;
        vy_1D[nmax + n] = 0.0;
        vz_1D[nmax + n] = 0.0;

        dx[nmax + n] = dx[nmax];
        xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
        dx0[nmax + n] = dx0[nmax];
        xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmax + n] = inflow_density;
                rhoup[nmax + n] = inflow_density;
                vxdown[nmax + n] = inflow_velocity;
                vxup[nmax + n] = inflow_velocity;
                vydown[nmax + n] = 0.0;
                vyup[nmax + n] = 0.0;
                vzdown[nmax + n] = 0.0;
                vzup[nmax + n] = 0.0;
            }
            if (dimension > 2) {
                rhofront[nmax + n] = inflow_density;
                rhoback[nmax + n] = inflow_density;
                vxfront[nmax + n] = inflow_velocity;
                vxback[nmax + n] = inflow_velocity;
                vyfront[nmax + n] = 0.0;
                vyback[nmax + n] = 0.0;
                vzfront[nmax + n] = 0.0;
                vzback[nmax + n] = 0.0;
            }
        }
    }


    return 0;
}

/*!
 periodic boundary conditions
 */
int left_boundary_periodic(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup,
        double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback, int viscosity_on_off,
        int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmin - n] = rho_1D[nmax + 1 - n];
        pre_1D[nmin - n] = pre_1D[nmax + 1 - n];
        eng_1D[nmin - n] = eng_1D[nmax + 1 - n];

        vx_1D[nmin - n] = vx_1D[nmax + 1 - n];
        vy_1D[nmin - n] = vy_1D[nmax + 1 - n];
        vz_1D[nmin - n] = vz_1D[nmax + 1 - n];


        dx[nmin - n] = dx[nmax + 1 - n];
        xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
        dx0[nmin - n] = dx0[nmax + 1 - n];
        xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmin - n] = rhodown[nmax + 1 - n];
                rhoup[nmin - n] = rhoup[nmax + 1 - n];
                vxdown[nmin - n] = vxdown[nmax + 1 - n];
                vxup[nmin - n] = vxup[nmax + 1 - n];
                vydown[nmin - n] = vydown[nmax + 1 - n];
                vyup[nmin - n] = vyup[nmax + 1 - n];
                vzdown[nmin - n] = vzdown[nmax + 1 - n];
                vzup[nmin - n] = vzup[nmax + 1 - n];
            }
            if (dimension > 2) {
                rhofront[nmin - n] = rhofront[nmax + 1 - n];
                rhoback[nmin - n] = rhoback[nmax + 1 - n];
                vxfront[nmin - n] = vxfront[nmax + 1 - n];
                vxback[nmin - n] = vxback[nmax + 1 - n];
                vyfront[nmin - n] = vyfront[nmax + 1 - n];
                vyback[nmin - n] = vyback[nmax + 1 - n];
                vzfront[nmin - n] = vzfront[nmax + 1 - n];
                vzback[nmin - n] = vzback[nmax + 1 - n];
            }
        }
    }


    return 0;
}

/*!
 periodic boundary conditions
 */
int right_boundary_periodic(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup,
        double *vxfront, double *vxback, double *vydown, double *vyup, double *vyfront,
        double *vyback, double *vzdown, double *vzup, double *vzfront, double *vzback,
        int viscosity_on_off, int dimension) {

    int n;

    for (n = 1; n < 7; n++) {
        rho_1D[nmax + n] = rho_1D[nmin + n - 1];
        pre_1D[nmax + n] = pre_1D[nmin + n - 1];
        eng_1D[nmax + n] = eng_1D[nmin + n - 1];

        vx_1D[nmax + n] = vx_1D[nmin + n - 1];
        vy_1D[nmax + n] = vy_1D[nmin + n - 1];
        vz_1D[nmax + n] = vz_1D[nmin + n - 1];

        dx[nmax + n] = dx[nmin + n - 1];
        xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
        dx0[nmax + n] = dx0[nmin + n - 1];
        xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];

        if (viscosity_on_off == 1) {
            if (dimension > 1) {
                rhodown[nmax + n] = rhodown[nmin + n - 1];
                rhoup[nmax + n] = rhoup[nmin + n - 1];
                vxdown[nmax + n] = vxdown[nmin + n - 1];
                vxup[nmax + n] = vxup[nmin + n - 1];
                vydown[nmax + n] = vydown[nmin + n - 1];
                vyup[nmax + n] = vyup[nmin + n - 1];
                vzdown[nmax + n] = vzdown[nmin + n - 1];
                vzup[nmax + n] = vzup[nmin + n - 1];
            }
            if (dimension > 2) {
                rhofront[nmax + n] = rhofront[nmin + n - 1];
                rhoback[nmax + n] = rhoback[nmin + n - 1];
                vxfront[nmax + n] = vxfront[nmin + n - 1];
                vxback[nmax + n] = vxback[nmin + n - 1];
                vyfront[nmax + n] = vyfront[nmin + n - 1];
                vyback[nmax + n] = vyback[nmin + n - 1];
                vzfront[nmax + n] = vzfront[nmin + n - 1];
                vzback[nmax + n] = vzback[nmin + n - 1];
            }
        }
    }

    return 0;
}
