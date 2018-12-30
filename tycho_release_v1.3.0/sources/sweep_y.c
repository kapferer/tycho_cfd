/*
 * sweep_y.c
 *
 * Author: Wolfgang Kapferer Wolfgang
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "variables_global.h"
#include "prototypes.h"

/*!
 Inline function to calculate the maximum of two doubles
 */
inline double max_sweep_y(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 The sweeping in all three directions (if 3D). This routine does the y-sweep and 
 is parallelized using OpenMP for shared memory machines.
 */
int sweep_y(int x, int y, int z, int flag) {
    int i, j, k, n, nmin, nmax;

    // direction setter 0==x, 1==y, 2==z
    int direction = 1;
    //1D variables
    double *rho_1D, *pre_1D, *eng_1D, *vx_1D, *vy_1D, *vz_1D, *marker_1D;
    //1D variables for viscosity
    double *rhodown, *rhoup, *rhofront, *rhoback;
    double *vxdown, *vxup, *vxfront, *vxback;
    double *vydown, *vyup, *vyfront, *vyback;
    double *vzdown, *vzup, *vzfront, *vzback;
    //for pressure calculations on the solid
    double *pressure_solid_1D;
    double *xa, *dx;
    double *dx0, *xa0;
    double *a_coef, *ai_coef, *b_coef, *bi_coef, *c_coef, *ci_coef;
    double *d_x;
    double *diffa;
    //parabola
    double *da, *ar;
    double *pl, *p6, *rl, *r6;
    double *u6, *ul, *vl, *v6;
    double *wl, *w6, *el, *e6;
    double *ql, *q6, *dp, *du;
    double *dv, *dw, *dq, *de;
    double *scratch1, *scratch2, *scratch3;
    double *dr, *deltaa;
    //for states
    double *plft, *prgh, *ulft, *urgh, *rlft, *rrgh, *Cdtdx, *fCdtdx;
    // for flattening
    double *steep, *flat;
    double **para;
    //riemann solver
    double *clft, *crgh, *plfti, *prghi, *pmid, *pmold, *wlft, *wrgh, *zlft, *zrgh;
    double *umidl, *umidr, *umid;
    //evolve
    double *dm, *dtbdm, *upmid, *xa1, *xa2, *xa3;
    double *vx_1D_old;
    double *dvol, *dvol0, *dvol1;
    //remap
    double *delta;
    double *fluxr, *fluxu, *fluxv, *fluxw, *fluxe, *fluxq;
    double *dm0;
    double *e_int_1D;

    dom_state state, state_next;

#ifdef _OPENMP
#pragma omp parallel default(none) \
                     private(i, j, k, n, nmin, nmax, state, state_next, rho_1D, pre_1D,\
                             eng_1D, vx_1D, vy_1D, vz_1D, marker_1D, \
                             pressure_solid_1D, dx0, xa0, xa, dx, \
                             a_coef, ai_coef, b_coef, bi_coef, c_coef, ci_coef, \
                             d_x, diffa, da, ar, pl, p6, rl, r6, \
                             u6, ul, vl, v6, wl, w6, \
                             el, e6, ql, q6, dp, du, dr, deltaa, \
                             dv, dw, dq, de, scratch1, scratch2, scratch3, \
                             plft, prgh, ulft, urgh, rlft, rrgh, Cdtdx, fCdtdx, \
                             steep, flat, para, clft, crgh, plfti, prghi, pmid, \
                             pmold, wlft, wrgh, zlft, zrgh, umidl, umidr, umid, \
                             dm, dtbdm, upmid, xa1, xa2, xa3, vx_1D_old, \
                             dvol, dvol0, dvol1, delta, fluxr, fluxu, fluxv, \
                             fluxw, fluxe, fluxq, dm0, e_int_1D, rhodown, rhoup, \
                             rhofront, rhoback, vxdown, vxup, vxfront, vxback, \
                             vydown, vyup, vyfront, vyback, vzdown, vzup, vzfront, \
                             vzback) \
                     shared(x, y, z, rho, dom,  \
                            pre, vx, vy, vz, marker, zya, zdy, \
                            pressure_on_solid, smallp, smallr, small, Gamma, max_array_length, flag, \
                            direction, obstacle_density, obstacle_temperature, advection, with_obstacles, \
                            viscosity_on_off, bound, inflow_density, inflow_velocity, dimension, \
                            rho_visc, vx_visc, vy_visc, vz_visc)
    {
#endif

        //==============================================================================
        init_ppm(&rho_1D, &pre_1D, &eng_1D, &vx_1D, &vy_1D, &vz_1D, &marker_1D, &dx0, &xa0, &xa, &dx,
                &a_coef, &ai_coef, &b_coef, &bi_coef, &c_coef, &ci_coef, &d_x, &da, &ar,
                &dp, &dr, &du, &pl, &p6, &rl, &r6, &ul, &u6, &vl, &v6, &wl, &w6, &el, &e6, &ql,
                &q6, &deltaa, &dv, &dw, &dq, &de, &scratch1, &scratch2, &scratch3, &diffa,
                &plft, &prgh, &ulft, &urgh, &rlft, &rrgh, &Cdtdx, &fCdtdx, &clft, &crgh,
                &plfti, &prghi, &pmid, &pmold, &wlft, &wrgh, &zlft, &zrgh, &umidl, &umidr,
                &umid, &dm, &dtbdm, &upmid, &xa1, &xa2, &xa3, &vx_1D_old, &e_int_1D, &dvol,
                &dvol0, &dvol1, &delta, &fluxr, &fluxu, &fluxv, &fluxw, &fluxe, &fluxq,
                &dm0, &steep, &flat, &para, &pressure_solid_1D, &rhodown, &rhoup, &rhofront,
                &rhoback, &vxdown, &vxup, &vxfront, &vxback, &vydown, &vyup, &vyfront,
                &vyback, &vzdown, &vzup, &vzfront, &vzback, dimension);

        //==============================================================================

        for (k = 0; k < z; k++) {

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif

            for (i = 0; i < x; i++) {

                int lefter;

                nmin = 0;
                nmax = 0;
                lefter = 0;

                for (j = 0; j < y; j++) {

                    n = j + 6 - nmin;

                    xa[n] = zya[j];
                    dx[n] = zdy[j];
                    xa0[n] = zya[j];
                    dx0[n] = zdy[j];

                    // Copy parts of a line from 3D domain into 1D domain
                    if (dom[i][j][k] == 0) {
                        state = DOM_FLUID;

                        rho_1D[n] = rho[i][j][k];
                        pre_1D[n] = pre[i][j][k];
                        vx_1D[n] = vy[i][j][k];
                        vy_1D[n] = vz[i][j][k];
                        vz_1D[n] = vx[i][j][k];

                        if (advection == 1) marker_1D[n] = marker[i][j][k];

                        if (with_obstacles == 1) pressure_solid_1D[n] = pressure_on_solid[i][j][k];

                        pre_1D[n] = max_sweep_y(smallp, pre_1D[n]);
                        eng_1D[n] = pre_1D[n] / (rho_1D[n] * Gamma) + 0.5 * ((pow(vx_1D[n], 2))+(pow(vy_1D[n], 2))+(pow(vz_1D[n], 2)));

                        // if viscosity is set on
                        if (viscosity_on_off == 1) {
                            if (dimension > 1) {
                                //====================================================
                                //first the i-1 case no obstacle
                                if ((i > 0) && (dom[i - 1][j][k] == 0) && (i < x - 1)) {
                                    rhodown[n] = rho[i - 1][j][k];
                                    vxdown[n] = vy[i - 1][j][k];
                                    vydown[n] = vz[i - 1][j][k];
                                    vzdown[n] = vx[i - 1][j][k];
                                }
                                if ((i > 0) && (dom[i - 1][j][k] != 0) && (i < x - 1)) {
                                    rhodown[n] = 0.0;
                                    vxdown[n] = 0.0;
                                    vydown[n] = 0.0;
                                    vzdown[n] = 0.0;
                                }
                                //at j==0 no obstacle is checked before entering viscosity part
                                if (i == 0) {
                                    if (bound.left == 0) rhodown[n] = rho_visc[i][j][k];
                                    if (bound.left == 1) rhodown[n] = rho_visc[i][j][k];
                                    if (bound.left == 2) rhodown[n] = small;
                                    if (bound.left == 3) rhodown[n] = rho_visc[i][j][k];
                                    if (bound.left == 4) rhodown[n] = inflow_density;
                                    if (bound.left == 5) rhodown[n] = rho_visc[x - 1][j][k];

                                    if (bound.left == 0) vxdown[n] = vy_visc[i][j][k];
                                    if (bound.left == 1) vxdown[n] = -vy_visc[i][j][k];
                                    if (bound.left == 2) vxdown[n] = small;
                                    if (bound.left == 3) vxdown[n] = vy_visc[i][j][k];
                                    if (bound.left == 4) vxdown[n] = 0.0;
                                    if (bound.left == 5) vxdown[n] = vy_visc[x - 1][j][k];

                                    if (bound.left == 0) vydown[n] = vz_visc[i][j][k];
                                    if (bound.left == 1) vydown[n] = -vz_visc[i][j][k];
                                    if (bound.left == 2) vydown[n] = small;
                                    if (bound.left == 3) {
                                        if (vy_visc[i][j][k] > 0.0) vydown[n] = small;
                                        if (vy_visc[i][j][k] <= 0.0) vydown[n] = vz_visc[i][j][k];
                                    }

                                    if (bound.left == 4) vydown[n] = 0.0;
                                    if (bound.left == 5) vydown[n] = vz_visc[x - 1][j][k];

                                    if (bound.left == 0) vzdown[n] = vx_visc[i][j][k];
                                    if (bound.left == 1) vzdown[n] = -vx_visc[i][j][k];
                                    if (bound.left == 2) vzdown[n] = small;
                                    if (bound.left == 3) vzdown[n] = vx_visc[i][j][k];
                                    if (bound.left == 4) vzdown[n] = inflow_velocity;
                                    if (bound.left == 5) vzdown[n] = vz_visc[x - 1][j][k];

                                    rhoup[n] = rho_visc[i + 1][j][k];
                                    vxup[n] = vx_visc[i + 1][j][k];
                                    vyup[n] = vy_visc[i + 1][j][k];
                                    vzup[n] = vz_visc[i + 1][j][k];
                                }
                                //====================================================

                                //====================================================
                                //now the case j = y - 1 no obstacle
                                if ((i < x - 1) && (dom[i + 1][j][k] == 0) && (i > 0)) {
                                    rhoup[n] = rho_visc[i + 1][j][k];
                                    vxup[n] = vy_visc[i + 1][j][k];
                                    vyup[n] = vz_visc[i + 1][j][k];
                                    vzup[n] = vx_visc[i + 1 ][j][k];
                                }
                                if ((i < x - 1) && (dom[i + 1][j][k] != 0) && (i > 0)) {
                                    rhoup[n] = 0.0;
                                    vxup[n] = 0.0;
                                    vyup[n] = 0.0;
                                    vzup[n] = 0.0;
                                }
                                //at j==0 no obstacle is checked before entering viscosity part
                                if (i == x - 1) {
                                    if (bound.right == 0) rhoup[n] = rho_visc[i][j][k];
                                    if (bound.right == 1) rhoup[n] = rho_visc[i][j][k];
                                    if (bound.right == 2) rhoup[n] = small;
                                    if (bound.right == 3) rhoup[n] = rho_visc[i][j][k];
                                    if (bound.right == 4) rhoup[n] = inflow_density;
                                    if (bound.right == 5) rhoup[n] = rho_visc[0][j][k];

                                    if (bound.right == 0) vxup[n] = vy_visc[i][j][k];
                                    if (bound.right == 1) vxup[n] = -vy_visc[i][j][k];
                                    if (bound.right == 2) vxup[n] = small;
                                    if (bound.right == 3) vxup[n] = vy_visc[i][j][k];
                                    if (bound.right == 4) vxup[n] = 0.0;
                                    if (bound.right == 5) vxup[n] = vx_visc[0][j][k];

                                    if (bound.right == 0) vyup[n] = vz_visc[i][j][k];
                                    if (bound.right == 1) vyup[n] = -vz_visc[i][j][k];
                                    if (bound.right == 2) vyup[n] = small;
                                    if (bound.right == 3) {
                                        if (vy_visc[i][j][k] < 0.0) vyup[n] = small;
                                        if (vy_visc[i][j][k] >= 0.0) vyup[n] = vz_visc[i][j][k];
                                    }
                                    if (bound.right == 4) vyup[n] = 0.0;
                                    if (bound.right == 5) vyup[n] = vz_visc[0][j][k];

                                    if (bound.right == 0) vzup[n] = vx_visc[i][j][k];
                                    if (bound.right == 1) vzup[n] = -vx_visc[i][j][k];
                                    if (bound.right == 2) vzup[n] = small;
                                    if (bound.right == 3) vzup[n] = vx_visc[i][j][k];
                                    if (bound.right == 4) vzup[n] = inflow_velocity;
                                    if (bound.right == 5) vzup[n] = vx_visc[0][j][k];

                                    rhodown[n] = rho_visc[i - 1][j][k];
                                    vxdown[n] = vx_visc[i - 1][j][k];
                                    vydown[n] = vy_visc[i - 1][j][k];
                                    vzdown[n] = vz_visc[i - 1][j][k];
                                }
                                //====================================================
                            }
                            if (dimension > 2) {
                                //====================================================
                                //now the k > 0 case no obstacle
                                if ((k > 0) && (dom[i][j][k - 1] == 0) && (k < z - 1)) {
                                    rhofront[n] = rho_visc[i][j][k - 1];
                                    vxfront[n] = vy_visc[i][j][k - 1];
                                    vyfront[n] = vz_visc[i][j][k - 1];
                                    vzfront[n] = vx_visc[i][j][k - 1];
                                }
                                if ((k > 0) && (dom[i][j][k - 1] != 0) && (k < z - 1)) {
                                    rhofront[n] = 0.0;
                                    vxfront[n] = 0.0;
                                    vyfront[n] = 0.0;
                                    vzfront[n] = 0.0;
                                }
                                //at k == 0 no obstacle is check-1ed before entering viscosity part
                                if (k == 0) {
                                    if (bound.front == 0) rhofront[n] = rho_visc[i][j][k];
                                    if (bound.front == 1) rhofront[n] = rho_visc[i][j][k];
                                    if (bound.front == 2) rhofront[n] = small;
                                    if (bound.front == 3) rhofront[n] = rho_visc[i][j][k];
                                    if (bound.front == 4) rhofront[n] = inflow_density;
                                    if (bound.front == 5) rhofront[n] = rho_visc[i][j][z - 1];

                                    if (bound.front == 0) vxfront[n] = vy_visc[i][j][k];
                                    if (bound.front == 1) vxfront[n] = -vy_visc[i][j][k];
                                    if (bound.front == 2) vxfront[n] = small;
                                    if (bound.front == 3) vxfront[n] = vy_visc[i][j][k];
                                    if (bound.front == 4) vxfront[n] = 0.0;
                                    if (bound.front == 5) vxfront[n] = vy_visc[i][j][z - 1];

                                    if (bound.front == 0) vyfront[n] = vz_visc[i][j][k];
                                    if (bound.front == 1) vyfront[n] = -vz_visc[i][j][k];
                                    if (bound.front == 2) vyfront[n] = small;
                                    if (bound.front == 3) vyfront[n] = vz_visc[i][j][k];
                                    if (bound.front == 4) vyfront[n] = 0.0;
                                    if (bound.front == 5) vyfront[n] = vz_visc[i][j][z - 1];

                                    if (bound.front == 0) vzfront[n] = vx_visc[i][j][k];
                                    if (bound.front == 1) vzfront[n] = -vx_visc[i][j][k];
                                    if (bound.front == 2) vzfront[n] = small;
                                    if (bound.front == 3) {
                                        if (vz_visc[i][j][k] > 0.0) vzfront[n] = small;
                                        if (vz_visc[i][j][k] <= 0.0) vzfront[n] = vx_visc[i][j][k];
                                    }
                                    if (bound.front == 4) vzfront[n] = inflow_velocity;
                                    if (bound.front == 5) vzfront[n] = vx_visc[i][j][z - 1];

                                    rhoback[n] = rho_visc[i][j][k + 1];
                                    vxback[n] = vx_visc[i][j][k + 1];
                                    vyback[n] = vy_visc[i][j][k + 1];
                                    vzback[n] = vz_visc[i][j][k + 1];
                                }
                                //====================================================

                                //====================================================
                                //now the k < z-1 case no obstacle
                                if ((k < z - 1) && (dom[i][j][k + 1] == 0) && (k > 0)) {
                                    rhoback[n] = rho_visc[i][j][k + 1];
                                    vxback[n] = vy_visc[i][j][k + 1];
                                    vyback[n] = vz_visc[i][j][k + 1];
                                    vzback[n] = vx_visc[i][j][k + 1];
                                }
                                if ((k < z - 1) && (dom[i][j][k + 1] != 0) && (k > 0)) {
                                    rhoback[n] = 0.0;
                                    vxback[n] = 0.0;
                                    vyback[n] = 0.0;
                                    vzback[n] = 0.0;
                                }
                                //at k == z - 1 no obstacle is check-1ed before entering viscosity part
                                if (k == z - 1) {
                                    if (bound.back == 0) rhoback[n] = rho_visc[i][j][k];
                                    if (bound.back == 1) rhoback[n] = rho_visc[i][j][k];
                                    if (bound.back == 2) rhoback[n] = small;
                                    if (bound.back == 3) rhoback[n] = rho_visc[i][j][k];
                                    if (bound.back == 4) rhoback[n] = inflow_density;
                                    if (bound.back == 5) rhoback[n] = rho_visc[i][j][0];

                                    if (bound.back == 0) vxback[n] = vy_visc[i][j][k];
                                    if (bound.back == 1) vxback[n] = -vy_visc[i][j][k];
                                    if (bound.back == 2) vxback[n] = small;
                                    if (bound.back == 3) vxback[n] = vy_visc[i][j][k];
                                    if (bound.back == 4) vxback[n] = 0.0;
                                    if (bound.back == 5) vxback[n] = vy_visc[i][j][0];

                                    if (bound.back == 0) vyback[n] = vz_visc[i][j][k];
                                    if (bound.back == 1) vyback[n] = -vz_visc[i][j][k];
                                    if (bound.back == 2) vyback[n] = small;
                                    if (bound.back == 3) vyback[n] = vz_visc[i][j][k];
                                    if (bound.back == 4) vyback[n] = 0.0;
                                    if (bound.back == 5) vyback[n] = vz_visc[i][j][0];

                                    if (bound.back == 0) vzback[n] = vx_visc[i][j][k];
                                    if (bound.back == 1) vzback[n] = -vx_visc[i][j][k];
                                    if (bound.back == 2) vzback[n] = small;
                                    if (bound.back == 3) {
                                        if (vz_visc[i][j][k] < 0.0) vzback[n] = small;
                                        if (vz_visc[i][j][k] >= 0.0) vzback[n] = vx_visc[i][j][k];
                                    }
                                    if (bound.back == 4) vzback[n] = inflow_velocity;
                                    if (bound.back == 5) vzback[n] = vx_visc[i][j][0];

                                    rhofront[n] = rho_visc[i][j][k - 1];
                                    vxfront[n] = vx_visc[i][j][k - 1];
                                    vyfront[n] = vy_visc[i][j][k - 1];
                                    vzfront[n] = vz_visc[i][j][k - 1];
                                }
                                //====================================================
                            }
                        }

                    } else if (dom[i][j][k] == 1) {
                        state = DOM_SOLID;

                        rho_1D[n] = obstacle_density;
                        pre_1D[n] = obstacle_temperature;
                        vx_1D[n] = 0.0;
                        vy_1D[n] = 0.0;
                        vz_1D[n] = 0.0;

                        if (advection == 1) marker_1D[n] = 0.0;

                        if (with_obstacles == 1) pressure_solid_1D[n] = 0.0;

                        pre_1D[n] = 0.0;
                        eng_1D[n] = 0.0;

                    }

                    if (j < y - 1) {
                        state_next = (dom[i][j + 1][k] == 0) ? DOM_FLUID : DOM_SOLID;
                    } else {
                        state_next = (state == DOM_FLUID) ? DOM_SOLID : DOM_FLUID;
                    }

                    if (state != state_next) {
                        if (state == DOM_FLUID) {
                            //*checks if we are at the most left or most right part of the domain
                            int bound_checker;
                            int jj;
                            int nminy2 = 6;
                            int nmaxy2 = (nmax - nmin) + 6;

                            //if the most left boundary is reached bound_checker is 0
                            //if the most right boundary is reached bound_checker is 0
                            if ((nmin == 0) || (nmax == (y - 1))) {
                                bound_checker = 0;
                            }

                            //Now we consider an obstacle in the row
                            if ((nmax - nmin) != (y - 1)) {
                                bound_checker = 1;
                                //is the left boundary at the edge of the computational domain
                                //lefter = 1
                                if (nmin == 0) lefter = 1;
                                //is the right boundary at the edge of the computational domain
                                //lefter =2;
                                if (nmax == (y - 1)) lefter = 2;
                                //if there is an obstacle left and right than
                                //lefter = 3
                                if ((nmin != 0) && (nmax != (y - 1))) lefter = 3;
                            }

                            ppm_step(i, j, k, direction, flag, nminy2, nmaxy2, a_coef, ai_coef, b_coef, bi_coef, c_coef, ci_coef,
                                    d_x, diffa, da, ar, pl, p6, rl, r6, u6, ul, vl, v6, wl, w6, el,
                                    e6, ql, q6, dp, du, dr, dv, dw, dq, de, scratch1, scratch2, scratch3,
                                    plft, prgh, ulft, urgh, rlft, rrgh, Cdtdx, fCdtdx, steep, flat,
                                    para, clft, crgh, plfti, prghi, pmid, pmold, wlft, wrgh, zlft,
                                    zrgh, umidl, umidr, umid, dm, dtbdm, upmid, xa1, xa2, xa3,
                                    vx_1D_old, dvol, dvol0, dvol1, delta, fluxr, fluxu, fluxv,
                                    fluxw, fluxe, fluxq, dm0, e_int_1D, rho_1D, pre_1D, eng_1D,
                                    vx_1D, vy_1D, vz_1D, marker_1D, pressure_solid_1D, dx0, xa0,
                                    xa, dx, bound_checker, lefter, rhodown, rhoup, rhofront,
                                    rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront,
                                    vyback, vzdown, vzup, vzfront, vzback, viscosity_on_off, dimension);

                            //put the solution back into the 3D arrays
                            for (jj = nmin; jj <= nmax; jj++) {
                                n = jj + 6 - nmin;

                                rho[i][jj][k] = rho_1D[n];
                                pre[i][jj][k] = pre_1D[n];
                                vy[i][jj][k] = vx_1D[n];
                                vz[i][jj][k] = vy_1D[n];
                                vx[i][jj][k] = vz_1D[n];

                                if (advection == 1) marker[i][jj][k] = marker_1D[n];

                                //for the pressure on the obstacle calculation
                                if (with_obstacles == 1) pressure_on_solid[i][jj][k] = pressure_solid_1D[n];
                            }

                        } else if (state == DOM_SOLID) {
                            nmin = nmax + 1;
                        }
                    }
                    nmax++;
                }
            }
        }
        //free section
        ppm_free(&rho_1D, &pre_1D, &eng_1D, &vx_1D, &vy_1D, &vz_1D, &marker_1D,
                &pressure_solid_1D, &dx0, &xa0, &xa, &dx, &a_coef, &ai_coef,
                &b_coef, &bi_coef, &c_coef, &ci_coef, &d_x, &da, &ar,
                &dp, &dr, &du, &pl, &p6, &rl, &r6, &ul, &u6, &vl, &v6, &wl, &w6, &el, &e6, &ql,
                &q6, &deltaa, &dv, &dw, &dq, &de, &scratch1, &scratch2, &scratch3, &diffa,
                &plft, &prgh, &ulft, &urgh, &rlft, &rrgh, &Cdtdx, &fCdtdx, &clft, &crgh,
                &plfti, &prghi, &pmid, &pmold, &wlft, &wrgh, &zlft, &zrgh, &umidl, &umidr,
                &umid, &dm, &dtbdm, &upmid, &xa1, &xa2, &xa3, &vx_1D_old, &e_int_1D, &dvol,
                &dvol0, &dvol1, &delta, &fluxr, &fluxu, &fluxv, &fluxw, &fluxe, &fluxq,
                &dm0, &steep, &flat, &para, &rhodown, &rhoup, &rhofront,
                &rhoback, &vxdown, &vxup, &vxfront, &vxback, &vydown, &vyup, &vyfront,
                &vyback, &vzdown, &vzup, &vzfront, &vzback, dimension);

#ifdef _OPENMP
    }
#endif

    return 0;
}