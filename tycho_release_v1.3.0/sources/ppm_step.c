/*
 * ppm_step.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 * The inner PPM Hydro scheme. First the relevant 1D distributions from the computational
 * domain is given to this function. The scheme is as followed: First fill the ghost cells according
 * to the defined boundary conditions. Solve the diffusion equation if needed.
 * Calculate the parabolic function coefficients, then the volumes and the flatten coefficients
 * due toe the steepness of the neighboring cells. Now the parabolic interpolated
 * quantities are determined. The quantities at the interpolated points are given.
 * The Riemann solver solves now the shock-tube problem. Then the hydrodynamic 
 * quantities are evolved and finally re-mapped on the computational grid.
 */
int ppm_step(int i, int j, int k, int direction, int flag, int nmin, int nmax, double *a_coef, double *ai_coef,
        double *b_coef, double *bi_coef, double *c_coef, double *ci_coef,
        double *d_x, double *diffa, double *da, double *ar, double *pl,
        double *p6, double *rl, double *r6, double *u6, double *ul,
        double *vl, double *v6, double *wl, double *w6, double *el,
        double *e6, double *ql, double *q6, double *dp, double *du,
        double *dr, double *dv, double *dw, double *dq, double *de,
        double *scratch1, double *scratch2, double *scratch3,
        double *plft, double *prgh, double *ulft, double *urgh,
        double *rlft, double *rrgh, double *Cdtdx, double *fCdtdx,
        double *steep, double *flat, double **para, double *clft,
        double *crgh, double *plfti, double *prghi, double *pmid,
        double *pmold, double *wlft, double *wrgh, double *zlft,
        double *zrgh, double *umidl, double *umidr, double *umid,
        double *dm, double *dtbdm, double *upmid, double *xa1, double *xa2, double *xa3,
        double *vx_1D_old, double *dvol, double *dvol0, double *dvol1,
        double *delta, double *fluxr, double *fluxu, double *fluxv, double *fluxw,
        double *fluxe, double *fluxq, double *dm0, double *e_int_1D, double *rho_1D,
        double *pre_1D, double *eng_1D, double *vx_1D, double *vy_1D,
        double *vz_1D, double *marker_1D, double *pressure_solid_1D,
        double *dx0, double *xa0, double *xa, double *dx, int bound_checker, int lefter,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown,
        double *vxup, double *vxfront, double *vxback, double *vydown, double *vyup,
        double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension) {

    int n;

    set_boundary(nmin, nmax, flag, bound_checker, lefter, rho_1D, eng_1D, pre_1D, vx_1D,
            vy_1D, vz_1D, pressure_solid_1D, xa0, dx0, xa, dx, rhodown, rhoup, rhofront,
            rhoback, vxdown, vxup, vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown,
            vzup, vzfront, vzback, viscosity_on_off, dimension);

    para_coef(nmin - 4, nmax + 5, flag, a_coef, dx,
            ai_coef, b_coef, bi_coef, c_coef,
            d_x, para, ci_coef);

    volume(nmin, nmax, dvol, dx, dvol0, dx0);

    flatten(nmin, nmax, pre_1D, vx_1D, steep, flat);

    parabola(nmin - 4, nmax + 4, pre_1D, dp, p6, pl, 0, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 4, nmax + 4, rho_1D, dr, r6, rl, 1, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);
    parabola(nmin - 4, nmax + 4, vx_1D, du, u6, ul, 2, diffa,
            da, para, ar, flat, scratch1, scratch2, scratch3);

    states(nmin, nmax, flag, pre_1D, rho_1D, dx, Cdtdx, fCdtdx, plft,
            pl, dp, p6, ulft, ul, du, u6, rlft, rl, dr,
            r6, prgh, urgh, rrgh);

    riemann(nmin - 3, nmax + 4, clft, crgh, rlft, rrgh, plfti,
            prghi, pmid, pmold, plft, wrgh, prgh, wlft,
            zlft, zrgh, umidl, umidr, umid, urgh, ulft);

    evolve(nmin, nmax, flag, rho_1D, dvol, dm, dtbdm,
            xa1, xa, dvol1, umid, upmid, pmid, xa2,
            dx, xa3, vx_1D_old, vx_1D, vy_1D, vz_1D,
            eng_1D, e_int_1D, pre_1D);

    // if viscosity is set on
    if (viscosity_on_off == 1) {
        viscosity(i, j, k, flag, nmin, nmax, rho_1D, vx_1D, vy_1D, vz_1D, pre_1D,
                e_int_1D, eng_1D, lefter, rhodown, rhoup, rhofront, rhoback, vxdown, vxup,
                vxfront, vxback, vydown, vyup, vyfront, vyback, vzdown, vzup,
                vzfront, vzback, dimension);
    }

    remap(nmin, nmax, flag, a_coef, dx,
            ai_coef, b_coef, bi_coef, c_coef,
            d_x, para, ci_coef, dr, r6, rl, diffa,
            da, ar, flat, scratch1, scratch2, scratch3,
            du, u6, ul, dv, v6, vl, w6, wl,
            dq, q6, ql, de, e6, el, xa, xa0, delta,
            fluxr, fluxu, fluxv, fluxw, dw, fluxe, fluxq,
            dm, rho_1D, dvol, dvol0, dm0, vx_1D, vy_1D,
            vz_1D, eng_1D, e_int_1D, pre_1D);


    if (advection == 1) {
        advect(nmin, nmax, flag, dx, vx_1D, marker_1D);
    }


    return 0;
}
