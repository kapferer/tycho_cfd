/*
 * init.c
 *
 * Author: Wolfgang Kapferer
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 all quantities used in the hydro computation are initiated with zeros
 */
int init_ppm(double **rho_1D, double **pre_1D, double **eng_1D, double **vx_1D, double **vy_1D,
        double **vz_1D, double **marker_1D, double **dx0, double **xa0, double **xa, double **dx, double **a_coef,
        double **ai_coef, double **b_coef, double **bi_coef, double **c_coef, double **ci_coef,
        double **d_x, double **da, double **ar, double **dp, double **dr, double **du,
        double **pl, double **p6, double **rl, double **r6, double **ul, double **u6,
        double **vl, double **v6, double **wl, double **w6, double **el, double **e6,
        double **ql, double **q6, double **deltaa, double **dv, double **dw, double **dq,
        double **de, double **scratch1, double **scratch2, double **scratch3, double **diffa,
        double **plft, double **prgh, double **ulft, double **urgh, double **rlft, double **rrgh,
        double **Cdtdx, double **fCdtdx, double **clft, double **crgh, double **plfti, double **prghi,
        double **pmid, double **pmold, double **wlft, double **wrgh, double **zlft, double **zrgh,
        double **umidl, double **umidr, double **umid, double **dm, double **dtbdm, double **upmid,
        double **xa1, double **xa2, double **xa3, double **vx_1D_old, double **e_int_1D, double **dvol,
        double **dvol0, double **dvol1, double **delta, double **fluxr, double **fluxu, double **fluxv,
        double **fluxw, double **fluxe, double **fluxq, double **dm0, double **steep, double **flat,
        double ***para, double **pressure_solid_1D, double **rhodown, double **rhoup, double **rhofront,
        double **rhoback, double **vxdown, double **vxup, double **vxfront, double **vxback, double **vydown,
        double **vyup, double **vyfront, double **vyback, double **vzdown, double **vzup, double **vzfront,
        double **vzback, int dimension) {

    int i;

    // the arrays which hold the hydro quantities
    *rho_1D = calloc((max_array_length + 12), sizeof ** rho_1D);
    *pre_1D = calloc((max_array_length + 12), sizeof ** pre_1D);
    *eng_1D = calloc((max_array_length + 12), sizeof ** eng_1D);
    *vx_1D = calloc((max_array_length + 12), sizeof ** vx_1D);
    *vy_1D = calloc((max_array_length + 12), sizeof ** vy_1D);
    *vz_1D = calloc((max_array_length + 12), sizeof ** vz_1D);

    *marker_1D = calloc((max_array_length + 12), sizeof ** marker_1D);
    *pressure_solid_1D = calloc((max_array_length + 12), sizeof ** pressure_solid_1D);

    //for viscosity calculations
    if (dimension > 1) {
        *rhodown = calloc((max_array_length + 12), sizeof ** rhodown);
        *rhoup = calloc((max_array_length + 12), sizeof ** rhoup);
        *vxdown = calloc((max_array_length + 12), sizeof ** vxdown);
        *vxup = calloc((max_array_length + 12), sizeof ** vxup);
        *vydown = calloc((max_array_length + 12), sizeof ** vydown);
        *vyup = calloc((max_array_length + 12), sizeof ** vyup);
        *vzdown = calloc((max_array_length + 12), sizeof ** vzdown);
        *vzup = calloc((max_array_length + 12), sizeof ** vzup);
    }
    if (dimension > 2) {
        *rhofront = calloc((max_array_length + 12), sizeof ** rhofront);
        *rhoback = calloc((max_array_length + 12), sizeof ** rhoback);
        *vxfront = calloc((max_array_length + 12), sizeof ** vxfront);
        *vxback = calloc((max_array_length + 12), sizeof ** vxback);
        *vyfront = calloc((max_array_length + 12), sizeof ** vyfront);
        *vyback = calloc((max_array_length + 12), sizeof ** vyback);
        *vzfront = calloc((max_array_length + 12), sizeof ** vzfront);
        *vzback = calloc((max_array_length + 12), sizeof ** vzback);
    }

    *dx0 = calloc((max_array_length + 12), sizeof ** dx0);
    *xa0 = calloc((max_array_length + 12), sizeof ** xa0);
    *xa = calloc((max_array_length + 12), sizeof ** xa);
    *dx = calloc((max_array_length + 12), sizeof ** dx);
    //allocates for the ppm routine
    //for the parabolic interpolation constants
    *a_coef = calloc((max_array_length + 12), sizeof ** a_coef);
    *ai_coef = calloc((max_array_length + 12), sizeof ** ai_coef);
    *b_coef = calloc((max_array_length + 12), sizeof ** b_coef);
    *bi_coef = calloc((max_array_length + 12), sizeof ** bi_coef);
    *c_coef = calloc((max_array_length + 12), sizeof ** c_coef);
    *ci_coef = calloc((max_array_length + 12), sizeof ** ci_coef);
    *d_x = calloc((max_array_length + 12), sizeof ** d_x);
    //for the parabolic interpolation
    *da = calloc((max_array_length + 12), sizeof ** da);
    *ar = calloc((max_array_length + 12), sizeof ** ar);
    *dp = calloc((max_array_length + 12), sizeof ** dp);
    *dr = calloc((max_array_length + 12), sizeof ** dr);
    *du = calloc((max_array_length + 12), sizeof ** du);
    *pl = calloc((max_array_length + 12), sizeof ** pl);
    *p6 = calloc((max_array_length + 12), sizeof ** p6);
    *rl = calloc((max_array_length + 12), sizeof ** rl);
    *r6 = calloc((max_array_length + 12), sizeof ** r6);
    *ul = calloc((max_array_length + 12), sizeof ** ul);
    *u6 = calloc((max_array_length + 12), sizeof ** u6);
    *vl = calloc((max_array_length + 12), sizeof ** vl);
    *v6 = calloc((max_array_length + 12), sizeof ** v6);
    *wl = calloc((max_array_length + 12), sizeof ** wl);
    *w6 = calloc((max_array_length + 12), sizeof ** w6);
    *el = calloc((max_array_length + 12), sizeof ** el);
    *e6 = calloc((max_array_length + 12), sizeof ** e6);
    *ql = calloc((max_array_length + 12), sizeof ** ql);
    *q6 = calloc((max_array_length + 12), sizeof ** q6);
    *deltaa = calloc((max_array_length + 12), sizeof ** deltaa);
    *dv = calloc((max_array_length + 12), sizeof ** dv);
    *dw = calloc((max_array_length + 12), sizeof ** dw);
    *dq = calloc((max_array_length + 12), sizeof ** dq);
    *de = calloc((max_array_length + 12), sizeof ** de);
    *scratch1 = calloc((max_array_length + 12), sizeof ** scratch1);
    *scratch2 = calloc((max_array_length + 12), sizeof ** scratch2);
    *scratch3 = calloc((max_array_length + 12), sizeof ** scratch3);
    *diffa = calloc((max_array_length + 12), sizeof ** diffa);
    //for the states
    *plft = calloc((max_array_length + 12), sizeof ** plft);
    *prgh = calloc((max_array_length + 12), sizeof ** prgh);
    *ulft = calloc((max_array_length + 12), sizeof ** ulft);
    *urgh = calloc((max_array_length + 12), sizeof ** urgh);
    *rlft = calloc((max_array_length + 12), sizeof ** rlft);
    *rrgh = calloc((max_array_length + 12), sizeof ** rrgh);
    *Cdtdx = calloc((max_array_length + 12), sizeof ** Cdtdx);
    *fCdtdx = calloc((max_array_length + 12), sizeof ** fCdtdx);
    //for the riemann solver
    *clft = calloc((max_array_length + 12), sizeof ** clft);
    *crgh = calloc((max_array_length + 12), sizeof ** crgh);
    *plfti = calloc((max_array_length + 12), sizeof ** plfti);
    *prghi = calloc((max_array_length + 12), sizeof ** prghi);
    *pmid = calloc((max_array_length + 12), sizeof ** pmid);
    *pmold = calloc((max_array_length + 12), sizeof ** pmold);
    *wlft = calloc((max_array_length + 12), sizeof ** wlft);
    *wrgh = calloc((max_array_length + 12), sizeof ** wrgh);
    *zlft = calloc((max_array_length + 12), sizeof ** zlft);
    *zrgh = calloc((max_array_length + 12), sizeof ** zrgh);
    *umidl = calloc((max_array_length + 12), sizeof ** umidl);
    *umidr = calloc((max_array_length + 12), sizeof ** umidr);
    *umid = calloc((max_array_length + 12), sizeof ** umid);
    // for the evolution 
    *dm = calloc((max_array_length + 12), sizeof ** dm);
    *dtbdm = calloc((max_array_length + 12), sizeof ** dtbdm);
    *upmid = calloc((max_array_length + 12), sizeof ** upmid);
    *xa1 = calloc((max_array_length + 12), sizeof ** xa1);
    *xa2 = calloc((max_array_length + 12), sizeof ** xa2);
    *xa3 = calloc((max_array_length + 12), sizeof ** xa3);
    *vx_1D_old = calloc((max_array_length + 12), sizeof ** vx_1D_old);
    *e_int_1D = calloc((max_array_length + 12), sizeof ** e_int_1D);
    *dvol = calloc((max_array_length + 12), sizeof ** dvol);
    *dvol0 = calloc((max_array_length + 12), sizeof ** dvol0);
    *dvol1 = calloc((max_array_length + 12), sizeof ** dvol1);
    // for remapping 
    *delta = calloc((max_array_length + 12), sizeof ** delta);
    *fluxr = calloc((max_array_length + 12), sizeof ** fluxr);
    *fluxu = calloc((max_array_length + 12), sizeof ** fluxu);
    *fluxv = calloc((max_array_length + 12), sizeof ** fluxv);
    *fluxw = calloc((max_array_length + 12), sizeof ** fluxw);
    *fluxe = calloc((max_array_length + 12), sizeof ** fluxe);
    *fluxq = calloc((max_array_length + 12), sizeof ** fluxq);
    *dm0 = calloc((max_array_length + 12), sizeof ** dm0);
    // steep and flat coefficients for the flatten function
    *steep = calloc((max_array_length + 12), sizeof ** steep);
    *flat = calloc((max_array_length + 12), sizeof ** flat);
    //for the para_coefficient
    *para = calloc((max_array_length + 12), sizeof **para);
    for (i = 0; i < (max_array_length + 12); i++) {
        *((*para) + i) = calloc(5, sizeof ***para);
    }

    return 0;
}

/*!
 all quantities within the hydro part are freed in this function
 */
int ppm_free(double **rho_1D, double **pre_1D, double **eng_1D, double **vx_1D, double **vy_1D,
        double **vz_1D, double **marker_1D, double **pressure_solid_1D,
        double **dx0, double **xa0, double **xa, double **dx, double **a_coef,
        double **ai_coef, double **b_coef, double **bi_coef, double **c_coef, double **ci_coef,
        double **d_x, double **da, double **ar, double **dp, double **dr, double **du,
        double **pl, double **p6, double **rl, double **r6, double **ul, double **u6,
        double **vl, double **v6, double **wl, double **w6, double **el, double **e6,
        double **ql, double **q6, double **deltaa, double **dv, double **dw, double **dq,
        double **de, double **scratch1, double **scratch2, double **scratch3, double **diffa,
        double **plft, double **prgh, double **ulft, double **urgh, double **rlft, double **rrgh,
        double **Cdtdx, double **fCdtdx, double **clft, double **crgh, double **plfti, double **prghi,
        double **pmid, double **pmold, double **wlft, double **wrgh, double **zlft, double **zrgh,
        double **umidl, double **umidr, double **umid, double **dm, double **dtbdm, double **upmid,
        double **xa1, double **xa2, double **xa3, double **vx_1D_old, double **e_int_1D, double **dvol,
        double **dvol0, double **dvol1, double **delta, double **fluxr, double **fluxu, double **fluxv,
        double **fluxw, double **fluxe, double **fluxq, double **dm0, double **steep, double **flat,
        double ***para, double **rhodown, double **rhoup, double **rhofront,
        double **rhoback, double **vxdown, double **vxup, double **vxfront, double **vxback, double **vydown,
        double **vyup, double **vyfront, double **vyback, double **vzdown, double **vzup, double **vzfront,
        double **vzback, int dimension) {
    int i;

    // the arrays which hold the hydro quantities
    free(*rho_1D);
    free(*pre_1D);
    free(*eng_1D);
    free(*vx_1D);
    free(*vy_1D);
    free(*vz_1D);

    free(*marker_1D);
    free(*pressure_solid_1D);

    //free the 1D quantities for the viscosity calculations
    if (dimension > 1) {
        free(*rhodown);
        free(*rhoup);
        free(*vxdown);
        free(*vxup);
        free(*vydown);
        free(*vyup);
        free(*vzdown);
        free(*vzup);
    }
    if (dimension > 2) {
        free(*rhofront);
        free(*rhoback);
        free(*vxfront);
        free(*vxback);
        free(*vyfront);
        free(*vyback);
        free(*vzfront);
        free(*vzback);
    }

    free(*dx0);
    free(*xa0);
    free(*xa);
    free(*dx);
    //allocates for the ppm routine
    //for the parabolic interpolation constants
    free(*a_coef);
    free(*ai_coef);
    free(*b_coef);
    free(*bi_coef);
    free(*c_coef);
    free(*ci_coef);
    free(*d_x);
    free(*da);
    free(*ar);
    free(*dp);
    free(*dr);
    free(*du);
    free(*pl);
    free(*p6);
    free(*rl);
    free(*r6);
    free(*ul);
    free(*u6);
    free(*vl);
    free(*v6);
    free(*wl);
    free(*w6);
    free(*el);
    free(*e6);
    free(*ql);
    free(*q6);
    free(*deltaa);
    free(*dv);
    free(*dw);
    free(*dq);
    free(*de);
    free(*scratch1);
    free(*scratch2);
    free(*scratch3);
    free(*diffa);
    //for states
    free(*plft);
    free(*prgh);
    free(*ulft);
    free(*urgh);
    free(*rlft);
    free(*rrgh);
    free(*Cdtdx);
    free(*fCdtdx);
    //Riemann solver
    free(*clft);
    free(*crgh);
    free(*plfti);
    free(*prghi);
    free(*pmid);
    free(*pmold);
    free(*wlft);
    free(*wrgh);
    free(*zlft);
    free(*zrgh);
    free(*umidl);
    free(*umidr);
    free(*umid);
    //evolve
    free(*dm);
    free(*dtbdm);
    free(*upmid);
    free(*xa1);
    free(*xa2);
    free(*xa3);
    free(*vx_1D_old);
    free(*e_int_1D);
    free(*dvol);
    free(*dvol0);
    free(*dvol1);
    //remap
    free(*delta);
    free(*fluxr);
    free(*fluxu);
    free(*fluxv);
    free(*fluxw);
    free(*fluxe);
    free(*fluxq);
    free(*dm0);
    //steep and flat for flatten
    free(*steep);
    free(*flat);
    //for para_coef
    for (i = 0; i < (max_array_length + 12); i++) {
        free(*((*para) + i));
    }
    free(*para);

    return 0;
}

int init_diffusion(double **temperature_1D, double **temperature_1D_future, double **mass, int **dom_1D) {

    *temperature_1D = calloc((max_array_length + 4), sizeof **temperature_1D);
    *temperature_1D_future = calloc((max_array_length + 4), sizeof **temperature_1D_future);
    *mass = calloc((max_array_length + 4), sizeof **mass);
    *dom_1D = calloc((max_array_length + 4), sizeof **dom_1D);

    return 0;
}

int free_diffusion(double **temperature_1D, double **temperature_1D_future, double **mass, int **dom_1D) {

    free(*temperature_1D);
    free(*temperature_1D_future);
    free(*mass);
    free(*dom_1D);

    return 0;
}

int init_viscosity(double **viscosity_array, double **vx_flux, double **vy_flux, double **vz_flux, int **dom_1D) {

    *viscosity_array = calloc((max_array_length + 2), sizeof **viscosity_array);
    *vx_flux = calloc((max_array_length + 2), sizeof **vx_flux);
    *vy_flux = calloc((max_array_length + 2), sizeof **vy_flux);
    *vz_flux = calloc((max_array_length + 2), sizeof **vz_flux);
    *dom_1D = calloc((max_array_length + 2), sizeof **dom_1D);


    return 0;
}

int free_viscosity(double **viscosity_array, double **vx_flux, double **vy_flux, double **vz_flux, int **dom_1D) {

    free(*viscosity_array);
    free(*vx_flux);
    free(*vy_flux);
    free(*vz_flux);
    free(*dom_1D);

    return 0;
}