/*
 * prototypes.h
 * all the functions used in TYCHO
 *
 * Author: Wolfgang Kapferer
 */

#ifndef PROTOTYPES_H_
#define PROTOTYPES_H_

int callocate_arrays_global(int x, int y, int z, int max_array_length);
int callocate_arrays_sweeps(int max_array_length);
int callocate_arrays_ppm(int max_array_length);
int free_sweep(void);
int free_ppm(void);

int initiate_grid(int x, double xmin, double xmax, int y, double ymin, double ymax, int z, double zmin, double zmax);
int initiate_domain(int x, int y, int z);
int initiate_domain_dom(int x, int y, int z);
int initiate_domain_marker(int x, int y, int z);

int hydro_sweeps(int x, int y, int z, int direction);

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
        double **vzback, int dimension);

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
        double **vzback, int dimension);

int init_diffusion(double **temperature_1D, double **temperature_1D_future, double **mass, int **dom_1D);
int init_diffusion(double **temperature_1D, double **temperature_1D_future, double **mass, int **dom_1D);

int make_ic(int x, int y, int z);
int make_sod_ic(int x, int y, int z);
int make_kh_instabilities(int x, int y, int z);

int read_ic(int x, int y, int z);
int read_restart(int x, int y, int z);
int read_dom(int x, int y, int z);

int read_marker_file(int x, int y, int z);

int start_file_reader(char filename[]);

int calculate_pressure(int x, int y, int z, double temperature);

long long int ntohll(const long long int data);

int write_ic_tyc(int x, int y, int z);
int write_ic_vtk(int x, int y, int z);
int write_ic_amira(int x, int y, int z);
int write_ic_ifrit(int x, int y, int z);

int write_tyc(int x, int y, int z, int counter);
int write_vtk(int x, int y, int z, int counter);
int write_amira(int x, int y, int z, int counter);
int write_ifrit(int x, int y, int z, int counter);

int write_restart_tyc(int x, int y, int z);

int dt_calc(int x, int y, int z);
int dt_calc_first(int x, int y, int z);

int sweep_x(int x, int y, int z, int flag);
int sweep_y(int x, int y, int z, int flag);
int sweep_z(int x, int y, int z, int flag);

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
        double *dm, double *dtbdm, double *upmid, double *xa1,
        double *xa2, double *xa3, double *vx_1D_old, double *dvol,
        double *dvol0, double *dvol1, double *delta, double *fluxr,
        double *fluxu, double *fluxv, double *fluxw, double *fluxe,
        double *fluxq, double *dm0, double *e_int_1D, double *rho_1D,
        double *pre_1D, double *eng_1D, double *vx_1D, double *vy_1D,
        double *vz_1D, double *marker_1D, double *pressure_solid_1D,
        double *dx0, double *xa0, double *xa, double *dx, int bound_checker, int lefter,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown,
        double *vxup, double *vxfront, double *vxback, double *vydown, double *vyup,
        double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int set_boundary(int nmin, int nmax, int flag, int bound_checker, int lefter,
        double *rho_1D, double *eng_1D, double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, double *rhodown,
        double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback, double *vzdown,
        double *vzup, double *vzfront, double *vzback, int viscosity_on_off, int dimension);

int left_boundary_zero_gradient(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int right_boundary_zero_gradient(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int left_boundary_reflecting(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int right_boundary_reflecting(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int left_boundary_small_padding_on_obstacle(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback, double *vzdown,
        double *vzup, double *vzfront, double *vzback, int viscosity_on_off, int dimension);

int right_boundary_reflecting_on_obstacle(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback, double *vzdown,
        double *vzup, double *vzfront, double *vzback, int viscosity_on_off, int dimension);

int left_boundary_reflecting_on_obstacle(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *pressure_solid_1D, double *xa0, double *dx0, double *xa, double *dx, int flag,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback, double *vxdown, double *vxup, double *vxfront,
        double *vxback, double *vydown, double *vyup, double *vyfront, double *vyback, double *vzdown,
        double *vzup, double *vzfront, double *vzback, int viscosity_on_off, int dimension);

int left_boundary_small_padding(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int right_boundary_small_padding(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int left_boundary_outflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int right_boundary_outflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int left_boundary_inflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int right_boundary_inflow(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int left_boundary_periodic(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int right_boundary_periodic(int nmin, int nmax, double *rho_1D, double *eng_1D,
        double *pre_1D, double *vx_1D, double *vy_1D, double *vz_1D,
        double *xa0, double *dx0, double *xa, double *dx, int flag, double *rhodown, double *rhoup, double *rhofront,
        double *rhoback, double *vxdown, double *vxup, double *vxfront, double *vxback, double *vydown,
        double *vyup, double *vyfront, double *vyback, double *vzdown, double *vzup, double *vzfront,
        double *vzback, int viscosity_on_off, int dimension);

int para_coef(int nmin, int nmax, int flag1, double *a_coef, double *dx,
        double *ai_coef, double *b_coef, double *bi_coef, double *c_coef,
        double *d_x, double **para, double *ci_coef);

int volume(int nmin, int nmax, double *vol, double *dx, double *vol0, double *dx0);

int flatten(int nmin, int nmax, double *pre_1D, double *vx_1D,
        double *steep, double *flat);

int parabola(int nmin, int nmax, double *a, double *deltaa,
        double *a6, double *al, int flag1, double *diffa,
        double *da, double **para, double *ar, double *flat,
        double *scratch1, double *scratch2, double *scratch3);

int states(int nmin, int nmax, int flag, double *pre_1D, double *rho_1D,
        double *dx, double *Cdtdx, double *fCdtdx, double *plft,
        double *pl, double *dp, double *p6, double *ulft, double *ul,
        double *du, double *u6, double *rlft, double *rl, double *dr,
        double *r6, double *prgh, double *urgh, double *rrgh);

int riemann(int nmin, int nmax, double *clft, double *crgh,
        double *rlft, double *rrgh, double *plfti,
        double *prghi, double *pmid, double *pmold,
        double *plft, double *wrgh, double *prgh, double *wlft,
        double *zlft, double *zrgh, double *umidl, double *umidr,
        double *umid, double *urgh, double *ulft);

int evolve(int nmin, int nmax, int flag, double *rho_1D,
        double *dvol, double *dm, double *dtbdm, double *xa1,
        double *xa, double *dvol1, double *umid, double *upmid,
        double *pmid, double *xa2, double *dx, double *xa3,
        double *vx_1D_old, double *vx_1D, double *vy_1D,
        double *vz_1D, double *eng_1D, double *e_int_1D,
        double *pre_1D);

int remap(int nmin, int nmax, int flag, double *a_coef, double *dx,
        double *ai_coef, double *b_coef, double *bi_coef, double *c_coef,
        double *d_x, double **para, double *ci_coef, double *dr, double *r6, double *rl, double *diffa,
        double *da, double *ar, double *flat, double *scratch1, double *scratch2, double *scratch3,
        double *du, double *u6, double *ul, double *dv, double *v6, double *vl, double *w6, double *wl,
        double *dq, double *q6, double *ql, double *de, double *e6, double *el, double *xa, double *xa0, double *delta,
        double *fluxr, double *fluxu, double *fluxv, double *fluxw, double *dw, double *fluxe, double *fluxq,
        double *dm, double *rho_1D, double *dvol, double *dvol0, double *dm0, double * vx_1D, double *vy_1D,
        double *vz_1D, double *eng_1D, double *e_int_1D, double *pre_1D);

int advect(int nmin, int nmax, int flag, double *dx, double *vx_1D, double *marker_1D);

int diffusion(int x, int y, int z);

double alpha_heat_transfer(int x, int y, int z, int direction);

int copy_arrays_for_viscosity(int x, int y, int z);

int viscosity(int i, int j, int k, int flag, int nmin, int nmax,
        double *rho_1D, double * vx_1D, double *vy_1D, double *vz_1D,
        double *pre_1D, double *e_int_1D, double *eng_1D, int lefter,
        double *rhodown, double *rhoup, double *rhofront, double *rhoback,
        double *vxdown, double *vxup, double *vxfront, double *vxback,
        double *vydown, double *vyup, double *vyfront, double *vyback,
        double *vzdown, double *vzup, double *vzfront, double *vzback,
        int dimension);

double kinematic_viscosity(double temperature, double density);

double dynamic_viscosity(double temperature);

int wind(int x, int y, int z);

int upper_atmos(int x, int y, int z);

//For data exploration
int velocity_field_analyser(int x, int y, int z);

int pressure_on_solid_calc(int x, int y, int z);

//The wind tunnel corrector
int boundary_velo_corrector(int x, int y, int z);

//An explosion
int insert_pressure(int x, int y, int z);

//to reset the dB-Map
int reset_pressure_integrated(int x, int y, int z);

//to prepare the dB-Map
int prepare_the_dB_map(int x, int y, int z);
int pre_old_copy(int x, int y, int z);

#endif /* PROTOTYPES_H_ */
