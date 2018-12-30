/*
 * dt_calc.c
 *
 * Author: Wolfgang Kapferer
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
inline double max_dt(double a, double b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;
}

/*!
 Inline function to calculate the minimum of two doubles
 */
inline double min_dt(double a, double b) {
    double tmp;

    if (a < b) tmp = a;
    if (a == b)tmp = a;
    if (a > b) tmp = b;

    return tmp;
}

/*!
 This function calculates the hydrodynamical timestep 
 for the simulation. It is based on the signal 
 velocity in the medium and the courant factor,
 specified in the parameterfile.
 */
int dt_calc(int x, int y, int z) {
    int i, j, k;
    double rdt1, tmpdy, tmpdz, tmp;
    double xvel, yvel, zvel, svel;
    double max_tmp1, max_tmp2, max_tmp3, min_tmp1, dtx, dt3;
    double kin_viscosity;
    double s_visc;
    double temperature;

    rdt1 = 0.0;
    svel = 0.0;
    kin_viscosity = 0.0;
    s_visc = 0.0;

    /*
     The signal velocity is calculated and compared to the
     gas velocity in each cell. The routine searches for
     the maximum in the whole computational domain.
     */
    if (dimension == 1) {
        k = 0;
        j = 0;

        for (i = 0; i < x; i++) {
            if (dom[i][j][k] == 0) {
                svel = sqrt(Gamma1 * pre[i][j][k] / rho[i][j][k]) / zdx[i];
                xvel = fabs(vx[i][j][k]) / zdx[i];
                tmp = max_dt(svel, xvel);
                rdt1 = max_dt(tmp, rdt1);

                if (viscosity_on_off == 1) {
                    temperature = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                    kin_viscosity = kinematic_viscosity(temperature, rho[i][j][k]);
                    s_visc = kin_viscosity / (tmp * tmp);
                    rdt1 = max_dt(s_visc, rdt1);
                }
            }
        }
    }


#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, tmpdy, tmp, svel, xvel, \
                        yvel, max_tmp2, max_tmp3, temperature, \
                        kin_viscosity,s_visc) \
                shared(x, y, z, dom, rho, pre, Gamma1, vx, \
                       vy, zdx, zdy, dimension, rdt1, \
                       viscosity_on_off, gasconstant)
    {
#endif
        if (dimension == 2) {
            k = 0;

            for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif    
                for (j = 0; j < y; j++) {
                    if (dom[i][j][k] == 0) {
                        tmpdy = zdy[j];
                        tmp = min_dt(zdx[i], tmpdy);
                        svel = sqrt(Gamma1 * pre[i][j][k] / rho[i][j][k]) / tmp;
                        xvel = fabs(vx[i][j][k]) / zdx[i];
                        yvel = fabs(vy[i][j][k]) / tmpdy;
                        max_tmp2 = max_dt(svel, xvel);
                        max_tmp3 = max_dt(max_tmp2, yvel);
                        rdt1 = max_dt(max_tmp3, rdt1);

                        if (viscosity_on_off == 1) {
                            temperature = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                            kin_viscosity = kinematic_viscosity(temperature, rho[i][j][k]);
                            s_visc = kin_viscosity / (tmp * tmp);
                            rdt1 = max_dt(s_visc, rdt1);
                        }
                    }
                }
            }

        }

#ifdef _OPENMP
    }
#endif


#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, tmpdy, tmp, tmpdz, svel, xvel, \
                        yvel, zvel, max_tmp1, max_tmp2, max_tmp3, \
                        min_tmp1, temperature, \
                        kin_viscosity,s_visc) \
                shared(x, y, z, dom, rho, pre, Gamma1, vx, \
                       vy, vz, zdx, zdy, zdz,dimension, rdt1, \
                       viscosity_on_off, gasconstant)
    {
#endif

        if (dimension == 3) {

            for (i = 0; i < x; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif    
                for (j = 0; j < y; j++) {
                    for (k = 0; k < z; k++) {
                        if (dom[i][j][k] == 0) {
                            tmpdy = zdy[j];
                            tmpdz = zdz[k];
                            min_tmp1 = min_dt(tmpdy, tmpdz);
                            tmp = min_dt(zdx[i], min_tmp1);
                            svel = sqrt(Gamma1 * pre[i][j][k] / rho[i][j][k]) / tmp;
                            xvel = fabs(vx[i][j][k]) / zdx[i];
                            yvel = fabs(vy[i][j][k]) / tmpdy;
                            zvel = fabs(vz[i][j][k]) / tmpdz;
                            max_tmp1 = max_dt(xvel, yvel);
                            max_tmp2 = max_dt(zvel, svel);
                            max_tmp3 = max_dt(max_tmp1, max_tmp2);
                            rdt1 = max_dt(max_tmp3, rdt1);

                            if (viscosity_on_off == 1) {
                                temperature = pre[i][j][k] / (gasconstant * rho[i][j][k]);
                                kin_viscosity = kinematic_viscosity(temperature, rho[i][j][k]);
                                s_visc = kin_viscosity / (tmp * tmp);
                                rdt1 = max_dt(s_visc, rdt1);
                            }
                        }
                    }
                }
            }
        }

#ifdef _OPENMP
    }
#endif

    // here the CFL condition is applied
    dtx = courant / rdt1;

    // In order to avoid strong time step changes
    dt3 = 1.2 * olddt;
    dt = min_dt(dt3, dtx);

    if (dt == 0.0) {
        printf("A time_simstep of zero ---> STOP\n");
        exit(42);
    }
    
    return 0;
}
