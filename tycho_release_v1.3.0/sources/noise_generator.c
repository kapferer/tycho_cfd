/*
 * noise_generator.c
 *
 * Author: Wolfgang Kapferer
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 A simple noise generator, i.e. periodic pressure changes in
 soundemitter == 1 areas
 */
int noise_generator(int x, int y, int z) {

    double local_pressure_left, local_pressure_right;
    int one_boundary;
    int left, right;
    int counter_mean_pressure;
    int i, j, k;
    double tmp;
    double perfect_sinus_emitter;
    
    perfect_sinus_emitter = 1.0 / sqrt(2);

    mean_pressure = 0.0;
    counter_mean_pressure = 0;
    //==========================================================Here the mean pressure in the domain===============================================

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, local_pressure_left, \
                local_pressure_right, one_boundary, \
                left, right) \
                shared(x, y, z, pre, dom, \
                pressure_integrated, mean_pressure, counter_mean_pressure)
    {
#endif      

        for (i = 0; i < x; i++) {
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {

                    if (dom[i][j][k] == 0) {

#ifdef _OPENMP
#pragma omp critical
#endif
                        mean_pressure += pre[i][j][k];
#ifdef _OPENMP
#pragma omp end critical
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
                        counter_mean_pressure++;
#ifdef _OPENMP
#pragma omp end critical
#endif

                    }
                }
            }
        }

#ifdef _OPENMP
    }
#endif        

    mean_pressure = mean_pressure / counter_mean_pressure;
    if (with_one_pulse == 0) pressure = mean_pressure + sound_pressure_level * sin(2 * M_PI * sound_frequency * time_sim);

    //just one pulse    
    if ((with_one_pulse == 1) && (counter == 0)) {
        pressure = fabs(mean_pressure + sound_pressure_level * sin(0.5 * M_PI));
    }

    if (dimension == 3) {

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, pre_old, pre, \
                soundemitter, dom, pressure, \
                pressure_on_solid, with_one_pulse, counter)
        {
#endif
            for (i = 1; i < x - 1; i++) {

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif

                for (j = 1; j < y - 1; j++) {
                    for (k = 1; k < z - 1; k++) {

                        if ((with_one_pulse == 1) && (counter == 0)) {
                            if (soundemitter[i - 1][j - 1][k - 1] != 1) pre[i - 1][j - 1][k - 1] = pressure;
                            if (soundemitter[i - 1][j ][k - 1] != 1) pre[i - 1][j ][k - 1] = pressure;
                            if (soundemitter[i - 1][j + 1][k - 1] != 1) pre[i - 1][j + 1][k - 1] = pressure;
                            if (soundemitter[i - 1][j - 1][k] != 1) pre[i - 1][j - 1][k] = pressure;
                            if (soundemitter[i - 1][j][k] != 1) pre[i - 1][j][k] = pressure;
                            if (soundemitter[i - 1][j + 1][k] != 1) pre[i - 1][j + 1][k] = pressure;
                            if (soundemitter[i - 1][j - 1][k + 1] != 1) pre[i - 1][j - 1][k + 1] = pressure;
                            if (soundemitter[i - 1][j ][k + 1] != 1) pre[i - 1][j][k + 1] = pressure;
                            if (soundemitter[i - 1][j + 1][k + 1] != 1) pre[i - 1][j + 1][k + 1] = pressure;

                            if (soundemitter[i][j - 1][k - 1] != 1) pre[i][j - 1][k - 1] = pressure;
                            if (soundemitter[i][j][k - 1] != 1) pre[i][j][k - 1] = pressure;
                            if (soundemitter[i][j + 1][k - 1] != 1) pre[i][j + 1][k - 1] = pressure;
                            if (soundemitter[i][j - 1][k] != 1) pre[i][j - 1][k] = pressure;
                            if (soundemitter[i][j + 1][k] != 1) pre[i][j + 1][k - 1] = pressure;
                            if (soundemitter[i][j - 1][k + 1] != 1) pre[i][j - 1][k + 1] = pressure;
                            if (soundemitter[i][j][k + 1] != 1) pre[i][j][k + 1] = pressure;
                            if (soundemitter[i][j + 1][k + 1] != 1) pre[i][j + 1][k + 1] = pressure;

                            if (soundemitter[i + 1][j - 1][k - 1] != 1) pre[i + 1][j - 1][k - 1] = pressure;
                            if (soundemitter[i + 1][j][k - 1] != 1) pre[i + 1][j][k - 1] = pressure;
                            if (soundemitter[i + 1][j + 1][k - 1] != 1) pre[i + 1][j + 1][k - 1] = pressure;
                            if (soundemitter[i + 1][j - 1][k] != 1) pre[i + 1][j - 1][k] = pressure;
                            if (soundemitter[i + 1][j][k] != 1) pre[i + 1][j][k] = pressure;
                            if (soundemitter[i + 1][j + 1][k] != 1) pre[i + 1][j + 1][k] = pressure;
                            if (soundemitter[i + 1][j - 1][k + 1] != 1) pre[i + 1][j - 1][k + 1] = pressure;
                            if (soundemitter[i + 1][j][k + 1] != 1) pre[i + 1][j][k + 1] = pressure;
                            if (soundemitter[i + 1][j + 1][k + 1] != 1) pre[i + 1][j + 1][k + 1] = pressure;
                        }
                        if (with_one_pulse == 0) {
                            if (soundemitter[i - 1][j - 1][k - 1] != 1) pre[i - 1][j - 1][k - 1] = pressure;
                            if (soundemitter[i - 1][j ][k - 1] != 1) pre[i - 1][j ][k - 1] = pressure;
                            if (soundemitter[i - 1][j + 1][k - 1] != 1) pre[i - 1][j + 1][k - 1] = pressure;
                            if (soundemitter[i - 1][j - 1][k] != 1) pre[i - 1][j - 1][k] = pressure;
                            if (soundemitter[i - 1][j][k] != 1) pre[i - 1][j][k] = pressure;
                            if (soundemitter[i - 1][j + 1][k] != 1) pre[i - 1][j + 1][k] = pressure;
                            if (soundemitter[i - 1][j - 1][k + 1] != 1) pre[i - 1][j - 1][k + 1] = pressure;
                            if (soundemitter[i - 1][j ][k + 1] != 1) pre[i - 1][j][k + 1] = pressure;
                            if (soundemitter[i - 1][j + 1][k + 1] != 1) pre[i - 1][j + 1][k + 1] = pressure;

                            if (soundemitter[i][j - 1][k - 1] != 1) pre[i][j - 1][k - 1] = pressure;
                            if (soundemitter[i][j][k - 1] != 1) pre[i][j][k - 1] = pressure;
                            if (soundemitter[i][j + 1][k - 1] != 1) pre[i][j + 1][k - 1] = pressure;
                            if (soundemitter[i][j - 1][k] != 1) pre[i][j - 1][k] = pressure;
                            if (soundemitter[i][j + 1][k] != 1) pre[i][j + 1][k - 1] = pressure;
                            if (soundemitter[i][j - 1][k + 1] != 1) pre[i][j - 1][k + 1] = pressure;
                            if (soundemitter[i][j][k + 1] != 1) pre[i][j][k + 1] = pressure;
                            if (soundemitter[i][j + 1][k + 1] != 1) pre[i][j + 1][k + 1] = pressure;

                            if (soundemitter[i + 1][j - 1][k - 1] != 1) pre[i + 1][j - 1][k - 1] = pressure;
                            if (soundemitter[i + 1][j][k - 1] != 1) pre[i + 1][j][k - 1] = pressure;
                            if (soundemitter[i + 1][j + 1][k - 1] != 1) pre[i + 1][j + 1][k - 1] = pressure;
                            if (soundemitter[i + 1][j - 1][k] != 1) pre[i + 1][j - 1][k] = pressure;
                            if (soundemitter[i + 1][j][k] != 1) pre[i + 1][j][k] = pressure;
                            if (soundemitter[i + 1][j + 1][k] != 1) pre[i + 1][j + 1][k] = pressure;
                            if (soundemitter[i + 1][j - 1][k + 1] != 1) pre[i + 1][j - 1][k + 1] = pressure;
                            if (soundemitter[i + 1][j][k + 1] != 1) pre[i + 1][j][k + 1] = pressure;
                            if (soundemitter[i + 1][j + 1][k + 1] != 1) pre[i + 1][j + 1][k + 1] = pressure;
                        }

                        if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                            if (dom[i - 1][j - 1][k - 1] != 1) pressure_on_solid[i - 1][j - 1][k - 1] = fabs(pre[i - 1][j - 1][k - 1] - pre_old[i - 1][j - 1][k - 1]);
                            if (dom[i - 1][j ][k - 1] != 1) pressure_on_solid[i - 1][j ][k - 1] = fabs(pre[i - 1][j ][k - 1] - pre_old[i - 1][j ][k - 1]);
                            if (dom[i - 1][j + 1][k - 1] != 1) pressure_on_solid[i - 1][j + 1][k - 1] = fabs(pre[i - 1][j + 1][k - 1] - pre_old[i - 1][j + 1][k - 1]);
                            if (dom[i - 1][j - 1][k] != 1) pressure_on_solid[i - 1][j - 1][k] = fabs(pre[i - 1][j - 1][k] - pre_old[i - 1][j - 1][k]);
                            if (dom[i - 1][j][k] != 1) pressure_on_solid[i - 1][j][k] = fabs(pre[i - 1][j][k] - pre_old[i - 1][j][k]);
                            if (dom[i - 1][j + 1][k] != 1) pressure_on_solid[i - 1][j + 1][k] = fabs(pre[i - 1][j + 1][k] - pre_old[i - 1][j + 1][k]);
                            if (dom[i - 1][j - 1][k + 1] != 1) pressure_on_solid[i - 1][j - 1][k + 1] = fabs(pre[i - 1][j - 1][k + 1] - pre_old[i - 1][j - 1][k + 1]);
                            if (dom[i - 1][j ][k + 1] != 1) pressure_on_solid[i - 1][j][k + 1] = fabs(pre[i - 1][j][k + 1] - pre_old[i - 1][j][k + 1]);
                            if (dom[i - 1][j + 1][k + 1] != 1) pressure_on_solid[i - 1][j + 1][k + 1] = fabs(pre[i - 1][j + 1][k + 1] - pre_old[i - 1][j + 1][k + 1]);

                            if (dom[i][j - 1][k - 1] != 1) pressure_on_solid[i][j - 1][k - 1] = fabs(pre[i][j - 1][k - 1] - pre_old[i][j - 1][k - 1]);
                            if (dom[i][j][k - 1] != 1) pressure_on_solid[i][j][k - 1] = fabs(pre[i][j][k - 1] - pre_old[i][j][k - 1]);
                            if (dom[i][j + 1][k - 1] != 1) pressure_on_solid[i][j + 1][k - 1] = fabs(pre[i][j + 1][k - 1] - pre_old[i][j + 1][k - 1]);
                            if (dom[i][j - 1][k] != 1) pressure_on_solid[i][j - 1][k] = fabs(pre[i][j - 1][k] - pre_old[i][j - 1][k]);
                            if (dom[i][j + 1][k] != 1) pressure_on_solid[i][j + 1][k - 1] = fabs(pre[i][j + 1][k - 1] - pre_old[i][j + 1][k - 1]);
                            if (dom[i][j - 1][k + 1] != 1) pressure_on_solid[i][j - 1][k + 1] = fabs(pre[i][j - 1][k + 1] - pre_old[i][j - 1][k + 1]);
                            if (dom[i][j][k + 1] != 1) pressure_on_solid[i][j][k + 1] = fabs(pre[i][j][k + 1] - pre_old[i][j][k + 1]);
                            if (dom[i][j + 1][k + 1] != 1) pressure_on_solid[i][j + 1][k + 1] = fabs(pre[i][j + 1][k + 1] - pre_old[i][j + 1][k + 1]);

                            if (dom[i + 1][j - 1][k - 1] != 1) pressure_on_solid[i + 1][j - 1][k - 1] = fabs(pre[i + 1][j - 1][k - 1] - pre_old[i + 1][j - 1][k - 1]);
                            if (dom[i + 1][j][k - 1] != 1) pressure_on_solid[i + 1][j][k - 1] = fabs(pre[i + 1][j][k - 1] - pre_old[i + 1][j][k - 1]);
                            if (dom[i + 1][j + 1][k - 1] != 1) pressure_on_solid[i + 1][j + 1][k - 1] = fabs(pre[i + 1][j + 1][k - 1] - pre_old[i + 1][j + 1][k - 1]);
                            if (dom[i + 1][j - 1][k] != 1) pressure_on_solid[i + 1][j - 1][k] = fabs(pre[i + 1][j - 1][k] - pre_old[i + 1][j - 1][k]);
                            if (dom[i + 1][j][k] != 1) pressure_on_solid[i + 1][j][k] = fabs(pre[i + 1][j][k] - pre_old[i + 1][j][k]);
                            if (dom[i + 1][j + 1][k] != 1) pressure_on_solid[i + 1][j + 1][k] = fabs(pre[i + 1][j + 1][k] - pre_old[i + 1][j + 1][k]);
                            if (dom[i + 1][j - 1][k + 1] != 1) pressure_on_solid[i + 1][j - 1][k + 1] = fabs(pre[i + 1][j - 1][k + 1] - pre_old[i + 1][j - 1][k + 1]);
                            if (dom[i + 1][j][k + 1] != 1) pressure_on_solid[i + 1][j][k + 1] = fabs(pre[i + 1][j][k + 1] - pre_old[i + 1][j][k + 1]);
                            if (dom[i + 1][j + 1][k + 1] != 1) pressure_on_solid[i + 1][j + 1][k + 1] = fabs(pre[i + 1][j + 1][k + 1] - pre_old[i + 1][j + 1][k + 1]);
                        }
                    }
                }
            }

#ifdef _OPENMP
        }
#endif

        //==========================================================Here the acoustic energy exchange===========================================================

        local_pressure_right = 0.0;
        local_pressure_left = 0.0;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, local_pressure_left, \
                local_pressure_right, one_boundary, \
                left, right) \
                shared(x, y, z, pre_old, pre, \
                soundemitter, dom, pressure, \
                pressure_on_solid, obstalce_absorption_coefficient, \
                pressure_integrated, mean_pressure, counter, \
                sound_reflexion)
        {
#endif

            for (j = 1; j < y - 1; j++) {

                left = -1;
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
                for (i = 1; i < x - 1; i++) {
                    for (k = 1; k < z - 1; k++) {

                        //now we simulate the obstacle acoustic-energy exchange 
                        if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                            if (dom[i - 1][j][k] != 1) {
                                local_pressure_left = pressure_on_solid[i - 1][j][k];
                                left = i - 1;
                            }
                            if ((dom[i + 1][j][k] != 1) && (left > 0)) {
                                local_pressure_right = pressure_on_solid[i + 1][j][k];
                                right = i + 1;
                                if (fabs(local_pressure_left) > fabs(local_pressure_right)) {
                                    pre[right][j][k] += obstalce_absorption_coefficient*local_pressure_left;
                                    pre[left][j][k] = pre_old[left][j][k] + (pressure_on_solid[left][j][k] * sound_reflexion);
                                }
                                if (fabs(local_pressure_left) < fabs(local_pressure_right)) {
                                    pre[left][j][k] += obstalce_absorption_coefficient*local_pressure_right;
                                    pre[right][j][k] = pre_old[right][j][k] + (pressure_on_solid[right][j][k] * sound_reflexion);
                                }
                            }
                        }
                    }
                }
            }

            for (i = 1; i < x - 1; i++) {

                left = -1;
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
                for (j = 1; j < y - 1; j++) {
                    for (k = 1; k < z - 1; k++) {

                        //now we simulate the obstacle acoustic-energy exchange 
                        if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                            if (dom[i][j - 1][k] != 1) {
                                local_pressure_left = pressure_on_solid[i][j - 1][k];
                                left = j - 1;
                            }

                            if ((dom[i][j + 1][k] != 1) && (left > 0)) {
                                local_pressure_right = pressure_on_solid[i][j + 1][k];
                                right = j + 1;
                                if (fabs(local_pressure_left) > fabs(local_pressure_right)) {
                                    pre[i][right][k] += obstalce_absorption_coefficient*local_pressure_left;
                                    pre[i][left][k] = pre_old[i][left][k] + (pressure_on_solid[i][left][k] * sound_reflexion);

                                }
                                if (fabs(local_pressure_left) < fabs(local_pressure_right)) {
                                    pre[i][left][k] += obstalce_absorption_coefficient*local_pressure_right;
                                    pre[i][right][k] = pre_old[i][right][k] + (pressure_on_solid[i][right][k] * sound_reflexion);
                                }
                            }
                        }
                    }
                }
            }

            for (i = 1; i < x - 1; i++) {
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
                for (j = 1; j < y - 1; j++) {

                    left = -1;
                    for (k = 1; k < z - 1; k++) {

                        //now we simulate the obstacle acoustic-energy exchange 
                        if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                            if (dom[i][j][k - 1] != 1) {
                                local_pressure_left = pressure_on_solid[i][j][k - 1];
                                left = k - 1;
                            }
                            if ((dom[i][j][k + 1] != 1) && (left > 0)) {
                                local_pressure_right = pressure_on_solid[i][j][k + 1];
                                right = k + 1;
                                if (fabs(local_pressure_left) > fabs(local_pressure_right)) {
                                    pre[i][j][right] += obstalce_absorption_coefficient*local_pressure_left;
                                    pre[i][j][left] = pre_old[i][j][left] + (pressure_on_solid[i][j][left] * sound_reflexion);
                                }
                                if (fabs(local_pressure_left) < fabs(local_pressure_right)) {
                                    pre[i][j][left] += obstalce_absorption_coefficient*local_pressure_right;
                                    pre[i][j][right] = pre_old[i][j][right] + (pressure_on_solid[i][j][right] * sound_reflexion);
                                }
                            }
                        }
                    }
                }
            }
#ifdef _OPENMP
        }
#endif 

        //==========================================================Here the acoustic energy exchange END========================================================

    }

    if (dimension == 2) {

        k = 0;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, pre_old, pre, \
                soundemitter, dom, pressure, \
                pressure_on_solid, mean_pressure, \
                with_one_pulse, counter)
        {
#endif

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif          
            for (i = 1; i < x - 1; i++) {
                for (j = 1; j < y - 1; j++) {

                    if (soundemitter[i][j][k] == 1) {

                        if ((with_one_pulse == 1) && (counter == 0)) {
                            if (soundemitter[i - 1][j - 1][k] != 1) pre[i - 1][j - 1][k] = pressure;
                            if (soundemitter[i - 1][j][k] != 1) pre[i - 1][j][k] = pressure;
                            if (soundemitter[i - 1][j + 1][k] != 1) pre[i - 1][j + 1][k] = pressure;

                            if (soundemitter[i][j - 1][k] != 1) pre[i][j - 1][k] = pressure;
                            if (soundemitter[i][j + 1][k] != 1) pre[i][j + 1][k] = pressure;


                            if (soundemitter[i + 1][j - 1][k] != 1) pre[i + 1][j - 1][k] = pressure;
                            if (soundemitter[i + 1][j][k] != 1) pre[i + 1][j][k] = pressure;
                            if (soundemitter[i + 1][j + 1][k] != 1) pre[i + 1][j + 1][k] = pressure;
                        }
                        if (with_one_pulse == 0) {
                            if (soundemitter[i - 1][j - 1][k] != 1) pre[i - 1][j - 1][k] = pressure;
                            if (soundemitter[i - 1][j][k] != 1) pre[i - 1][j][k] = pressure;
                            if (soundemitter[i - 1][j + 1][k] != 1) pre[i - 1][j + 1][k] = pressure;

                            if (soundemitter[i][j - 1][k] != 1) pre[i][j - 1][k] = pressure;
                            if (soundemitter[i][j + 1][k] != 1) pre[i][j + 1][k] = pressure;


                            if (soundemitter[i + 1][j - 1][k] != 1) pre[i + 1][j - 1][k] = pressure;
                            if (soundemitter[i + 1][j][k] != 1) pre[i + 1][j][k] = pressure;
                            if (soundemitter[i + 1][j + 1][k] != 1) pre[i + 1][j + 1][k] = pressure;
                        }
                    }

                    if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {

                        if (dom[i - 1][j - 1][k] != 1) pressure_on_solid[i - 1][j - 1][k] = pre[i - 1][j - 1][k] - pre_old[i - 1][j - 1][k];
                        if (dom[i - 1][j][k] != 1) pressure_on_solid[i - 1][j][k] = pre[i - 1][j][k] - pre_old[i - 1][j][k];
                        if (dom[i - 1][j + 1][k] != 1) pressure_on_solid[i - 1][j + 1][k] = pre[i - 1][j + 1][k] - pre_old[i - 1][j + 1][k];

                        if (dom[i][j - 1][k] != 1) pressure_on_solid[i][j - 1][k] = pre[i][j - 1][k] - pre_old[i][j - 1][k];
                        if (dom[i][j + 1][k] != 1) pressure_on_solid[i][j + 1][k] = pre[i][j + 1][k] - pre_old[i][j + 1][k];


                        if (dom[i + 1][j - 1][k] != 1) pressure_on_solid[i + 1][j - 1][k] = pre[i + 1][j - 1][k] - pre_old[i + 1][j - 1][k];
                        if (dom[i + 1][j][k] != 1) pressure_on_solid[i + 1][j][k] = pre[i + 1][j][k] - pre_old[i + 1][j][k];
                        if (dom[i + 1][j + 1][k] != 1) pressure_on_solid[i + 1][j + 1][k] = pre[i + 1][j + 1][k] - pre_old[i + 1][j + 1][k];
                    }

                }
            }
#ifdef _OPENMP
        }
#endif



        //==========================================================Here the acoustic energy exchange===========================================================
        //==========================================================Here the acoustic energy exchange===========================================================


        local_pressure_right = 0.0;
        local_pressure_left = 0.0;
        k = 0;


#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, local_pressure_left, \
                local_pressure_right, \
                left, right) \
                shared(x, y, z, pre_old, pre, \
                soundemitter, dom, pressure_on_solid, \
                obstalce_absorption_coefficient, \
                sound_reflexion, counter)
        {
#endif

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif

            for (j = 1; j < y - 1; j++) {
                left = -1;

                for (i = 1; i < x - 1; i++) {
                    //now we simulate the obstacle acoustic-energy exchange 
                    if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                        if (dom[i - 1][j][k] != 1) {
                            local_pressure_left = pressure_on_solid[i - 1][j][k];
                            left = i - 1;
                        }

                        if ((dom[i + 1][j][k] != 1) && (left > 0)) {
                            local_pressure_right = pressure_on_solid[i + 1][j][k];
                            right = i + 1;

                            //if (counter > 5) printf("left: %g   rhight: %g\n", local_pressure_left, local_pressure_right);

                            if (fabs(local_pressure_left) > fabs(local_pressure_right)) {
                                pre[right][j][k] += obstalce_absorption_coefficient*local_pressure_left;
                                pre[left][j][k] = pre_old[left][j][k] + (pressure_on_solid[left][j][k] * sound_reflexion);
                            }
                            if (fabs(local_pressure_left) < fabs(local_pressure_right)) {
                                pre[left][j][k] += obstalce_absorption_coefficient*local_pressure_right;
                                pre[right][j][k] = pre_old[right][j][k] + (pressure_on_solid[right][j][k] * sound_reflexion);
                            }
                        }
                    }
                }
            }
#ifdef _OPENMP
        }
#endif 

        local_pressure_right = 0.0;
        local_pressure_left = 0.0;
        k = 0;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, local_pressure_left, \
                local_pressure_right, \
                left, right) \
                shared(x, y, z, pre_old, pre, \
                soundemitter, dom, pressure_on_solid, \
                obstalce_absorption_coefficient, \
                sound_reflexion)
        {
#endif

#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
            for (i = 1; i < x - 1; i++) {
                left = -1;
                for (j = 1; j < y - 1; j++) {
                    //now we simulate the obstacle acoustic-energy exchange 
                    if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                        if (dom[i][j - 1][k] != 1) {
                            local_pressure_left = pressure_on_solid[i][j - 1][k];
                            left = j - 1;
                        }
                        if ((dom[i][j + 1][k] != 1) && (left > 0)) {
                            local_pressure_right = pressure_on_solid[i][j + 1][k];
                            right = j + 1;
                            if (fabs(local_pressure_left) > fabs(local_pressure_right)) {
                                pre[i][right][k] += obstalce_absorption_coefficient*local_pressure_left;
                                pre[i][left][k] = pre_old[i][left][k] + (pressure_on_solid[i][left][k] * sound_reflexion);
                            }
                            if (fabs(local_pressure_left) < fabs(local_pressure_right)) {
                                pre[i][left][k] += obstalce_absorption_coefficient*local_pressure_right;
                                pre[i][right][k] = pre_old[i][right][k] + (pressure_on_solid[i][right][k] * sound_reflexion);
                            }
                        }
                    }
                }
            }
#ifdef _OPENMP
        }
#endif


        //==========================================================Here the acoustic energy exchange END========================================================


    }


    local_pressure_right = 0.0;
    local_pressure_left = 0.0;
    k = 0;
    j = 0;

    if (dimension == 1) {
        for (i = 1; i < x - 1; i++) {
            if ((with_one_pulse == 1) && (counter == 0)) {
                if (soundemitter[i - 1][j][k] != 1) pre[i - 1][j][k] = pressure;
                if (soundemitter[i + 1][j][k] != 1) pre[i + 1][j][k] = pressure;
            }
            if (with_one_pulse == 0) {
                if (soundemitter[i - 1][j][k] != 1) pre[i - 1][j][k] = pressure;
                if (soundemitter[i + 1][j][k] != 1) pre[i + 1][j][k] = pressure;
            }
            if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                if (dom[i - 1][j][k] != 1) pressure_on_solid[i - 1][j][k] = fabs(pre[i - 1][j][k] - pre_old[i - 1][j][k]);
                if (dom[i + 1][j][k] != 1) pressure_on_solid[i + 1][j][k] = fabs(pre[i + 1][j][k] - pre_old[i + 1][j][k]);
            }
        }

        for (i = 1; i < x - 1; i++) {
            left = -1;
            //now we simulate the obstacle acoustic-energy exchange 
            if ((dom[i][j][k] == 1) && (soundemitter[i][j][k] == 0)) {
                if (dom[i - 1][j][k] != 1) {
                    local_pressure_left = pressure_on_solid[i - 1][j][k];
                    left = i - 1;
                }
                if ((dom[i + 1][j][k] != 1) && (left > 0)) {
                    local_pressure_right = pressure_on_solid[i + 1][j][k];
                    right = i + 1;
                    if (fabs(local_pressure_left) > fabs(local_pressure_right)) {
                        pre[right][j][k] += obstalce_absorption_coefficient*local_pressure_left;
                        pre[left][j][k] = pre_old[left][j][k] + (pressure_on_solid[left][j][k] * sound_reflexion);
                    }
                    if (fabs(local_pressure_left) < fabs(local_pressure_right)) {
                        pre[left][j][k] += obstalce_absorption_coefficient*local_pressure_right;
                        pre[right][j][k] = pre_old[right][j][k] + (pressure_on_solid[right][j][k] * sound_reflexion);
                    }
                }
            }
        }
    }


    //==========================================================Here parts of the dB Map are calculated============================================================

    //here the preparation for the periodic dB - Map is done
    if (with_one_pulse == 0) {

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, local_pressure_left, \
                local_pressure_right, one_boundary, \
                left, right) \
                shared(x, y, z, pre, dom, \
                pressure_integrated, mean_pressure, dt)
        {
#endif      

            for (i = 0; i < x; i++) {
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
                for (j = 0; j < y; j++) {
                    for (k = 0; k < z; k++) {

                        if (dom[i][j][k] == 0) {
                            pressure_integrated[i][j][k] += pow((mean_pressure - pre[i][j][k]), 2) * dt;
                        }
                    }
                }
            }

#ifdef _OPENMP
        }
#endif
    }

    //here the preparation for the pulse dB - Map is done
    if (with_one_pulse == 1) {

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k, local_pressure_left, \
                local_pressure_right, one_boundary, \
                left, right, tmp) \
                shared(x, y, z, pre, pre_old, dom, \
                pressure_integrated, mean_pressure, dt, \
                perfect_sinus_emitter)
        {
#endif 

            for (i = 0; i < x; i++) {
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
                for (j = 0; j < y; j++) {
                    for (k = 0; k < z; k++) {
                        if (dom[i][j][k] == 0) {

                            tmp = fabs(pre_old[i][j][k] - pre[i][j][k]) * perfect_sinus_emitter;

                            if (tmp < 2.0e-5) tmp = 0.0;
                            if (tmp > pressure_integrated[i][j][k]) pressure_integrated[i][j][k] = tmp;
                        }
                    }
                }
            }
#ifdef _OPENMP
        }
#endif
    }

    return 0;

}

int reset_pressure_integrated(int x, int y, int z) {

    int i, j, k;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, pressure_integrated)
    {
#endif      

        for (i = 0; i < x; i++) {
#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif

            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {

                    pressure_integrated[i][j][k] = 0.0;
                }
            }
        }

#ifdef _OPENMP
    }
#endif        

    return 0;

}

int prepare_the_dB_map(int x, int y, int z) {

    int i, j, k;

#ifdef _OPENMP   
#pragma omp parallel default(none) \
                private(j, i, k) \
                shared(x, y, z, pressure_integrated, \
                dB_map, mean_pressure, dt_integrated, \
                with_one_pulse, dB_Map_ready)
    {
#endif      

        //just to make sure 0.0 dB-Maps are not written out
        dB_Map_ready = 0;


#ifdef _OPENMP               
#pragma omp for schedule(static)
#endif
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                for (k = 0; k < z; k++) {

                    if (pressure_integrated[i][j][k] > 0.0) {
                        if (with_one_pulse == 0) dB_map[i][j][k] = 20 * log10((sqrt((1.0 / dt_integrated) * pressure_integrated[i][j][k])) / 2.0E-5);
                        //Assumption perfect sinus sound emitter
                        if (with_one_pulse == 1) dB_map[i][j][k] = 20 * log10(pressure_integrated[i][j][k] / 2.0E-5);

                        dB_Map_ready = 1;

                    } else {
                        dB_map[i][j][k] = 0.0;
                    }
                }
            }
        }

#ifdef _OPENMP
    }
#endif        

    return 0;
}