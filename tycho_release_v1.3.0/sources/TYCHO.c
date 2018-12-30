/*
 Name        : TYCHO.c
 Author      : Wolfgang Kapferer
 Version     : 1.3.0
 Copyright   : free software
 Release     : October 2013
 Description :
 TYCHO is a multidimensional (1D/2D and 3D) compressible hydrodynamics code written in C and parallelized with OpenMP.
 A Lagrangian remap version of the Piecewise Parabolic Method developed by Paul Woodward and Phil Colella (1984) is
 applied. The code is based on the freely available VH-1 Package VH-1 (http://wonka.physics.ncsu.edu/pub/VH-1/index.php). 
 The simulation package is focused on wind-flow experiments with obstacles in it and advection of marker fields for investigating
 obstacle-gas interactions. In addition momenta and their direction on obstacles are calculated. Thermal diffusion on obstacles, 
 thermal-exchange between obstacles and the surrounding gas and a Sutherland's-law viscosity routine can be switched on in.
 Gravity is included as a constant background-field and a stratified atmosphere can be set up, if needed.
 Tycho is now able to calculate sound-pressure maps of the computational from a simulation with a sound emitters.
 

 TYCHO is freely available to everyone. You are welcome to download it and do whatever you want with it.
 
 
 ATTENTION THE CODE USES SI UNITS
 ATTENTION THE CODE USES SI UNITS
 ATTENTION THE CODE USES SI UNITS
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 Header Files with the global variables  	
 */
#include "variables_global.h"
#include "prototypes.h"

/*
To determine the maximum of two integers
 */
inline int max_array(int a, int b) {
    double tmp;

    if (a < b) tmp = b;
    if (a == b)tmp = b;
    if (a > b) tmp = a;

    return tmp;

}

/*! 
 * This is the main driver of the simulation program.
 * In this function the simulation is initiated,
 * initial conditions are read in or calculated and
 * the main simulation loop is set up. The hydrodynamic
 * simulation is realized with 1D sweeps in all direction 
 * and  forward/ backward. In addition diffusion and
 * obstacle physics routines are called, when stated
 * in the parameter file.
 */

int main(int argc, char *argv[]) {

    int tmp1;
    int i, j, k;

#ifndef _OPENMP
    // For performance analysis issues.
    clock_t cputime_sim1, cputime_sim2;
    clock_t cputime_sim1_eta, cputime_sim2_eta;
#endif   

#ifdef _OPENMP
    // For performance analysis issues.
    double cputime_sim1, cputime_sim2;
    double cputime_sim1_eta, cputime_sim2_eta;
#endif   


    specific_heat_capacity_gas = 0.0;
    specific_heat_capacity_obstacle = 0.0;

    double seconds, seconds_eta;

    // File for loading simulation parameter from parameter file.
    char filename[200];

    restart = 0;

    // How many second has a minute ;-)
    minute_in_sec = 60.0;

    //for sound-emitter simulations some initial values
    with_sound_emitter = 0;
    sound_pressure_level = 0.0;
    sound_frequency = 0.0;
    sound_reflexion = 1.0;
    obstalce_absorption_coefficient = 0.0;
    dt_integrated = 0.0;
    dB_Map_ready = 0;
    // Parameterfilename as argument.
    if (argc > 1) {
        sprintf(filename, "%s", argv[1]);
        start_file_reader(filename);
    }

    if (argc == 1) {
        printf("-TYCHO SIMULATION--\n");
        printf("-Linux  Version  --\n");
        printf("-V1.3.0 Oct. 2013--\n");
        printf("-------------------\n");
        printf("Please specify a\n");
        printf("parameter File\n");
        printf("-------------------\n");
        exit(11);
    }

    // restart logic
    if (argc == 3) {
        restart = 1;
        printf("---------------This is a TYCHO Simulation--------------\n");
        printf("---------------------Linux Version --------------------\n");
        printf("---------------Version 1.3.0 released Oct. 2013--------\n");
        printf("The code will restart from last checkpoint\n");
    }

    /*
    Output on the standardout at programme start
     */
    printf("---------------This is a TYCHO Simulation--------------\n");
    printf("---------------------Linux ----------------------------\n");
    printf("---------------Version 1.3.0 released Oct. 2013--------\n");
    printf("The resolution of the simulation is\n%i %i %i\n", x, y, z);
    printf("The domain extend in meter in x,y and z is\n%g %g %g\n", xmax,
            ymax, zmax);
    //to avoid zero time steps
    if ((xmax == 0.0) && (xmax == 0.0) && (zmax == 0.0)) {
        printf("Attention the extend of your computational domain is a singularity!\n");
        printf("We better stop.\n");
        exit(0);
    }
    printf("Your boundary conditions are\n%i %i %i %i %i %i\n",
            bound.down, bound.up, bound.left, bound.right, bound.front, bound.back);
    printf("Your simulation endtime in minutes is\n%g\n", endtime_sim / minute_in_sec);
    if (gravity_on_off == 0) printf("Gravity is off\n");
    if (gravity_on_off == 1) printf("Gravity is on\n");
    if (wind_on_off == 1) printf("A wind marker file will be processed\n");
    if (with_obstacles == 0) printf("No Obstacles within computational domain\n");
    if (with_sound_emitter == 1) printf("With sound emitter\n");
    if (with_one_pulse == 1) printf("With a single sound-shock emitted\n");
    if (advection == 1) printf("A marker field will be advected with the HYDRO solution\n");
    if (viscosity_on_off == 1) printf("Viscosity is switched on\n");
    if (obstacle_heat_conductivity != 0.0) printf("Heat Diffusion is switched on\n");
#ifdef _OPENMP
    printf("You are using %i threads in parallel\n", number_of_threads);
#endif
    printf("-------------------------------------------------------\n");

    //To get output on the standard out
    fflush(stdout);

#ifdef _OPENMP
    // Sets the number of threads specified in the parameter-file.
    omp_set_num_threads(number_of_threads);
#endif

    // Get the maximum number of cells in all three dimension.
    tmp1 = max_array(x, y);
    max_array_length = max_array(tmp1, z);

    // The minimum and maximum index of the cells without ghost cells.
    nminx = 6;
    nmaxx = nminx + x - 1;
    nminy = 6;
    nmaxy = nminy + y - 1;
    nminz = 6;
    nmaxz = nminz + z - 1;

    //Cell Volume
    if (dimension == 3) cell_volume = (xmax / x)*(ymax / y)*(zmax / z);
    //Cell Area
    if (dimension == 2) cell_volume = (xmax / x)*(ymax / y);
    //Cell length
    if (dimension == 1) cell_volume = (xmax / x);

    // Let's start the simulation at the beginning of time.
    time_sim = 0;

    // If Gravitational acceleration the famous constant is defined here.
    if (gravity_on_off == 1) {
        grav_acc = -9.81;
    } else {
        grav_acc = 0.0;
    }

    // Gamma for gas physics.
    Gamma = Gamma1 - 1.0;

    // Small numbers of error handling.
    small = 1.0E-50;
    smallr = 1.0E-50;
    smallp = 1.0E-50;

    //C1_visc and S_visc constant for determine the viscosity with Sutherland's law
    C1_visc = 1.456E-6;
    S_visc = 110.4;

    // The point of origin.
    xmin = 0;
    ymin = 0;
    zmin = 0;

    counter = 0;
    counter_restart = 0;

    //
    intial_soundspeed = 0.0;

    // Dynamic allocation for global variables and arrays.
    callocate_arrays_global(x, y, z, max_array_length);

    //To get output on the standard out
    fflush(stdout);

    // Initiate variables for the hydro solver.
    initiate_grid(x, xmin, xmax, y, ymin, ymax, z, zmin, zmax);

    //To get output on the standard out
    fflush(stdout);

    // Make Initial Conditions or read the from a file in or read only the obstacle domain.
    if (restart == 0) {

        if (advection == 1) {
            read_marker_file(x, y, z);
            initiate_domain_marker(x, y, z);
        }

        if (wind_on_off == 1) {
            read_wind_file(x, y, z);
        }

        if (make_ics == 0) make_ic(x, y, z);
        if (make_ics == 1) read_ic(x, y, z);
        if (make_ics == 2) {
            make_ic(x, y, z);
            read_dom(x, y, z);
        }
        if (make_ics == 3) make_sod_ic(x, y, z);
        if (make_ics == 4) make_kh_instabilities(x, y, z);

        // Initiate variables for the computational domain if obstacles are involved.
        if ((with_obstacles == 1) && (make_ics == 0)) initiate_domain(x, y, z);
        if ((with_obstacles == 1) && (make_ics == 2)) initiate_domain_dom(x, y, z);

        //read in the sound emitter file
        if (with_sound_emitter == 1) read_soundemitter(x, y, z);

        // The initial conditions are written into files.
        if ((filetype == 0) && ((make_ics == 0) || (make_ics == 2) || (make_ics == 3) || (make_ics == 4))) write_ic_tyc(x, y, z);
        if ((filetype == 1) && ((make_ics == 0) || (make_ics == 2) || (make_ics == 3) || (make_ics == 4))) write_ic_vtk(x, y, z);
        if ((filetype == 2) && ((make_ics == 0) || (make_ics == 2) || (make_ics == 3) || (make_ics == 4))) write_ic_amira(x, y, z);
        if ((filetype == 3) && ((make_ics == 0) || (make_ics == 2) || (make_ics == 3) || (make_ics == 4))) write_ic_ifrit(x, y, z);

    }

    if (restart == 1) {

        read_restart(x, y, z);
        printf("Restart files are read.\n");


        //To get output on the standard out
        fflush(stdout);
    }

    // Here the actual hydro solver starts.
    // calculate a first time-step
    dt_calc_first(x, y, z);

    //In the case of a harmonic soundemitter, we have to calculate a new extend, 
    //because of the frequency 100 cells per lambda is quite okay
    //if the extent of the computational domain is large enough this given extend
    //will not be changed
    //if just one soundpulse is emitted (i.e. a weak shock), this is not necessary
    if ((with_sound_emitter == 1) && (with_one_pulse == 0)) {

        double x_max_new, y_max_new, z_max_new;

        x_max_new = xmax;
        y_max_new = ymax;
        z_max_new = zmax;

        double wavelength_divisor = 100.0;
        x_max_new = (intial_soundspeed / sound_frequency) / wavelength_divisor * x;
        if (dimension > 1) y_max_new = (intial_soundspeed / sound_frequency) / wavelength_divisor * y;
        if (dimension > 2) z_max_new = (intial_soundspeed / sound_frequency) / wavelength_divisor * z;

        if ((x_max_new < xmax) || (y_max_new < ymax) || (z_max_new < zmax)) {
            xmax = x_max_new;
            ymax = y_max_new;
            zmax = z_max_new;
            printf("New extend due to harmonic soundsource\nextend in x-direction: %g\nextend in y-direction: %g\nextend in z-direction: %g\n", xmax, ymax, zmax);
            printf("Intial_soundspeed %g [m/s]\n", intial_soundspeed);
            printf("-------------------------------------------------------\n");

            dt_calc_first(x, y, z);
        }
    }

    printf("The beginning time-step: %g [s]\n", dt);
    printf("-------------------------------------------------------\n");

    //To get output on the standard out
    fflush(stdout);

    // For the shock detection we set the variable shock_or_not to zero.
    shock_or_not = 0;

    //the present pressure are copied into the pre_old array
    if ((with_sound_emitter == 1) && (restart == 0)) reset_pressure_integrated(x, y, z);
    if (with_sound_emitter == 1) pre_old_copy(x, y, z);
    if (with_sound_emitter == 1) noise_generator(x, y, z);
    // For performance analysis issues.
    cputime_sim1 = clock();
#ifdef _OPENMP
    cputime_sim1 = omp_get_wtime();
#endif

    // The main Hydro loop, which runs until counter is smaller then a certain integer.
    do {

        // For performance analysis issues.
        cputime_sim1_eta = clock();
#ifdef _OPENMP
        cputime_sim1_eta = omp_get_wtime();
#endif

        //Noise source, i.e. periodic pressure changes
        if (with_sound_emitter == 1) noise_generator(x, y, z);

        //the present pressure are coped into the pre_old array
        if (with_sound_emitter == 1) pre_old_copy(x, y, z);

        // The time-step is now copied to the old time-step.
        olddt = dt;

        // The time-step due to the CFL criteria is determined.
        dt_calc(x, y, z);

        // the arrays where the pressure on the solid is stored 
        // are reseted here
        pressure_on_solid_reset(x, y, z);

        /*
        The Hydro Sweeps forward: In 1D only in x, in 2D in x and y, and
        finally in 3D in the x,y and z direction.
         */
        hydro_sweeps(x, y, z, 0);

        //if heat diffusion is set on
        if (obstacle_heat_conductivity != 0.0) diffusion(x, y, z);

        //Noise source, i.e. periodic pressure changes
        if (with_sound_emitter == 1) noise_generator(x, y, z);

        //the present pressure are coped into the pre_old array
        if (with_sound_emitter == 1) pre_old_copy(x, y, z);
        // The time-step due to the CFL criteria is determined.
        dt_calc(x, y, z);

        // the arrays where the pressure on the solid is stored 
        // are reseted here
        pressure_on_solid_reset(x, y, z);

        /*
         The Hydro Sweeps backward: In 1D only in x, in 2D in x and y, and
         finally in 3D in the x,y and z direction.
         */
        hydro_sweeps(x, y, z, 1);

        // if heat diffusion is set on
        if (obstacle_heat_conductivity != 0.0) diffusion(x, y, z);

        // Velocity Field Analyser
        velocity_field_analyser(x, y, z);

        //Noise source, i.e. periodic pressure changes
        if (with_sound_emitter == 1) noise_generator(x, y, z);

        //the present pressure are coped into the pre_old array
        if (with_sound_emitter == 1) pre_old_copy(x, y, z);

        // Calculate the total pressure on the solid in the computational domain.
        if (with_obstacles == 1) pressure_on_solid_calc(x, y, z);

        // What time is it?
        time_sim = time_sim + 2 * dt;

        // For performance analysis issues.
        cputime_sim2_eta = clock();
#ifdef _OPENMP
        cputime_sim2_eta = omp_get_wtime();
#endif

        seconds_eta = (cputime_sim2_eta - cputime_sim1_eta) / (double) CLOCKS_PER_SEC;
        seconds_eta /= 2.0;

#ifdef _OPENMP
        seconds_eta = (cputime_sim2_eta - cputime_sim1_eta);
        seconds_eta /= 2.0;
#endif
        //check if we have to reset the pressure_integrated array
        //and store the stuff in the dba-map;
        dt_integrated += dt;

        if ((dt_integrated >= 1.0 / sound_frequency) && (with_sound_emitter == 1) && (with_one_pulse == 0)) {

            //make a db map;
            prepare_the_dB_map(x, y, z);
            //reset stuff for the generation of the next dB map
            dt_integrated = 0.0;
            reset_pressure_integrated(x, y, z);
        }

        if ((with_sound_emitter == 1) && (with_one_pulse == 1)) {
            //make a db map;
            prepare_the_dB_map(x, y, z);
        }

        // Write output files at a given frequency
        if (time_sim > output_frequency * counter) {
            if (filetype == 0) write_tyc(x, y, z, counter);
            if (filetype == 1) write_vtk(x, y, z, counter);
            if (filetype == 2) write_amira(x, y, z, counter);
            if (filetype == 3) write_ifrit(x, y, z, counter);

            // For performance analysis issues.
            cputime_sim2 = clock();
#ifdef _OPENMP
            cputime_sim2 = omp_get_wtime();
#endif

            seconds = (cputime_sim2 - cputime_sim1) / (double) CLOCKS_PER_SEC;

#ifdef _OPENMP
            seconds = (cputime_sim2 - cputime_sim1);
#endif

            printf("%3.1f percent done\n", 100.0 / (float) endtime_sim * (float) time_sim);
            printf("time:%g [s], dt:%g [s]\nshock:%i\ncounter:%i\n", time_sim, dt, shock_or_not, counter);

            if (counter > 0) {
                if ((((endtime_sim - time_sim) / dt) * seconds_eta) < 60.0) {
                    printf("Estimated time to simulation end: %4.2f seconds\n", ((endtime_sim - time_sim) / dt) * seconds_eta);
                }
                if (((((endtime_sim - time_sim) / dt) * seconds_eta) > 60.0) && (((((endtime_sim - time_sim) / dt) * seconds_eta) <= 3600.0))) {
                    printf("Estimated time to simulation end: %4.2f minutes\n", ((endtime_sim - time_sim) / dt) * seconds_eta / 60.0);
                }
                if (((((endtime_sim - time_sim) / dt) * seconds_eta) > 3600.0) && ((((endtime_sim - time_sim) / dt) * seconds_eta) <= 3600.0 * 24.0)) {
                    printf("Estimated time to simulation end: %5.2f hours\n", ((endtime_sim - time_sim) / dt) * seconds_eta / 3600.0);
                }
                if ((((endtime_sim - time_sim) / dt) * seconds_eta) > 3600.0 * 24.0) {
                    printf("Estimated time to simulation end: %5.2f days\n", ((endtime_sim - time_sim) / dt) * seconds_eta / (3600.0 * 24.0));
                }
                if (seconds < 60.0) printf("Time since last output: %4.2g [s]\n", seconds);
                if ((seconds > 60.0) && (seconds < 3600)) printf("Time since last output: %4.2g [min]\n", seconds / 60.0);
                if (seconds > 3600.0) printf("Time since last output: %4.2g [h]\n", seconds / 3600.0);
            }
            printf("-------------------------------------------------------\n");

            if (time_sim > restart_frequency * counter_restart) {
                if (filetype != 0) {
                    write_restart_tyc(x, y, z);
                    printf("Restart files are written.\n");
                    counter_restart++;
                }
            }

            counter++;

            //To get output on the standard out
            fflush(stdout);
            // For performance analysis issues.

            cputime_sim1 = clock();
#ifdef _OPENMP
            cputime_sim1 = omp_get_wtime();
#endif
        }


    } while (time_sim <= endtime_sim);

    printf("\n--------------------------");
    printf("\n--------------------------\n");
    printf("Program aborted normally\n");
    printf("End-time of simulation reached\n");
    printf("\n--------------------------");
    printf("\n--------------------------\n");

    return 0;
}
