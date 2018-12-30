/*
 * start_file_reader.c
 *
 * Author      : Wolfgang Kapferer
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 The read routine for the parameterfile
 */
int start_file_reader(char filename[]) {

    char data[55][200];
    int l;
    FILE *fd;
    char buffer[200];

    fd = fopen(filename, "r");
    if (fd == NULL) {
        printf("The Parameter File you specified does not exist.\n");
        exit(10);
    }

    l = 0;
    while (feof(fd) == 0) {
        fgets(data[l], 200, fd);
        if ((data[l][0] == '#') || (data[l][0] == '\n')) continue;
        data[l][strlen(data[l]) - 1] = '\0';
        l++;
    }
    fclose(fd);



    filetype = atoi(data[0]);

    make_ics = atoi(data[1]);

    sprintf(dens_ic, "%s", data[2]);

    sprintf(temp_ic, "%s", data[3]);

    sprintf(vel_ic, "%s", data[4]);

    sprintf(wind_ic, "%s", data[5]);

    sprintf(obst_ic, "%s", data[6]);
    
    sprintf(sound_ic, "%s", data[7]);

    sprintf(marker_ic, "%s", data[8]);

    intial_velocity_file = atoi(data[9]);

    sprintf(output_dir, "%s", data[10]);

    strat_const_atmos = atoi(data[11]);

    x = atoi(data[12]);
    y = atoi(data[13]);
    z = atoi(data[14]);

    dimension = atoi(data[15]);

    xmax = atof(data[16]);
    ymax = atof(data[17]);
    zmax = atof(data[18]);

    bound.down = atoi(data[19]);
    bound.up = atoi(data[20]);
    bound.left = atoi(data[21]);
    bound.right = atoi(data[22]);
    bound.front = atoi(data[23]);
    bound.back = atoi(data[24]);

    inflow_velocity = atof(data[25]);
    inflow_density = atof(data[26]);
    inflow_temperature = atof(data[27]);
    starting_flow = atoi(data[28]);

    endtime_sim = atof(data[29]);
    endtime_sim = endtime_sim*minute_in_sec;

    output_frequency = atof(data[30]);

    restart_frequency = atof(data[31]);

    courant = atof(data[32]);

    gravity_on_off = atoi(data[33]);

    gasconstant = atof(data[34]);

    Gamma1 = atof(data[35]);

    wind_on_off = atoi(data[36]);

    wind_speed = atof(data[37]);

    number_of_threads = atoi(data[38]);

    with_obstacles = atoi(data[39]);

    obstacle_density = atof(data[40]);

    obstacle_temperature = atof(data[41]);
    obstacle_temperature = obstacle_density * gasconstant * obstacle_temperature;

    obstacle_heat_conductivity = atof(data[42]);

    specific_heat_capacity_gas = atof(data[43]);

    specific_heat_capacity_obstacle = atof(data[44]);

    advection = atoi(data[45]);

    marker_density = atof(data[46]);

    viscosity_on_off = atoi(data[47]);
    
    with_sound_emitter = atoi(data[48]);

    sound_pressure_level = atof(data[49]);

    sound_pressure_level = pow(10, (sound_pressure_level / 20))*2.0E-5 * sqrt(2.0);
    
    with_one_pulse = atoi(data[50]);

    sound_frequency = atof(data[51]);
    
    sound_reflexion = atof(data[52]);

    obstalce_absorption_coefficient = atof(data[53]);
    obstalce_absorption_coefficient = pow(10,(obstalce_absorption_coefficient/20));
    obstalce_absorption_coefficient = 1.0/obstalce_absorption_coefficient;
    
    //printf("%g\n", obstalce_absorption_coefficient);

    //for debug
    //printf("%i\n", l);

    if (l != 54) {
        printf("Simulation stops, because of too few data in the start file.\n");
        exit(99);
    }

    return 0;
}
