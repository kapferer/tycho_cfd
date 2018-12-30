/*
 * variables.h
 *
 * Author: Wolfgang Kapferer
 */

#ifndef VARIABLES_global_
#define VARIABLES_global_

double grav_acc;

/* Data types */
typedef enum {
    DOM_FLUID = 0,
    DOM_SOLID
} dom_state;


/* Welcome to variables home*/

//restarting the code or not
int restart;

// the resolution
int x, y, z, dimension;

int max_array_length;

// for the grid the nmin, nmax
int nminx, nmaxx;
int nminy, nmaxy;
int nminz, nmaxz;

int shock_or_not, counter, counter_restart;

double Gamma, Gamma1;

//error catch
double smallp, smallr, small;

double xmin, xmax, ymin, ymax, zmin, zmax;

//volume of one grid cell
double cell_volume;

//timing
double dt, dt_no_shock, dt_shock;
double time_sim, endtime_sim;
double minute_in_sec;
double dt_integrated;

//the old time_sim step;
double olddt;

//the cfl parameter
double courant;

//the 3D arrays for hydro
double ***rho, ***pre, ***eng, ***vx, ***vy, ***vz, ***e_in;

//the 3D arrays for viscosity calculation
double ***rho_visc, ***vx_visc, ***vy_visc, ***vz_visc;

//an advectable gas
double ***marker;

//wind marker
double ***wind_marker;

//the 3D array for the domain
int ***dom;

//the 3D array for the sound_emitter
int ***soundemitter;

//pressure on solid and old pressure array
double ***pressure_on_solid, ***pre_old;

//for the dB-map
double ***pressure_integrated, ***dB_map;

// init grid
double *zdx, *zxc, *zxa;
double *zdy, *zyc, *zya;
double *zdz, *zzc, *zza;

double spacing;

//boundary condition switch
typedef struct {
    int up;
    int down;

    int left;
    int right;

    int front;
    int back;

} boundary;

boundary bound;

char dens_ic[200], temp_ic[200], vel_ic[200], obst_ic[200];
char sound_ic[200], marker_ic[200], wind_ic[200];
char output_dir[200];

int intial_velocity_file;

double output_frequency;
double restart_frequency;

int gravity_on_off;

int viscosity_on_off;

double C1_visc;
double S_visc;

int wind_on_off;
double wind_speed;
int wind_direction;

int make_ics;

int filetype;

double inflow_velocity;
double inflow_density;
double inflow_temperature;
int starting_flow;

double gasconstant;

int number_of_threads;

//if obstacles are included
int with_obstacles;

double obstacle_density, obstacle_temperature, obstacle_heat_conductivity;

//for the heat-transfer calculation at the boundary of obstacles
double specific_heat_capacity_gas;
double specific_heat_capacity_obstacle;

int strat_const_atmos;

//if advection of a maker field is desired
int advection;

double marker_density;

//for the sound generator the initial pressure at the position of the sound source
int with_sound_emitter;
double pre_initial_sound;
double sound_pressure_level;
int with_one_pulse;
double sound_frequency;
double sound_reflexion;
double obstalce_absorption_coefficient;
double mean_pressure;
double pressure;

double intial_soundspeed;

int dB_Map_ready;

#endif /* VARIABLES_global_ */
