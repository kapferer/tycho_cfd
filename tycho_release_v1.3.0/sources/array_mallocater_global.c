/*
 * array_callocater_global.c
 *
 * Author: Wolfgang Kapferer 
 */

#include <stdio.h>
#include <stdlib.h>
#include "variables_global.h"

/*!
 The main hydro arrays allocation routine.
 Here density, pressure, velocities and
 the computational domain arrays are dynamically
 allocated. The size is determined from
 the resolution parameters x,y and z in
 the specified parameter file.
 */
int callocate_arrays_global(int x, int y, int z) {

    int i, j;

    // allocate density
    rho = calloc(x, sizeof *rho);

    for (i = 0; i < x; i++) {
        rho[i] = calloc(y, sizeof **rho);
        for (j = 0; j < y; j++) {
            rho[i][j] = calloc(z, sizeof ***rho);
        }
    }

    printf("Density array memory allocated\n");

    // allocate pressure
    pre = calloc(x, sizeof * pre);

    for (i = 0; i < x; i++) {
        pre[i] = calloc(y, sizeof **pre);
        for (j = 0; j < y; j++) {
            pre[i][j] = calloc(z, sizeof ***pre);
        }
    }

    printf("Pressure array memory allocated\n");

    // allocate domain
    dom = calloc(x, sizeof *dom);

    for (i = 0; i < x; i++) {
        dom[i] = calloc(y, sizeof **dom);
        for (j = 0; j < y; j++) {
            dom[i][j] = calloc(z, sizeof ***dom);
        }
    }

    // allocate pressure on solid
    pressure_on_solid = calloc(x, sizeof *pressure_on_solid);

    for (i = 0; i < x; i++) {
        pressure_on_solid[i] = calloc(y, sizeof **pressure_on_solid);
        for (j = 0; j < y; j++) {
            pressure_on_solid[i][j] = calloc(z, sizeof ***pressure_on_solid);
        }
    }


    printf("Domain array and pressure on solid memory allocated\n");

    // allocate energy
    eng = calloc(x, sizeof * eng);

    for (i = 0; i < x; i++) {
        eng[i] = calloc(y, sizeof **eng);
        for (j = 0; j < y; j++) {
            eng[i][j] = calloc(z, sizeof ***eng);
        }
    }

    e_in = calloc(x, sizeof * e_in);

    for (i = 0; i < x; i++) {
        e_in[i] = calloc(y, sizeof **e_in);
        for (j = 0; j < y; j++) {
            e_in[i][j] = calloc(z, sizeof ***e_in);
        }
    }

    printf("Energy and Temperature array memory allocated\n");

    // allocate velocities
    vx = calloc(x, sizeof * vx);

    for (i = 0; i < x; i++) {
        vx[i] = calloc(y, sizeof **vx);
        for (j = 0; j < y; j++) {
            vx[i][j] = calloc(z, sizeof ***vx);
        }
    }

    vy = calloc(x, sizeof * vy);

    for (i = 0; i < x; i++) {
        vy[i] = calloc(y, sizeof **vy);
        for (j = 0; j < y; j++) {
            vy[i][j] = calloc(z, sizeof ***vy);
        }
    }

    vz = calloc(x, sizeof * vz);

    for (i = 0; i < x; i++) {
        vz[i] = calloc(y, sizeof **vz);
        for (j = 0; j < y; j++) {
            vz[i][j] = calloc(z, sizeof ***vz);
        }
    }

    printf("Velocities array memory allocated\n");

    if (advection == 1) {
        // allocate marker
        marker = calloc(x, sizeof *marker);

        for (i = 0; i < x; i++) {
            marker[i] = calloc(y, sizeof **marker);
            for (j = 0; j < y; j++) {
                marker[i][j] = calloc(z, sizeof ***marker);
            }
        }

        printf("Marker array memory allocated\n");
    }

    if (wind_on_off == 1) {
        // allocate marker
        wind_marker = calloc(x, sizeof *wind_marker);

        for (i = 0; i < x; i++) {
            wind_marker[i] = calloc(y, sizeof **wind_marker);
            for (j = 0; j < y; j++) {
                wind_marker[i][j] = calloc(z, sizeof ***wind_marker);
            }
        }

        printf("Wind_marker array memory allocated\n");
    }



    //only if viscosity is switched on 
    if (viscosity_on_off == 1) {
        // allocate density
        rho_visc = calloc(x, sizeof *rho_visc);

        for (i = 0; i < x; i++) {
            rho_visc[i] = calloc(y, sizeof **rho_visc);
            for (j = 0; j < y; j++) {
                rho_visc[i][j] = calloc(z, sizeof ***rho_visc);
            }
        }

        // allocate velocities
        vx_visc = calloc(x, sizeof * vx_visc);

        for (i = 0; i < x; i++) {
            vx_visc[i] = calloc(y, sizeof **vx_visc);
            for (j = 0; j < y; j++) {
                vx_visc[i][j] = calloc(z, sizeof ***vx_visc);
            }
        }

        vy_visc = calloc(x, sizeof * vy_visc);

        for (i = 0; i < x; i++) {
            vy_visc[i] = calloc(y, sizeof **vy_visc);
            for (j = 0; j < y; j++) {
                vy_visc[i][j] = calloc(z, sizeof ***vy_visc);
            }
        }

        vz_visc = calloc(x, sizeof * vz_visc);

        for (i = 0; i < x; i++) {
            vz_visc[i] = calloc(y, sizeof **vz_visc);
            for (j = 0; j < y; j++) {
                vz_visc[i][j] = calloc(z, sizeof ***vz_visc);
            }
        }

        printf("viscosity arrays are allocated\n");
    }

    if (with_sound_emitter == 1) {
        // allocate soundemitter
        soundemitter = calloc(x, sizeof *soundemitter);
        for (i = 0; i < x; i++) {
            soundemitter[i] = calloc(y, sizeof **soundemitter);
            for (j = 0; j < y; j++) {
                soundemitter[i][j] = calloc(z, sizeof ***soundemitter);
            }
        }
        // allocate old pressure
        pre_old = calloc(x, sizeof *pre_old);

        for (i = 0; i < x; i++) {
            pre_old[i] = calloc(y, sizeof **pre_old);
            for (j = 0; j < y; j++) {
                pre_old[i][j] = calloc(z, sizeof ***pre_old);
            }
        }

        // allocate pressure_integrated
        pressure_integrated = calloc(x, sizeof *pressure_integrated);

        for (i = 0; i < x; i++) {
            pressure_integrated[i] = calloc(y, sizeof **pressure_integrated);
            for (j = 0; j < y; j++) {
                pressure_integrated[i][j] = calloc(z, sizeof ***pressure_integrated);
            }
        }

        // allocate dB_map
        dB_map = calloc(x, sizeof *dB_map);

        for (i = 0; i < x; i++) {
            dB_map[i] = calloc(y, sizeof **dB_map);
            for (j = 0; j < y; j++) {
                dB_map[i][j] = calloc(z, sizeof ***dB_map);
            }
        }

        printf("dB-Map allocated\n");

        printf("Soundemitter, pressure old array, pressure integrated and dB-Map memory allocated\n");
    }

    // allocate arrays for grid initiate
    zxa = calloc(x, sizeof * zxa);
    zxc = calloc(x, sizeof * zxc);
    zdx = calloc(x, sizeof * zdx);

    zya = calloc(y, sizeof * zya);
    zyc = calloc(y, sizeof * zyc);
    zdy = calloc(y, sizeof * zdy);

    zza = calloc(z, sizeof * zza);
    zzc = calloc(z, sizeof * zzc);
    zdz = calloc(z, sizeof * zdz);


    return 0;
}