#ifndef PARAMS_H
#define PARAMS_H

#include "definitions.hpp"

struct Parameters {

    // Total number of particles    
    static uint64_t Np;

    // grid size
    static uint64_t Nx;
    static uint64_t Ny;
    static uint64_t Nz;

    // real spatial dimensions
    static double grid_xmin;
    static double grid_xmax;

    static double grid_ymin;
    static double grid_ymax;

    static double grid_zmin;
    static double grid_zmax;

    static double grid_dx;
    static double grid_dy;
    static double grid_dz;

    // periodicity
    static bool Nx_wrap;
    static bool Ny_wrap;
    static bool Nz_wrap;

    // AMR info
    static uint64_t N_neighb;
    static uint64_t max_ref_lvl;

    static double dt;

};



#endif

