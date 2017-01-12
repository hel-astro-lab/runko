#ifndef PARAMS_H
#define PARAMS_H

#include "vector"

#include "definitions.hpp"
#include "common.h"

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
    static bool Nx_wrap;    ///< wrap grid in X dir
    static bool Ny_wrap;    ///< wrap grid in y dir
    static bool Nz_wrap;    ///< wrap grid in z dir

    // AMR info
    static uint64_t N_neighb; ///<  size of the ghost cell neighborhood
    static uint64_t max_ref_lvl; ///< maximum AMR refinement level

    static double dt; ///< Current time-step

    static std::vector<CellID> localCells; ///< Cached copy of spatial cell IDs on this process.

};



#endif

