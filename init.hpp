/* Physical and simulation related parameters

*/

#ifndef INIT_HPP
#define INIT_HPP

// total particle number
const uint64_t Np = 1e1;


// grid size
#ifdef DEBUG
const uint64_t Nx = 6;
const uint64_t Ny = 1;
const uint64_t Nz = 1;
#else
const uint64_t Nx = 100;
const uint64_t Ny = 20;
const uint64_t Nz = 20;
#endif

// real spatial dimensions
double grid_xmin = 0.0;
double grid_xmax = 10.0;

double grid_ymin = 0.0;
double grid_ymax = 1.0;

double grid_zmin = 0.0;
double grid_zmax = 1.0;

// periodicity
const bool Nx_wrap = true;
const bool Ny_wrap = false;
const bool Nz_wrap = false;

// AMR info
const uint64_t N_neighb = 1;
const uint64_t max_ref_lvl = 0;

// Neighborhoods
//-------------------------------------------------- 
const int Cp1_shift  = 10;
const int Cm1_shift  = 11;
// const int four_current_shift = 11; // not implemented
// const int nodal_field_shift  = 13;


// physical stuff
//-------------------------------------------------- 
const double c = 10.0;
const double q = 1.0;
const double e = 1.0;
const double pi = 3.14159265359;


const double grid_dx = (grid_xmax - grid_xmin)/Nx;
const double grid_dy = (grid_ymax - grid_ymin)/Ny;
const double grid_dz = (grid_zmax - grid_zmin)/Nz;
const double dt = 0.99/sqrt(1.0/grid_dx/grid_dx + 1.0/grid_dy/grid_dy + 1.0/grid_dz/grid_dz)/c;

const double me = 1.0;
const double mp = 16.0;



#endif
