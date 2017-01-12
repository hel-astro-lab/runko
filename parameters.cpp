#include "parameters.hpp"

typedef Parameters P;


uint64_t P::Np = 1e1;

#ifdef DEBUG
uint64_t P::Nx = 6;
uint64_t P::Ny = 1;
uint64_t P::Nz = 1;
#else
uint64_t P::Nx = 100;
uint64_t P::Ny = 20;
uint64_t P::Nz = 20;
#endif

// real spatial dimensions
double P::grid_xmin = 0.0;
double P::grid_xmax = 10.0;

double P::grid_ymin = 0.0;
double P::grid_ymax = 1.0;

double P::grid_zmin = 0.0;
double P::grid_zmax = 1.0;

double P::grid_dx = (P::grid_xmax - P::grid_xmin)/P::Nx;
double P::grid_dy = (P::grid_ymax - P::grid_ymin)/P::Ny;
double P::grid_dz = (P::grid_zmax - P::grid_zmin)/P::Nz;

bool P::Nx_wrap = true;
bool P::Ny_wrap = false;
bool P::Nz_wrap = false;

double P::dt = 0.99/std::sqrt(1.0/P::grid_dx/P::grid_dx + 1.0/P::grid_dy/P::grid_dy + 1.0/P::grid_dz/P::grid_dz)/c;


uint64_t P::N_neighb = 1;
uint64_t P::max_ref_lvl = 0;

