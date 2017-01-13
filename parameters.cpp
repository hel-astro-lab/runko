#include "parameters.hpp"

typedef Parameters P;


uint64_t P::Np = 1e6;

uint64_t P::Nx = 10;
uint64_t P::Ny = 10;
uint64_t P::Nz = 10;


// real spatial dimensions
double P::grid_xmin = 0.0;
double P::grid_xmax = 1.0;

double P::grid_ymin = 0.0;
double P::grid_ymax = 1.0;

double P::grid_zmin = 0.0;
double P::grid_zmax = 1.0;

double P::grid_dx = (P::grid_xmax - P::grid_xmin)/P::Nx;
double P::grid_dy = (P::grid_ymax - P::grid_ymin)/P::Ny;
double P::grid_dz = (P::grid_zmax - P::grid_zmin)/P::Nz;

bool P::Nx_wrap = true;
bool P::Ny_wrap = true;
bool P::Nz_wrap = true;

double P::dt = 0.99/std::sqrt(1.0/P::grid_dx/P::grid_dx + 1.0/P::grid_dy/P::grid_dy + 1.0/P::grid_dz/P::grid_dz)/c;


uint64_t P::N_neighb = 1;
uint64_t P::max_ref_lvl = 0;

std::vector<CellID> P::localCells;

