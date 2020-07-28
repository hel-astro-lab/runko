#include <iostream>
#include <cmath>

#include "tile.h"
#include <nvtx3/nvToolsExt.h> 



namespace fields {

void copy_vert_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nz, int Ny, int halo, int ito, int ifro, int in);

}