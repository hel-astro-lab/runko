#include <iostream>
#include <cmath>

#include "tile.h"
#include <nvtx3/nvToolsExt.h> 



namespace fields {

void copy_vert_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nz, int Ny, int halo, int ito, int ifro, int in);
void copy_horz_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nx, int Nz, int halo, int jto, int jfro, int jn);
void copy_face_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nx, int Ny, int halo, int kto, int kfro, int kn);


void copy_x_pencil_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nx, int halo, int jto, int jfro, int kto, int kfro, int jn, int kn);
void copy_y_pencil_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Ny, int halo, int ito, int ifro, int kto, int kfro, int in, int kn);
void copy_z_pencil_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nz, int halo, int ito, int ifro, int jto, int jfro, int in, int jn);

void copy_point_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int halo, int ito, int ifro, int jto, int jfro, int kto, int kfro, int in, int jn, int kn);

}