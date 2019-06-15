#include "driftcurrent.h"
#include "../../em-fields/tile.h"

#include <cmath>


/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::solve(ffe::Tile<2>& tile)
{
  fields::YeeLattice& mesh = tile.get_yee();

  Realf C = 1.0 * tile.cfl;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Ex
    mesh.ex(i,j,k) += 
      + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k) 
           -mesh.by(i-1,j,   k) + mesh.by(i,j,k));

  }
}

