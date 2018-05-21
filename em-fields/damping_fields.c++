#include "damping_fields.h"


fields::PlasmaCellDamped::PlasmaCellDamped(
    size_t i, size_t j,
    int o,
    size_t NxG, size_t NyG,
    size_t NxMesh, size_t NyMesh, size_t NzMesh
    ) : 
  corgi::Cell(i, j, o, NxG, NyG),
  PlasmaCell(i, j, o, NxG, NyG, NxMesh, NyMesh, NzMesh)
{ }


void fields::PlasmaCellDamped::pushE() {

  // this->pushE1d();
  this->pushE2d_damped();
  // this->pushE3d();
}


/// 2D E pusher
void fields::PlasmaCellDamped::pushE2d_damped() {

  fields::YeeLattice& mesh = getYee();

  //std::cout << "Calling DAMPED E update\n";

 
  int k = 0;
  for(int j=0; j<(int)NyMesh; j++) {
    for(int i=0; i<(int)NxMesh; i++) {

      // Ex
      mesh.ex(i,j,k) += 
        + yeeDt*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k)) / yeeDx;

      // Ey
      mesh.ey(i,j,k) += 
        + yeeDt*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k)) / yeeDx;

      // Ez
      mesh.ez(i,j,k) += 
        + yeeDt*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k)) / yeeDx
        + yeeDt*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k)) / yeeDx;

    }
  }
}


void fields::PlasmaCellDamped::depositCurrent() {
  fields::YeeLattice& mesh = getYee();

  //std::cout << "Calling DAMPED J update\n";

  //std::cout<<"dt:"<<yeeDt<<"  and vol:"<<yeeDx<<" .. " <<(yeeDx*yeeDy*yeeDz) <<"\n";
  
  Realf resistivity = 10.0;

  mesh.ex -= mesh.jx * yeeDt / resistivity;
  mesh.ey -= mesh.jy * yeeDt / resistivity;
  mesh.ez -= mesh.jz * yeeDt / resistivity;

}




