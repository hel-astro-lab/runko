#include <iostream>

#include "maxwell.h"
#include "tools/Mesh.h"


maxwell::PlasmaCell::PlasmaCell(
    size_t i, size_t j,
    int o,
    size_t NxG, size_t NyG,
    size_t NxMesh, size_t NyMesh
    ) :
  corgi::Cell(i, j, o, NxG, NyG),
  NxMesh(NxMesh), NyMesh(NyMesh), NzMesh(1)
{

  // append fresh Yee lattices into data container
  yee.push_back( maxwell::YeeLattice(NxMesh, NyMesh, NzMesh) );
  // yee.push_back( maxwell::YeeLattice(NxMesh, NyMesh, NzMesh) );



}

/* 
 * 1D version:
 *	ey(i,j,k)=ey(i,j,k) + c *(bz(im1,j,k)-bz(i,j,k))  
 *	ez(i,j,k)=ez(i,j,k) + c *(by(i,j,k)-by(im1,j,k)) 
 *
 * 2D version
 * ex(i,j,k)=ex(i,j,k)+const*(-bz(i,jm1,k)+bz(i,j,k))
 * ey(i,j,k)=ey(i,j,k)+const*(bz(im1,j,k)-bz(i,j,k))
 * ez(i,j,k)=ez(i,j,k)+const*(bx(i,jm1,k)-bx(i,j,k)-  by(im1,j,k)+by(i,j,k))
*/

/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */
void maxwell::PlasmaCell::pushE() {

  // this->pushE1d();
  this->pushE2d();
  // this->pushE3d();
}

/// 1D E pusher
void maxwell::PlasmaCell::pushE1d() {

  maxwell::YeeLattice& mesh = getYee();

  int k = 0;
  int j = 0;
  for(int i=0; i<(int)NxMesh; i++) {

    // Ex
    // NONE

    // Ey
    mesh.ey(i,j,k) += 
      + dt*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k)) / dz;

    // Ez
    mesh.ez(i,j,k) += 
      + dt*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k)) / dy;

  }

}

/// 2D E pusher
void maxwell::PlasmaCell::pushE2d() {

  maxwell::YeeLattice& mesh = getYee();


  int k = 0;
  for(int j=0; j<(int)NyMesh; j++) {
    for(int i=0; i<(int)NxMesh; i++) {

      // Ex
      mesh.ex(i,j,k) += 
        + dt*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k)) / dz;

      // Ey
      mesh.ey(i,j,k) += 
        + dt*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k)) / dz;

      // Ez
      mesh.ez(i,j,k) += 
        + dt*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k)) / dx
        + dt*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k)) / dy;

    }
  }

}


/// 3D E pusher
void maxwell::PlasmaCell::pushE3d() {

  maxwell::YeeLattice& mesh = getYee();

  for(int k=0; k<(int)NzMesh; k++) {
    for(int j=0; j<(int)NyMesh; j++) {
      for(int i=0; i<(int)NxMesh; i++) {

        // Ex
        mesh.ex(i,j,k) += 
          + dt*( mesh.by(i,j,  k-1) - mesh.by(i,j,k)) / dy
          + dt*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k)) / dz;

        // Ey
        mesh.ey(i,j,k) += 
          + dt*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k)) / dz
          + dt*(-mesh.bx(i,  j, k-1) + mesh.bx(i,j,k)) / dx;

        // Ez
        mesh.ez(i,j,k) += 
          + dt*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k)) / dx
          + dt*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k)) / dy;

      }
    }
  }

}




/// Deposit current into electric field
void maxwell::PlasmaCell::depositCurrent() {
  maxwell::YeeLattice& mesh = getYee();

  mesh.ex -= mesh.jx * dt* (dx*dy*dz);
  mesh.ey -= mesh.jy * dt* (dx*dy*dz);
  mesh.ez -= mesh.jz * dt* (dx*dy*dz);

}



/*
 * 1D version:
  by(i,j,k)=by(i,j,k)+const*(ez(ip1,j,k)-ez(i,j,k)) !-0*ex(i,j,k+1)+0*ex(i,j,k))
	bz(i,j,k)=bz(i,j,k)+const*(ey(i,j,k)-ey(ip1,j,k)) !+0*ex(i,j+1,k)-0*ex(i,j,k))
  
 * 2D version:
 * bx(i,j,k)=bx(i,j,k)+const*(-ez(i,jp1,k)+ez(i,j,k))
 * by(i,j,k)=by(i,j,k)+const*(ez(ip1,j,k)-ez(i,j,k))
 * bz(i,j,k)=bz(i,j,k)+const*(ex(i,jp1,k)-ex(i,j,k) -ey(ip1,j,k)+ey(i,j,k))
*/

/// Update B field with a half step
void maxwell::PlasmaCell::pushHalfB() {

  // this->pushHalfB1d();
  this->pushHalfB2d();
  // this->pushHalfB3d();
}

/// 1D B pusher
void maxwell::PlasmaCell::pushHalfB1d() {
  maxwell::YeeLattice& mesh = getYee();

  int k = 0;
  int j = 0;
  for(int i=0; i<(int)NxMesh; i++) {

    // Bx
    // NONE

    // By
    mesh.by(i,j,k) += 
      + dt*0.5*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) / dz;

    // Bz
    mesh.bz(i,j,k) += 
      + dt*0.5*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k)) / dy;
  }

}

/// 2D B pusher
void maxwell::PlasmaCell::pushHalfB2d() {
  maxwell::YeeLattice& mesh = getYee();

  int k = 0;
  for(int j=0; j<(int)NyMesh; j++) {
    for(int i=0; i<(int)NxMesh; i++) {

      // Bx
      mesh.bx(i,j,k) += 
        + dt*0.5*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k)) / dz;

      // By
      mesh.by(i,j,k) += 
        + dt*0.5*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) / dz;

      // Bz
      mesh.bz(i,j,k) += 
        + dt*0.5*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)) / dx
        + dt*0.5*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k)) / dy;

    }
  }

}


/// 3D B pusher
void maxwell::PlasmaCell::pushHalfB3d() {
  maxwell::YeeLattice& mesh = getYee();

  for(int k=0; k<(int)NzMesh; k++) {
    for(int j=0; j<(int)NyMesh; j++) {
      for(int i=0; i<(int)NxMesh; i++) {

        // Bx
        mesh.bx(i,j,k) += 
         + dt*0.5*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k)) / dy
         + dt*0.5*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k)) / dz;

        // By
        mesh.by(i,j,k) += 
         + dt*0.5*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) / dz
         + dt*0.5*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k)) / dx;

        // Bz
        mesh.bz(i,j,k) += 
         + dt*0.5*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)) / dx
         + dt*0.5*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k)) / dy;

      }
    }
  }

}




/// Get current time snapshot of Yee lattice
maxwell::YeeLattice& maxwell::PlasmaCell::getYee() {
  return yee.getRef();
}

/// Get new time snapshot of Yee lattice
maxwell::YeeLattice& maxwell::PlasmaCell::getNewYee() {
  return yee.getNewRef();
};


/// Quick helper function to copy everything inside Yee lattice 
void copyVertYee(
    maxwell::YeeLattice& lhs, 
    maxwell::YeeLattice& rhs, 
    int lhsI, int rhsI) {

  lhs.ex.copyVert(rhs.ex, lhsI, rhsI); 
  lhs.ey.copyVert(rhs.ey, lhsI, rhsI); 
  lhs.ez.copyVert(rhs.ez, lhsI, rhsI); 

  lhs.bx.copyVert(rhs.bx, lhsI, rhsI); 
  lhs.by.copyVert(rhs.by, lhsI, rhsI); 
  lhs.bz.copyVert(rhs.bz, lhsI, rhsI); 

  lhs.jx.copyVert(rhs.jx, lhsI, rhsI); 
  lhs.jy.copyVert(rhs.jy, lhsI, rhsI); 
  lhs.jz.copyVert(rhs.jz, lhsI, rhsI); 

  lhs.rh.copyVert(rhs.rh, lhsI, rhsI); 

}

/// Quick helper function to copy everything inside Yee lattice 
void copyHorzYee(
    maxwell::YeeLattice& lhs, 
    maxwell::YeeLattice& rhs, 
    int lhsJ, int rhsJ) {

  lhs.ex.copyHorz(rhs.ex, lhsJ, rhsJ); 
  lhs.ey.copyHorz(rhs.ey, lhsJ, rhsJ); 
  lhs.ez.copyHorz(rhs.ez, lhsJ, rhsJ); 
                                    
  lhs.bx.copyHorz(rhs.bx, lhsJ, rhsJ); 
  lhs.by.copyHorz(rhs.by, lhsJ, rhsJ); 
  lhs.bz.copyHorz(rhs.bz, lhsJ, rhsJ); 
                                    
  lhs.jx.copyHorz(rhs.jx, lhsJ, rhsJ); 
  lhs.jy.copyHorz(rhs.jy, lhsJ, rhsJ); 
  lhs.jz.copyHorz(rhs.jz, lhsJ, rhsJ); 
                                    
  lhs.rh.copyHorz(rhs.rh, lhsJ, rhsJ); 

}



/// Update Yee grid boundaries
// TODO: assumes implicitly 2D (x-y) arrays only by setting k=0 and then ignoring it
void maxwell::PlasmaCell::updateBoundaries(corgi::Node& node) {

  // target
  maxwell::YeeLattice& mesh = getYee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<maxwell::PlasmaCell>(
        node.getCellPtr( neighs(-1, 0) ));
  maxwell::YeeLattice& mleft = cleft->getYee();

  // copy from right side to left
  copyVertYee(mesh, mleft, -1, mleft.Nx-1); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<maxwell::PlasmaCell>(
        node.getCellPtr( neighs(+1, 0) ));
  maxwell::YeeLattice& mright = cright->getYee();
    
  // copy from left side to right
  copyVertYee(mesh, mright, mesh.Nx, 0); 


  // top 
  auto ctop = 
    std::dynamic_pointer_cast<maxwell::PlasmaCell>(
        node.getCellPtr( neighs(0, +1) ));
  maxwell::YeeLattice& mtop = ctop->getYee();

  // copy from bottom side to top
  copyHorzYee(mesh, mtop, mesh.Ny, 0); 


  // bottom
  auto cbot = 
    std::dynamic_pointer_cast<maxwell::PlasmaCell>(
        node.getCellPtr( neighs(0, -1) ));
  maxwell::YeeLattice& mbot = cbot->getYee();
    
  // copy from top side to bottom
  copyHorzYee(mesh, mbot, -1, mbot.Ny-1); 


  // diagonals/corners
  // --------------------------------------------------  
  // TODO get corners also; for now they are not needed



}



void maxwell::PlasmaCell::cycleYee() {
  yee.cycle();
}




