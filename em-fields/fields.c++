#include <iostream>

#include "fields.h"


fields::PlasmaCell::PlasmaCell(
    size_t i, size_t j,
    int o,
    size_t NxG, size_t NyG,
    size_t NxMesh, size_t NyMesh, size_t NzMesh
    ) :
  corgi::Cell(i, j, o, NxG, NyG),
  NxMesh(NxMesh), NyMesh(NyMesh), NzMesh(NzMesh)
{

  // append fresh Yee lattices into data container
  yee.push_back( fields::YeeLattice(NxMesh, NyMesh, NzMesh) );
  // yee.push_back( fields::YeeLattice(NxMesh, NyMesh, NzMesh) );



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
void fields::PlasmaCell::pushE() {

  // this->pushE1d();
  this->pushE2d();
  // this->pushE3d();
}

/// 1D E pusher
void fields::PlasmaCell::pushE1d() {

  fields::YeeLattice& mesh = getYee();

  int k = 0;
  int j = 0;
  for(int i=0; i<(int)NxMesh; i++) {

    // Ex
    // NONE

    // Ey
    mesh.ey(i,j,k) += 
      + yeeDt*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k)) / yeeDx;

    // Ez
    mesh.ez(i,j,k) += 
      + yeeDt*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k)) / yeeDx;

  }

}

/// 2D E pusher
void fields::PlasmaCell::pushE2d() {

  fields::YeeLattice& mesh = getYee();


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


/// 3D E pusher
void fields::PlasmaCell::pushE3d() {

  fields::YeeLattice& mesh = getYee();

  for(int k=0; k<(int)NzMesh; k++) {
    for(int j=0; j<(int)NyMesh; j++) {
      for(int i=0; i<(int)NxMesh; i++) {

        // Ex
        mesh.ex(i,j,k) += 
          + yeeDt*( mesh.by(i,j,  k-1) - mesh.by(i,j,k)) / yeeDx
          + yeeDt*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k)) / yeeDx;

        // Ey
        mesh.ey(i,j,k) += 
          + yeeDt*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k)) / yeeDx
          + yeeDt*(-mesh.bx(i,  j, k-1) + mesh.bx(i,j,k)) / yeeDx;

        // Ez
        mesh.ez(i,j,k) += 
          + yeeDt*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k)) / yeeDx
          + yeeDt*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k)) / yeeDx;

      }
    }
  }

}




/// Deposit current into electric field
void fields::PlasmaCell::depositCurrent() {
  fields::YeeLattice& mesh = getYee();

  //std::cout<<"dt:"<<yeeDt<<"  and vol:"<<yeeDx<<" .. " <<(yeeDx*yeeDx*yeeDx) <<"\n";

  mesh.ex -= mesh.jx;
  mesh.ey -= mesh.jy;
  mesh.ez -= mesh.jz;

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
void fields::PlasmaCell::pushHalfB() {

  // this->pushHalfB1d();
  this->pushHalfB2d();
  // this->pushHalfB3d();
}

/// 1D B pusher
void fields::PlasmaCell::pushHalfB1d() {
  fields::YeeLattice& mesh = getYee();

  int k = 0;
  int j = 0;
  for(int i=0; i<(int)NxMesh; i++) {

    // Bx
    // NONE

    // By
    mesh.by(i,j,k) += 
      + yeeDt*0.5*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) / yeeDx;

    // Bz
    mesh.bz(i,j,k) += 
      + yeeDt*0.5*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k)) / yeeDx;
  }

}

/// 2D B pusher
void fields::PlasmaCell::pushHalfB2d() {
  fields::YeeLattice& mesh = getYee();

  int k = 0;
  for(int j=0; j<(int)NyMesh; j++) {
    for(int i=0; i<(int)NxMesh; i++) {

      // Bx
      mesh.bx(i,j,k) += 
        + yeeDt*0.5*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k)) / yeeDx;

      // By
      mesh.by(i,j,k) += 
        + yeeDt*0.5*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) / yeeDx;

      // Bz
      mesh.bz(i,j,k) += 
        + yeeDt*0.5*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)) / yeeDx
        + yeeDt*0.5*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k)) / yeeDx;

    }
  }

}


/// 3D B pusher
void fields::PlasmaCell::pushHalfB3d() {
  fields::YeeLattice& mesh = getYee();

  for(int k=0; k<(int)NzMesh; k++) {
    for(int j=0; j<(int)NyMesh; j++) {
      for(int i=0; i<(int)NxMesh; i++) {

        // Bx
        mesh.bx(i,j,k) += 
         + yeeDt*0.5*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k)) / yeeDx
         + yeeDt*0.5*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k)) / yeeDx;

        // By
        mesh.by(i,j,k) += 
         + yeeDt*0.5*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) / yeeDx
         + yeeDt*0.5*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k)) / yeeDx;

        // Bz
        mesh.bz(i,j,k) += 
         + yeeDt*0.5*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)) / yeeDx
         + yeeDt*0.5*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k)) / yeeDx;

      }
    }
  }

}




/// Get current time snapshot of Yee lattice
fields::YeeLattice& fields::PlasmaCell::getYee(size_t i) {
  return yee.get(i);
}


/// Get analysis lattice of i:th species
fields::PlasmaMomentLattice& fields::PlasmaCell::getAnalysis(size_t i) {
  return analysis[i];
}

void fields::PlasmaCell::addAnalysisSpecies() {
  analysis.push_back( PlasmaMomentLattice(NxMesh, NyMesh, NzMesh) );
}


/// Quick helper function to copy everything inside Yee lattice 
void copyVertYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsI, int rhsI) {

  lhs.ex.copyVert(rhs.ex, lhsI, rhsI); 
  lhs.ey.copyVert(rhs.ey, lhsI, rhsI); 
  lhs.ez.copyVert(rhs.ez, lhsI, rhsI); 

  lhs.bx.copyVert(rhs.bx, lhsI, rhsI); 
  lhs.by.copyVert(rhs.by, lhsI, rhsI); 
  lhs.bz.copyVert(rhs.bz, lhsI, rhsI); 

}


/// Quick helper function to add everything inside Yee lattice 
void addVertYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsI, int rhsI) {

  lhs.jx.addVert(rhs.jx, lhsI, rhsI); 
  lhs.jy.addVert(rhs.jy, lhsI, rhsI); 
  lhs.jz.addVert(rhs.jz, lhsI, rhsI); 

}



/// Quick helper function to copy everything inside Yee lattice 
void copyHorzYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsJ, int rhsJ) {

  lhs.ex.copyHorz(rhs.ex, lhsJ, rhsJ); 
  lhs.ey.copyHorz(rhs.ey, lhsJ, rhsJ); 
  lhs.ez.copyHorz(rhs.ez, lhsJ, rhsJ); 
                                    
  lhs.bx.copyHorz(rhs.bx, lhsJ, rhsJ); 
  lhs.by.copyHorz(rhs.by, lhsJ, rhsJ); 
  lhs.bz.copyHorz(rhs.bz, lhsJ, rhsJ); 

}


/// Quick helper function to add everything inside Yee lattice 
void addHorzYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsJ, int rhsJ) {
                                    
  lhs.jx.addHorz(rhs.jx, lhsJ, rhsJ); 
  lhs.jy.addHorz(rhs.jy, lhsJ, rhsJ); 
  lhs.jz.addHorz(rhs.jz, lhsJ, rhsJ); 
}

/// Quick helper function to copy everything inside Yee lattice 
void copyFaceYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsK, int rhsK) {

  lhs.ex.copyFace(rhs.ex, lhsK, rhsK); 
  lhs.ey.copyFace(rhs.ey, lhsK, rhsK); 
  lhs.ez.copyFace(rhs.ez, lhsK, rhsK); 
                                    
  lhs.bx.copyFace(rhs.bx, lhsK, rhsK); 
  lhs.by.copyFace(rhs.by, lhsK, rhsK); 
  lhs.bz.copyFace(rhs.bz, lhsK, rhsK); 

}


/// Quick helper function to add everything inside Yee lattice 
void addFaceYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsK, int rhsK) {
                                    
  lhs.jx.addFace(rhs.jx, lhsK, rhsK); 
  lhs.jy.addFace(rhs.jy, lhsK, rhsK); 
  lhs.jz.addFace(rhs.jz, lhsK, rhsK); 
}


void copyZdirPencilYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsI, int lhsJ,
    int rhsI, int rhsJ) {

  lhs.ex.copyZdirPencil(rhs.ex, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.ey.copyZdirPencil(rhs.ey, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.ez.copyZdirPencil(rhs.ez, lhsI, lhsJ, rhsI, rhsJ); 

  lhs.bx.copyZdirPencil(rhs.bx, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.by.copyZdirPencil(rhs.by, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.bz.copyZdirPencil(rhs.bz, lhsI, lhsJ, rhsI, rhsJ); 

}

void addZdirPencilYee(
    fields::YeeLattice& lhs, 
    fields::YeeLattice& rhs, 
    int lhsI, int lhsJ,
    int rhsI, int rhsJ) {

  lhs.jx.addZdirPencil(rhs.ex, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.jy.addZdirPencil(rhs.ey, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.jz.addZdirPencil(rhs.ez, lhsI, lhsJ, rhsI, rhsJ); 

}




/// Update Yee grid boundaries
void fields::PlasmaCell::updateBoundaries(corgi::Node& node) {

  // target
  fields::YeeLattice& mesh = getYee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, 0) ));
  fields::YeeLattice& mleft = cleft->getYee();

  // copy from right side to left
  copyVertYee(mesh, mleft, -1, mleft.Nx-1); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, 0) ));
  fields::YeeLattice& mright = cright->getYee();
    
  // copy from left side to right
  copyVertYee(mesh, mright, mesh.Nx, 0); 

}


/// Update Yee grid boundaries
void fields::PlasmaCell::updateBoundaries2D(corgi::Node& node) {

  // target
  fields::YeeLattice& mesh = getYee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, 0) ));
  fields::YeeLattice& mleft = cleft->getYee();

  // copy from right side to left
  copyVertYee(mesh, mleft, -1, mleft.Nx-1); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, 0) ));
  fields::YeeLattice& mright = cright->getYee();
    
  // copy from left side to right
  copyVertYee(mesh, mright, mesh.Nx, 0); 


  // TODO: fix these: they produce saw-like oscillations
  // top 
  auto ctop = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(0, +1) ));
  fields::YeeLattice& mtop = ctop->getYee();

  //copy from bottom side to top
  copyHorzYee(mesh, mtop, mesh.Ny, 0); 


  // bottom
  auto cbot = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(0, -1) ));
  fields::YeeLattice& mbot = cbot->getYee();
    
  // copy from top side to bottom
  copyHorzYee(mesh, mbot, -1, mbot.Ny-1); 


  // --------------------------------------------------  
  // diagonals
  auto ctopleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, +1) ));
  fields::YeeLattice& mtopleft = ctopleft->getYee();
  copyZdirPencilYee(mesh, mtopleft, -1, mesh.Ny,
                                 mtopleft.Nx-1, 0);

  auto ctopright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, +1) ));
  fields::YeeLattice& mtopright = ctopright->getYee();
  copyZdirPencilYee(mesh, mtopright, mesh.Nx, mesh.Ny,
                                  0,0);

  auto cbotleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, -1) ));
  fields::YeeLattice& mbotleft = cbotleft->getYee();
  copyZdirPencilYee(mesh, mbotleft, -1, -1,
                          mbotleft.Nx-1, mbotleft.Ny-1);

  auto cbotright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, -1) ));
  fields::YeeLattice& mbotright = cbotright->getYee();
  copyZdirPencilYee(mesh, mbotright, mesh.Nx, -1,
                          0, mbotright.Ny-1);


  // --------------------------------------------------  
  // front
  // TODO: hack to deal with 2D corgi tiles
  //copyFaceYee(mesh, mesh, -1, mesh.Nz-1);

  // back
  //copyFaceYee(mesh, mesh, mesh.Nz, 0);
  // --------------------------------------------------  
  // TODO: x-pencils
  // TODO: y-pencils
  // TODO: corners


}

void fields::PlasmaCell::exchangeCurrents(corgi::Node& node) {

  // target
  fields::YeeLattice& mesh = getYee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, 0) ));
  fields::YeeLattice& mleft = cleft->getYee();

  // add from right side to left
  addVertYee(mesh, mleft, -1, mleft.Nx-1); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, 0) ));
  fields::YeeLattice& mright = cright->getYee();
    
  // add from left side to right
  addVertYee(mesh, mright, mesh.Nx, 0); 

}


/// Update currents on Yee grid boundaries
// TODO: assumes implicitly 2D (x-y) arrays only by setting k=0 and then ignoring it
void fields::PlasmaCell::exchangeCurrents2D(corgi::Node& node) {

  // target
  fields::YeeLattice& mesh = getYee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, 0) ));
  fields::YeeLattice& mleft = cleft->getYee();

  // add from right side to left
  addVertYee(mesh, mleft, -1, mleft.Nx-1); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, 0) ));
  fields::YeeLattice& mright = cright->getYee();
    
  // add from left side to right
  addVertYee(mesh, mright, mesh.Nx, 0); 

  // top 
  auto ctop = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(0, +1) ));
  fields::YeeLattice& mtop = ctop->getYee();

  //copy from bottom side to top
  addHorzYee(mesh, mtop, mesh.Ny, 0); 


  // bottom
  auto cbot = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(0, -1) ));
  fields::YeeLattice& mbot = cbot->getYee();
    
  // copy from top side to bottom
  addHorzYee(mesh, mbot, -1, mbot.Ny-1); 

  // --------------------------------------------------  
  // diagonals
  auto ctopleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, +1) ));
  fields::YeeLattice& mtopleft = ctopleft->getYee();
  addZdirPencilYee(mesh, mtopleft, -1, mesh.Ny,
                                 mtopleft.Nx-1, 0);

  auto ctopright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, +1) ));
  fields::YeeLattice& mtopright = ctopright->getYee();
  addZdirPencilYee(mesh, mtopright, mesh.Nx, mesh.Ny,
                                  0,0);

  auto cbotleft = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(-1, -1) ));
  fields::YeeLattice& mbotleft = cbotleft->getYee();
  addZdirPencilYee(mesh, mbotleft, -1, -1,
                                  mbotleft.Nx-1, mbotleft.Ny-1);

  auto cbotright = 
    std::dynamic_pointer_cast<fields::PlasmaCell>(
        node.getCellPtr( neighs(+1, -1) ));
  fields::YeeLattice& mbotright = cbotright->getYee();
  addZdirPencilYee(mesh, mbotright, mesh.Nx, -1,
                                  0, mbotright.Ny-1);


  // front
  // TODO: hack to deal with 2D corgi tiles
  //addFaceYee(mesh, mesh, -1, mesh.Nz-1);

  // back
  //addFaceYee(mesh, mesh, mesh.Nz, 0);



}

void fields::PlasmaCell::cycleYee() {
  yee.cycle();
}




