#include <iostream>

#include "maxwell.h"
#include "tools/Mesh.h"


maxwell::PlasmaCell::PlasmaCell(
    size_t i, size_t j,
    int o,
    size_t Nx, size_t Ny) :
  corgi::Cell(i, j, o, Nx, Ny),
  Nx(Nx), Ny(Ny), Nz(1)
{

  // append fresh Yee lattices into data container
  yee.push_back( maxwell::YeeLattice(Nx, Ny, Nz) );
  yee.push_back( maxwell::YeeLattice(Nx, Ny, Nz) );


}

/// Update E field with full step
void maxwell::PlasmaCell::pushE() {
  std::cout << "push E\n";
}

/// Update B field with a half step
void maxwell::PlasmaCell::pushHalfB() {
  std::cout << "push B\n";
}


/// Get current time snapshot of Yee lattice
maxwell::YeeLattice& maxwell::PlasmaCell::getYee() {
  return *yee.get();
}

/// Get new time snapshot of Yee lattice
maxwell::YeeLattice& maxwell::PlasmaCell::getNewYee() {
  return *yee.getNew();
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

  std::cout << "updating boundaries...\n";

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






