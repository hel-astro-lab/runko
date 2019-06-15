#include <iostream>

#include "tile.h"


namespace ffe {
  using namespace mpi4cpp;


/// Get current integration time snapshot of Yee lattice
template<std::size_t D>
fields::YeeLattice& Tile<D>::get_yee(size_t i) 
{
  std::cout << "getting yee from snapshotter \n";
  return this->yee_snapshots.get(i);
}

//--------------------------------------------------
// Specialize Yee Lattice insertion
template<>
void Tile<1>::add_yee_lattice() 
{
  //std::cout << "add 1D Yee \n";
  yee_snapshots.push_back( 
      fields::YeeLattice( mesh_lengths[0], 1, 1) 
      );
}

template<>
void Tile<2>::add_yee_lattice() 
{
  yee_snapshots.push_back( 
      fields::YeeLattice( mesh_lengths[0], mesh_lengths[1], 1) 
      );
}

template<>
void Tile<3>::add_yee_lattice() 
{
  yee_snapshots.push_back( 
      fields::YeeLattice( mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]) 
      );
}



/// 2D Force-free current calculator
template<>
void Tile<2>::compute_perp_current() 
{
  fields::YeeLattice& mesh = get_yee();

  double C = 0.5 * cfl;

  int k = 0;
  for(int j=0; j<static_cast<int>(mesh_lengths[1]); j++)
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Bx
    mesh.bx(i,j,k) += 
      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)
          -mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

  }

}


/// 3D Force-free current calculator
template<>
void Tile<3>::compute_perp_current() 
{
  fields::YeeLattice& mesh = get_yee();
  double C = 0.5 * cfl;

  for(int k=0; k<static_cast<int>(mesh_lengths[2]); k++) 
  for(int j=0; j<static_cast<int>(mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Bx
    mesh.bx(i,j,k) += 
      + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k))
      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k))
      + C*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k))
      + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

  }
}

/// 2D Force-free E.B = 0 enforcer
template<>
void Tile<2>::subtract_parallel_e() 
{
  //fields::YeeLattice& mesh = get_yee();
}



//--------------------------------------------------
// explicit template instantiation
template class Tile<2>;
//template class Tile<3>;

} // end of ns ffe
