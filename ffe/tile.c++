#include <iostream>

#include "tile.h"


namespace ffe {
  using namespace mpi4cpp;


template<std::size_t D>
void Tile<D>::rk3_update(
    real_short c1, 
    real_short c2, 
    real_short c3
    )
{
  fields::YeeLattice&    m  = this->get_yee();
  ffe::SkinnyYeeLattice& n  = this->Fn; 
  ffe::SkinnyYeeLattice& dm = this->dF; 
  //real_short dt = tile.cfl;

  for(int k=0; k<static_cast<int>(this->mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(this->mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(this->mesh_lengths[0]); i++) {

        // RK3 E update
        m.ex(i,j,k) = c1*n.ex(i,j,k) + c2*m.ex(i,j,k) + c3*dm.ex(i,j,k);
        m.ey(i,j,k) = c1*n.ey(i,j,k) + c2*m.ey(i,j,k) + c3*dm.ey(i,j,k);
        m.ez(i,j,k) = c1*n.ez(i,j,k) + c2*m.ez(i,j,k) + c3*dm.ez(i,j,k);

        // RK3 B update
        m.bx(i,j,k) = c1*n.bx(i,j,k) + c2*m.bx(i,j,k) + c3*dm.bx(i,j,k);
        m.by(i,j,k) = c1*n.by(i,j,k) + c2*m.by(i,j,k) + c3*dm.by(i,j,k);
        m.bz(i,j,k) = c1*n.bz(i,j,k) + c2*m.bz(i,j,k) + c3*dm.bz(i,j,k);

      } 
    } 
  }

}


template<std::size_t D>
void Tile<D>::copy_eb()
{
  fields::YeeLattice&    m = this->get_yee();
  ffe::SkinnyYeeLattice& n = this->Fn; 

  for(int k=0; k<static_cast<int>(this->mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(this->mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(this->mesh_lengths[0]); i++) {

        n.ex(i,j,k) = m.ex(i,j,k);
        n.ey(i,j,k) = m.ey(i,j,k);
        n.ez(i,j,k) = m.ez(i,j,k);

        n.bx(i,j,k) = m.bx(i,j,k);
        n.by(i,j,k) = m.by(i,j,k);
        n.bz(i,j,k) = m.bz(i,j,k);

      }
    }
  }

}



//--------------------------------------------------
// explicit template instantiation
//template class Tile<2>;
template class Tile<3>;

} // end of ns ffe
