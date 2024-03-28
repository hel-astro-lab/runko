#include <iostream>

#include "core/ffe/tile.h"
#include "external/iter/iter.h"

namespace ffe {
  using namespace mpi4cpp;


template<std::size_t D>
void Tile<D>::rk3_update(
    float c1, 
    float c2, 
    float c3
    )
{
  emf::Grids&    m  = this->get_grids();
  ffe::SlimGrids& n  = this->Fn; 
  ffe::SlimGrids& dm = this->dF; 
  //float dt = tile.cfl;


  UniIter::iterate3D(
    [=] DEVCALLABLE(
      int i, int j, int k,
      ffe::SlimGrids& dm,
      emf::Grids& m,
      ffe::SlimGrids& n)
      {
        // RK3 E update
        m.ex(i,j,k) = c1*n.ex(i,j,k) + c2*m.ex(i,j,k) + c3*dm.ex(i,j,k);
        m.ey(i,j,k) = c1*n.ey(i,j,k) + c2*m.ey(i,j,k) + c3*dm.ey(i,j,k);
        m.ez(i,j,k) = c1*n.ez(i,j,k) + c2*m.ez(i,j,k) + c3*dm.ez(i,j,k);

        // RK3 B update
        m.bx(i,j,k) = c1*n.bx(i,j,k) + c2*m.bx(i,j,k) + c3*dm.bx(i,j,k);
        m.by(i,j,k) = c1*n.by(i,j,k) + c2*m.by(i,j,k) + c3*dm.by(i,j,k);
        m.bz(i,j,k) = c1*n.bz(i,j,k) + c2*m.bz(i,j,k) + c3*dm.bz(i,j,k);

        // variable switch for 1) e > b and 2) j_par calcs.
        // Enables to calculate both of the above as independent
        // corrections because interpolation is done via m.ex
        // meshes and results are stored in dm.ex meshes:
        dm.ex(i,j,k) = m.ex(i,j,k);
        dm.ey(i,j,k) = m.ey(i,j,k);
        dm.ez(i,j,k) = m.ez(i,j,k);
      },
    static_cast<int>(mesh_lengths[2]),
    static_cast<int>(mesh_lengths[1]),
    static_cast<int>(mesh_lengths[0]),
    dm, m, n);

#ifdef GPU
    UniIter::sync();
#endif

}


template<std::size_t D>
void Tile<D>::copy_eb()
{
  emf::Grids&    m = this->get_grids();
  ffe::SlimGrids& n = this->Fn; 

  UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, emf::Grids& m, ffe::SlimGrids& n)
    {
        n.ex(i,j,k) = m.ex(i,j,k);
        n.ey(i,j,k) = m.ey(i,j,k);
        n.ez(i,j,k) = m.ez(i,j,k);

        n.bx(i,j,k) = m.bx(i,j,k);
        n.by(i,j,k) = m.by(i,j,k);
        n.bz(i,j,k) = m.bz(i,j,k);
  },
    static_cast<int>(mesh_lengths[2]),
    static_cast<int>(mesh_lengths[1]),
    static_cast<int>(mesh_lengths[0]),
    m, n);

#ifdef GPU
//nvtxRangePop();
  UniIter::sync();
#endif

}



//--------------------------------------------------
// explicit template instantiation
//template class Tile<2>;
template class Tile<3>;

} // end of ns ffe
