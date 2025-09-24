#include "core/vlv/tile.h"
#include "core/vlv/spatial-solvers/amr_spatial_solver.h"

namespace vlv {

template<std::size_t D>
void Tile<D>::step_location(corgi::Grid<D>& grid)
{
  AmrSpatialLagrangianSolver<float> ssol;
  ssol.solve(*this, grid);
}



// TODO: separate into own communication module/header
template<>
PlasmaBlock& Tile<1>::get_external_data(
    corgi::Grid<1>& grid,
    int ispc,
    int i
    )
{
  auto neigh_index   = this->neighs(i); 
  uint64_t neigh_cid = grid.id( neigh_index );
  auto& tile_neigh = dynamic_cast<Tile<1>&>( grid.get_tile(neigh_cid) );

  auto& species = tile_neigh.steps.get();

  return species[ispc];
}


/// cycle temporary and true current arrays
template<std::size_t D>
void Tile<D>::cycle_current() 
{
  //auto& gs = this->get_grids();
  // fix these accessing mat 
  // std::swap( gs.jx.mat, jx1.mat );
  // std::swap( gs.jy.mat, jy1.mat );
  // std::swap( gs.jz.mat, jz1.mat );
}




//--------------------------------------------------
// explicit template instantiation

template class Tile<1>;
//template class Tile<2>;
//template class Tile<3>;

} // end of ns vlv
