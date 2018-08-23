#include "tile.h"

#include "amr_spatial_solver.h"


template<std::size_t D>
void vlv::Tile<D>::stepLocation(vlv::Grid<D>& grid)
{
  vlv::AmrSpatialLagrangianSolver<Realf> ssol;
  ssol.solve(*this, grid);
}



// TODO: separate into own communication module/header
template<>
vlv::PlasmaBlock& vlv::Tile<1>::get_external_data(
    vlv::Grid<1>& grid,
    int ispc,
    int i
    )
{
  auto neigh_index   = this->neighs(i); 
  uint64_t neigh_cid = grid.id( neigh_index );
  auto& tile_neigh = dynamic_cast<vlv::Tile<1>&>( grid.getTile(neigh_cid) );

  auto& species = tile_neigh.steps.get();

  return species[ispc];
}






//--------------------------------------------------
// explicit template instantiation

template class vlv::Tile<1>;
//template class vlv::Tile<2>;
//template class vlv::Tile<3>;
