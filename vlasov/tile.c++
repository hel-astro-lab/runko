#include "tile.h"

#include "amr_spatial_solver.h"


void vlasov::VlasovTile::stepLocation(vlasov::Grid& grid)
{
  vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
  ssol.solve(*this, grid);
}


// TODO: separate into own communication module/header
vlasov::PlasmaBlock& vlasov::VlasovTile::get_external_data(
    int i, int j, int ispc,
    vlasov::Grid& grid)
{
  auto neigh_index   = neighs(i, j); 
  uint64_t neigh_cid = grid.tileId( std::get<0>(neigh_index), std::get<1>(neigh_index) );
  auto& tile_neigh = dynamic_cast<vlasov::VlasovTile&>( grid.getTile(neigh_cid) );

  auto& species = tile_neigh.steps.get();

  return species[ispc];
}








