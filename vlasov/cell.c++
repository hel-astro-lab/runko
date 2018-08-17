#include "cell.h"

#include "amr_spatial_solver.h"


void vlasov::VlasovCell::stepLocation(vlasov::Grid& grid)
{
  vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
  ssol.solve(*this, grid);
}


// TODO: separate into own communication module/header
vlasov::PlasmaBlock& vlasov::VlasovCell::get_external_data(
    int i, int j, int ispc,
    vlasov::Grid& grid)
{
  auto neigh_index   = neighs(i, j); 
  uint64_t neigh_cid = grid.cellId( std::get<0>(neigh_index), std::get<1>(neigh_index) );
  auto& cell_neigh = dynamic_cast<vlasov::VlasovCell&>( grid.getCell(neigh_cid) );

  auto& species = cell_neigh.steps.get();

  return species[ispc];
}








