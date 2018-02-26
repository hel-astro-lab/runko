#include "cell.h"

#include "amr_spatial_solver.h"


void vlasov::VlasovCell::stepLocation(vlasov::Grid& grid)
{
  vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
  ssol.solve(*this, grid);
}










