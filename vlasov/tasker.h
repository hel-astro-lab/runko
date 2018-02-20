#pragma once

#include <omp.h>

#include "amr_momentum_solver.h"
#include "amr_spatial_solver.h"


namespace vlasov{


void stepLocation( vlasov::Grid& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getCellIds() ){
#pragma omp task
        {
          vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
          vlasov::VlasovCell& cell 
            = dynamic_cast<vlasov::VlasovCell&>(grid.getCell( cid ));
          ssol.solve(cell, grid);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}


template<int D>
void stepVelocity( vlasov::Grid& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getCellIds() ){
#pragma omp task
        {
          vlasov::AmrMomentumLagrangianSolver<Realf,D> vsol;
          vlasov::VlasovCell& cell 
            = dynamic_cast<vlasov::VlasovCell&>(grid.getCell( cid ));
          vsol.solve(cell);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}


/// Update Yee lattice boundaries
/*
void updateBoundaries()
{

  for(auto cid : getCellIds() ){
    vlasov::VlasovCell& cell = dynamic_cast<vlasov::VlasovCell& >(getCell( cid ));
    cell.updateBoundaries( *this );
  }

}
*/




}// end of namespace
