#pragma once

#include <omp.h>

#include "amr_momentum_solver.h"
#include "amr_spatial_solver.h"
#include "amr_analyzator.h"


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
            
          vlasov::VlasovCell& cell 
            = dynamic_cast<vlasov::VlasovCell&>(grid.getCell( cid ));

          //vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
          //ssol.solve(cell, grid);
            
          cell.stepLocation(grid);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}

template<int D>
void stepInitial( vlasov::Grid& grid )
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
          vsol.solve(cell, -0.5);
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


void analyze( vlasov::Grid& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getCellIds() ){
#pragma omp task
        {
          vlasov::Analyzator<Realf> analyzator;
          vlasov::VlasovCell& cell 
            = dynamic_cast<vlasov::VlasovCell&>(grid.getCell( cid ));
          analyzator.analyze(cell);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}


}// end of namespace
