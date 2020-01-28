#pragma once

#include <string>

#include <omp.h>

#include "momentum-solvers/amr_momentum_solver.h"
#include "spatial-solvers/amr_spatial_solver.h"
#include "amr_analyzator.h"

namespace vlv{


inline void step_location( corgi::Grid<1>& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.get_tile_ids() ){
#pragma omp task
        {
            
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.get_tile( cid ));

          //vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
          //ssol.solve(tile, grid);
            
          tile.step_location(grid);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}

template<int V>
void initial_step( corgi::Grid<1>& grid )
{

  vlv::AmrMomentumLagrangianSolver<Realf,1,V> vsol;

  #pragma omp parallel 
  {
    #pragma omp single 
    {

      for(auto cid : grid.get_tile_ids() ){
        #pragma omp task firstprivate(vsol)
        {
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.get_tile( cid ));
          vsol.solve(tile, -0.5);
        }// end of omp task
      }

    }// end of omp single
  #pragma omp taskwait
  }// end of omp parallel

}




template<int V>
void step_velocity( corgi::Grid<1>& grid )
{
  vlv::AmrMomentumLagrangianSolver<Realf,1,V> vsol;

  #pragma omp parallel
  {
    #pragma omp single
    {

      for(auto cid : grid.get_tile_ids() ){
        #pragma omp task firstprivate(vsol)
        {
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.get_tile( cid ));
          vsol.solve(tile);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}

template<int V>
void step_velocity_with_gravity( 
    corgi::Grid<1>& grid,
    Realf g0,
    Realf Lx
    )
{

  vlv::GravityAmrMomentumLagrangianSolver<Realf,1,V> vsol(g0, Lx);

  #pragma omp parallel
  {
    #pragma omp single
    {

      for(auto cid : grid.get_tile_ids() ){
        #pragma omp task firstprivate(vsol)
        {
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.get_tile( cid ));
          vsol.solve(tile);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}


inline void analyze( corgi::Grid<1>& grid )
{
  vlv::Analyzator<Realf> analyzator;

  #pragma omp parallel
  {
    #pragma omp single
    {

      for(auto cid : grid.get_local_tiles()) {
        #pragma omp task firstprivate(analyzator)
        {
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.get_tile( cid ));
          analyzator.analyze(tile);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}


}// end of namespace vlv



