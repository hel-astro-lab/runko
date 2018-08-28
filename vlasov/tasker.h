#pragma once

#include <string>

#include <omp.h>

#include "momentum-solvers/amr_momentum_solver.h"
#include "spatial-solvers/amr_spatial_solver.h"
#include "amr_analyzator.h"
#include "../io/io.h"


namespace vlv{


inline void stepLocation( corgi::Node<1>& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getTileIds() ){
#pragma omp task
        {
            
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.getTile( cid ));

          //vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
          //ssol.solve(tile, grid);
            
          tile.stepLocation(grid);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}

template<int V>
void stepInitial( corgi::Node<1>& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getTileIds() ){
#pragma omp task
        {
          vlv::AmrMomentumLagrangianSolver<Realf,1,V> vsol;
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.getTile( cid ));
          vsol.solve(tile, -0.5);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}




template<int V>
void stepVelocity( corgi::Node<1>& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getTileIds() ){
#pragma omp task
        {
          vlv::AmrMomentumLagrangianSolver<Realf,1,V> vsol;
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.getTile( cid ));
          vsol.solve(tile);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}

template<int V>
void stepVelocityGravity( 
    corgi::Node<1>& grid,
    Realf g0,
    Realf Lx
    )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getTileIds() ){
#pragma omp task
        {
          vlv::GravityAmrMomentumLagrangianSolver<Realf,1,V> vsol(g0, Lx);
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.getTile( cid ));
          vsol.solve(tile);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}



/// Update Yee lattice boundaries
/*
void updateBoundaries()
{

  for(auto cid : getTileIds() ){
    vlasov::VlasovTile& tile = dynamic_cast<vlasov::VlasovTile& >(getTile( cid ));
    tile.updateBoundaries( *this );
  }

}
*/


inline void analyze( corgi::Node<1>& grid )
{

#pragma omp parallel
  {
#pragma omp single
    {

      for(auto cid : grid.getTileIds() ){
#pragma omp task
        {
          vlv::Analyzator<Realf> analyzator;
          auto& tile 
            = dynamic_cast<vlv::Tile<1>&>(grid.getTile( cid ));
          analyzator.analyze(tile);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}



template<size_t D>
inline void writeYee( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir
    )
{

  std::string prefix = dir + "fields-"; 
  prefix += std::to_string(grid.rank);
  h5io::Writer writer(prefix, lap);

  for(auto cid : grid.getTileIds() ){
    auto& tile 
      = dynamic_cast<fields::Tile<D>&>(grid.getTile( cid ));
    writer.writeYee(tile);
  }


}


template<size_t D>
inline void writeAnalysis( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir
    )
{

  std::string prefix = dir + "analysis-"; 
  prefix += std::to_string(grid.rank);
  h5io::Writer writer(prefix, lap);

  for(auto cid : grid.getTileIds() ){
    auto& tile 
      = dynamic_cast<fields::Tile<D>&>(grid.getTile( cid ));
    writer.writeAnalysis(tile);
  }


}


inline void writeMesh( 
    corgi::Node<1>& grid, 
    int lap,
    const std::string& dir 
    )
{

  std::string prefix = dir + "meshes-"; 
  prefix += std::to_string(grid.rank);
  h5io::Writer writer(prefix, lap);


  for(auto cid : grid.getTileIds() ){

    auto& tile 
      = dynamic_cast<vlv::Tile<1>&>(grid.getTile( cid ));


    // tile location index
    int i = std::get<0>(tile.index);
    int j = 0;
    int k = 0;

    // get reference to the current time step 
    auto& step0 = tile.steps.get(0);

    // loop over different particle species 
    int ispc = 0; // ith particle species
    for(auto&& block0 : step0) {

      auto Nx = int(block0.Nx),
           Ny = int(block0.Ny),
           Nz = int(block0.Nz);

      for(int s=0; s<Nz; s++) {
        for(int r=0; r<Ny; r++) {
          for(int q=0; q<Nx; q++) {
            const auto& M   = block0.block(q,r,s);   // f_i

            // information about location is encoded in:
            // i j k | q r s | ispc
            //
            // that is formatted into:
            // tile-i_j_k/loc-q_r_s/sp-ispc

            h5io::TileInfo tinfo;
            tinfo.prefix = "tile";

            tinfo.i = i;
            tinfo.j = j;
            tinfo.k = k;

            tinfo.q = q;
            tinfo.r = r;
            tinfo.s = s;

            tinfo.sp = ispc;

              
            writer.write(M, tinfo);

          } // q
        } // r
      } // s
      ispc++;
    } // end of species
  } // end of loop over tiles

}


}// end of namespace vlv
