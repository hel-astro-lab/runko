#pragma once

#include <string>

#include <omp.h>

#include "momentum-solvers/amr_momentum_solver.h"
#include "spatial-solvers/amr_spatial_solver.h"
#include "amr_analyzator.h"

#include "../io/writer.h"
#include "../io/reader.h"
#include "../io/write_tags.h"


namespace vlv{


inline void stepLocation( corgi::Node<1>& grid )
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
            
          tile.stepLocation(grid);
        }// end of omp task
      }


    }// end of omp single
  }// end of omp parallel

}

template<int V>
void stepInitial( corgi::Node<1>& grid )
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
void stepVelocity( corgi::Node<1>& grid )
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
void stepVelocityGravity( 
    corgi::Node<1>& grid,
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



/// Update Yee lattice boundaries
/*
void update_boundaries()
{

  for(auto cid : get_tile_ids() ){
    vlasov::VlasovTile& tile = dynamic_cast<vlasov::VlasovTile& >(get_tile( cid ));
    tile.update_boundaries( *this );
  }

}
*/


inline void analyze( corgi::Node<1>& grid )
{
  vlv::Analyzator<Realf> analyzator;

  #pragma omp parallel
  {
    #pragma omp single
    {

      for(auto cid : grid.get_tile_ids() ){
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



template<size_t D>
inline void writeYee( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir
    )
{

  std::string prefix = dir + "fields-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  for(auto cid : grid.get_tile_ids() ){
    const auto& tile 
      = dynamic_cast<fields::Tile<D>&>(grid.get_tile( cid ));
    writer.write(tile);
  }

  //writer.~Writer();
}


template<size_t D>
inline void writeAnalysis( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir
    )
{

  std::string prefix = dir + "analysis-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  for(auto cid : grid.get_tile_ids() ){
    const auto& tile 
      = dynamic_cast<fields::Tile<D>&>(grid.get_tile( cid ));
    writer.write2(tile);
  }


  //writer.~Writer();
}


template<size_t D>
inline void writeMesh( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir 
    )
{

  std::string prefix = dir + "meshes-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  for(auto cid : grid.get_tile_ids() ){
    const auto& tile 
      = dynamic_cast<vlv::Tile<D>&>(grid.get_tile( cid ));
    writer.write(tile);
  }

  //writer.~Writer();
}


template<size_t D>
inline void readYee( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir 
    )
{

  h5io::Reader reader(dir, lap);

  for(auto cid : grid.get_tile_ids() ){
    auto& tile 
      = dynamic_cast<fields::Tile<D>&>(grid.get_tile( cid ));
    reader.read(tile);
  }

}


template<size_t D>
inline void readMesh( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir 
    )
{
  h5io::Reader reader(dir, lap);

  for(auto cid : grid.get_tile_ids() ){
    auto& tile 
      = dynamic_cast<vlv::Tile<D>&>(grid.get_tile( cid ));
    reader.read(tile);
  }
}




}// end of namespace vlv
