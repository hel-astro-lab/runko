#pragma once

#include "namer.h"
#include "../pic/tile.h"
#include "write_tags.h"


template<size_t D>
bool
h5io::Writer::write(
  const pic::Tile<D>& tile
  ) 
{

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // group for tile
  auto gr1 = file["tile_"+numbering];

  // loop over different particle species 
  for(size_t ispc=0; ispc<tile.Nspecies(); ispc++) {

    // group for species + tile metainfo
    auto gr = gr1["sp-" + std::to_string(ispc)];
    gr["i"]  = static_cast<int>( std::get<0>(my_ind) );
    gr["j"]  = static_cast<int>( std::get<1>(my_ind) );
    gr["k"]  = static_cast<int>( std::get<2>(my_ind) );
    gr["sp"] = static_cast<int>( ispc );

    const pic::ParticleContainer& container = tile.get_const_container(ispc);

    // initialize pointers to particle arrays
    //int nparts = container.size();
    //double* loc[3], * vel[3];
    //for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );
    //for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

    //// loop and check particles
    //int n1 = 0;
    //int n2 = nparts;

    //for(int n=n1; n<n2; n++) {
    //  // grid coordinate location
    //  x0 = loc[0][n];
    //  y0 = loc[1][n];
    //  z0 = loc[2][n];
    //}

    gr["x"]   = container.loc(0);
    gr["y"]   = container.loc(1);
    gr["z"]   = container.loc(2);
    gr["vx"]  = container.vel(0);
    gr["vy"]  = container.vel(1);
    gr["vz"]  = container.vel(2);
    gr["wgt"] = container.wgt();

  } // end of species

  return true;
}

