#pragma once

#include "io/namer.h"
#include "core/pic/tile.h"

template<size_t D>
bool
h5io::Writer::write(
  const pic::Tile<D>& tile,
  ezh5::File& file
  ) 
{
  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // group for tile
  auto gr1 = file["tile_"+numbering];

  gr1["i"] = static_cast<int>( std::get<0>(my_ind) );
  gr1["j"] = static_cast<int>( std::get<1>(my_ind) );
  gr1["k"] = static_cast<int>( std::get<2>(my_ind) );

  // TODO
  // loop over different particle species 
  //for(int ispc=0; ispc<tile.Nspecies(); ispc++) {

  //  // group for species + tile metainfo
  //  auto gr = gr1["sp-" + std::to_string(ispc)];
  //  gr["sp"] = static_cast<int>( ispc );

  //  // todo set back to const ref at some point...
  //  pic::ParticleContainer<D> container = tile.get_const_container(ispc);

  //  gr["x"]   = container.loc(0);
  //  gr["y"]   = container.loc(1);
  //  gr["z"]   = container.loc(2);

  //  gr["vx"]  = container.vel(0);
  //  gr["vy"]  = container.vel(1);
  //  gr["vz"]  = container.vel(2);

  //  gr["wgt"] = container.wgt();

  //  gr["id"]  = container.id(0);
  //  gr["proc"]= container.id(1);

  //} // end of species

  return true;
}

