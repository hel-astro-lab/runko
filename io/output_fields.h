#pragma once

#include "namer.h"
#include "../em-fields/tile.h"
#include "write_tags.h"


/// Write PlasmaTile content into a hdf5 data group
template<size_t D>
bool 
h5io::Writer::write( 
  const fields::Tile<D>& tile
  )
{
  const auto& yee = tile.get_const_yee();

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // open individual group for the data
  auto gr = file["yee_"+ numbering];

  // tile location inside grid
  gr["i"] = static_cast<int>( std::get<0>(my_ind) );
  gr["j"] = static_cast<int>( std::get<1>(my_ind) );
  gr["k"] = static_cast<int>( std::get<2>(my_ind) );

  // size
  gr["Nx"] = static_cast<int>( yee.Nx );
  gr["Ny"] = static_cast<int>( yee.Ny );
  gr["Nz"] = static_cast<int>( yee.Nz );

  //--------------------------------------------------
  // Yee lattice quantities

  gr["jx"] = yee.jx.serialize();
  gr["jy"] = yee.jy.serialize();
  gr["jz"] = yee.jz.serialize();

  gr["ex"] = yee.ex.serialize();
  gr["ey"] = yee.ey.serialize();
  gr["ez"] = yee.ez.serialize();

  gr["bx"] = yee.bx.serialize();
  gr["by"] = yee.by.serialize();
  gr["bz"] = yee.bz.serialize();

  gr["rho"] = yee.rho.serialize();


  return true;
}



