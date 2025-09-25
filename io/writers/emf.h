#pragma once

#include "io/namer.h"
#include "core/emf/tile.h"


/// Write PlasmaTile content into a hdf5 data group
template<size_t D>
bool 
h5io::Writer::write( 
  const emf::Tile<D>& tile,
  ezh5::File& file
  )
{
  // TODO

  //const auto& gs = tile.get_const_grids();

  //// internal tile numbering 
  //auto my_ind = expand_indices( &tile );
  //string numbering = create_numbering(my_ind);

  //// open individual group for the data
  //auto gr = file["yee_"+ numbering];

  //// tile location inside grid
  //gr["i"] = static_cast<int>( std::get<0>(my_ind) );
  //gr["j"] = static_cast<int>( std::get<1>(my_ind) );
  //gr["k"] = static_cast<int>( std::get<2>(my_ind) );

  //// size
  //gr["Nx"] = static_cast<int>( gs.Nx );
  //gr["Ny"] = static_cast<int>( gs.Ny );
  //gr["Nz"] = static_cast<int>( gs.Nz );

  //--------------------------------------------------
  // Yee lattice quantities

  //gr["jx"] = gs.jx.serialize();
  //gr["jy"] = gs.jy.serialize();
  //gr["jz"] = gs.jz.serialize();

  //gr["ex"] = gs.ex.serialize();
  //gr["ey"] = gs.ey.serialize();
  //gr["ez"] = gs.ez.serialize();

  //gr["bx"] = gs.bx.serialize();
  //gr["by"] = gs.by.serialize();
  //gr["bz"] = gs.bz.serialize();

  //gr["rho"] = gs.rho.serialize();


  return true;
}



