#pragma once

#include "namer.h"
#include "../em-fields/tile.h"

/// Read PlasmaTile content from hdf5 to Tile
template<size_t D>
bool 
h5io::Reader::read( 
  fields::Tile<D>& tile
  )
{
  //--------------------------------------------------
  // open file
  int frank = locate_rank(tile); // find what rank had this tile

  std::string fname = 
    folder + "/fields-" + to_string(frank) + "_" + to_string(lap) + ".h5";
  //std::cout << "reading file " << fname << std::endl;
    
  ezh5::File file(fname, H5F_ACC_RDONLY);

  //--------------------------------------------------
  // open group

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // open individual group for the data
  auto gr = file["yee_"+ numbering];


  //--------------------------------------------------
  // read

  // container that we update
  auto& yee = tile.getYee();

  // tile location inside node
  int i,j,k;
  i << gr["i"]; 
  j << gr["j"]; 
  k << gr["k"]; 

  //std::cout << "tile i=" << std::get<0>(my_ind) << "vs read value: " << i << "\n";
  //std::cout << "tile j=" << std::get<1>(my_ind) << "vs read value: " << j << "\n";
  //std::cout << "tile k=" << std::get<2>(my_ind) << "vs read value: " << k << "\n";

  // size
  int Nx, Ny, Nz;
  Nx << gr["Nx"];
  Ny << gr["Ny"];
  Nz << gr["Nz"];

  //std::cout << "tile Nx=" << yee.Nx << "vs read value: " << Nx << "\n";
  //std::cout << "tile Ny=" << yee.Ny << "vs read value: " << Ny << "\n";
  //std::cout << "tile Nz=" << yee.Nz << "vs read value: " << Nz << "\n";

  // meshes
  std::vector<Realf> arr;

  arr.clear();
  arr << gr["jx"];
  yee.jx.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["jy"];
  yee.jy.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["jz"];
  yee.jz.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["ex"];
  yee.ex.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["ey"];
  yee.ey.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["ez"];
  yee.ez.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["bx"];
  yee.bx.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["by"];
  yee.by.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["bz"];
  yee.bz.unserialize(arr, Nx, Ny, Nz);

  arr.clear();
  arr << gr["rho"];
  yee.rho.unserialize(arr, Nx, Ny, Nz);


  // file handle is closed automatically here as it goes out-of-scope
  return true;
}






