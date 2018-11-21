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
  auto& yee = tile.get_yee();

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
  std::vector<Realf> arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, arr9, arr10;

  arr1  << gr["jx"];
  arr2  << gr["jy"];
  arr3  << gr["jz"];
  arr4  << gr["ex"];
  arr5  << gr["ey"];
  arr6  << gr["ez"];
  arr7  << gr["bx"];
  arr8  << gr["by"];
  arr9  << gr["bz"];
  arr10 << gr["rho"];

  yee.jx.unserialize(arr1,   Nx, Ny, Nz);
  yee.jy.unserialize(arr2,   Nx, Ny, Nz);
  yee.jz.unserialize(arr3,   Nx, Ny, Nz);
  yee.ex.unserialize(arr4,   Nx, Ny, Nz);
  yee.ey.unserialize(arr5,   Nx, Ny, Nz);
  yee.ez.unserialize(arr6,   Nx, Ny, Nz);
  yee.bx.unserialize(arr7,   Nx, Ny, Nz);
  yee.by.unserialize(arr8,   Nx, Ny, Nz);
  yee.bz.unserialize(arr9,   Nx, Ny, Nz);
  yee.rho.unserialize(arr10, Nx, Ny, Nz);

  // file handle is closed automatically here as it goes out-of-scope
  // TODO: seems that it does not...
  file.~File();

  return true;
}






