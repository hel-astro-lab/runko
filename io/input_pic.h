#pragma once

#include "namer.h"
#include "../pic/tile.h"


template<size_t D>
bool 
h5io::Reader::read( 
  pic::Tile<D>& tile
  )
{
  //--------------------------------------------------
  // open file
  int frank = locate_rank(tile); // find what rank had this tile

  std::string fname = 
    folder + "/particles-" + to_string(frank) + "_" + to_string(lap) + ".h5";
  //std::cout << "reading file " << fname << std::endl;
    
  ezh5::File file(fname, H5F_ACC_RDONLY);

  //--------------------------------------------------
  // open group

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // open individual group for the data
  auto gr1 = file["tile_"+ numbering];


  // tile location inside node
  int i,j,k,sp;
  for (size_t ispc=0; ispc<tile.Nspecies(); ispc++) {

    auto gr = gr1["sp-" + std::to_string(ispc)];

    pic::ParticleContainer& container = tile.get_container(ispc);

    i  << gr["i"]; 
    j  << gr["j"]; 
    k  << gr["k"]; 
    sp << gr["sp"]; 

    std::vector<double> arr1, arr2, arr3, arr4, arr5, arr6, arr7;
    arr1  << gr["x"];
    arr2  << gr["y"];
    arr3  << gr["z"];
    arr4  << gr["vx"];
    arr5  << gr["vy"];
    arr6  << gr["vz"];
    arr7  << gr["wgt"];

    //container.x = arr1;
    //...

  }



  file.~File();

  return true;
}







