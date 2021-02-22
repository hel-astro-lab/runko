#pragma once

#include <cmath>

#include "../namer.h"
#include "../../pic/tile.h"


template<size_t D>
bool 
h5io::Reader::read( 
  pic::Tile<D>& tile,
  ezh5::File& file
  )
{
  //--------------------------------------------------
  // open file
  //int frank = locate_rank(tile); // find what rank had this tile

  //std::string fname = 
  //  folder + "/particles-" + to_string(frank) + "_" + to_string(lap) + ".h5";
  ////std::cout << "reading file " << fname << std::endl;
  //  
  //ezh5::File file(fname, H5F_ACC_RDONLY);
    
  // tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  //maximum acceptable prtcl 4-velocity
  const real_prtcl max_vel = 1.0e8; 

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // open individual group for the tile
  auto gr1 = file["tile_"+ numbering];

  int i,j,k,sp;
  i  << gr1["i"]; 
  j  << gr1["j"]; 
  k  << gr1["k"]; 

  // loop over species and add (individually) all particles
  for (int ispc=0; ispc<tile.Nspecies(); ispc++) {
    pic::ParticleContainer<D>& container = tile.get_container(ispc);

    auto gr = gr1["sp-" + std::to_string(ispc)];
    sp << gr["sp"]; 

    // read into explicitly initialized arrays; otherwise, some -OX option
    // tries to optimize these away and we don't read anything.
    std::vector<real_prtcl> arr1, arr2, arr3, arr4, arr5, arr6, arr7;
    arr1  << gr["x"];
    arr2  << gr["y"];
    arr3  << gr["z"];
    arr4  << gr["vx"];
    arr5  << gr["vy"];
    arr6  << gr["vz"];
    arr7  << gr["wgt"];

    std::vector<int> iarr1, iarr2;
    iarr1  << gr["id"];
    iarr2  << gr["proc"];


    size_t nparts = arr1.size();

    // XXX: are these asserts needed?
    assert(arr1.size() == nparts);
    assert(arr2.size() == nparts);
    assert(arr3.size() == nparts);
    assert(arr4.size() == nparts);
    assert(arr5.size() == nparts);
    assert(arr6.size() == nparts);
    assert(arr7.size() == nparts);

    assert(iarr1.size() == nparts);
    assert(iarr2.size() == nparts);

    real_prtcl xloc, yloc, zloc;
    real_prtcl ux, uy, uz;
    for(size_t n=0; n<nparts; n++) {
      // generates key
      //container.add_particle(
      //  {arr1[n], arr2[n], arr3[n]}, 
      //  {arr4[n], arr5[n], arr6[n]}, arr7[n]);

      //--------------------------------------------------
      // check validity of read points

      xloc = arr1[n];  //gr["x"];
      yloc = arr2[n];  //gr["y"];
      zloc = arr3[n];  //gr["z"];

      ux   = arr4[n];  //gr["vx"];
      uy   = arr5[n];  //gr["vy"];
      uz   = arr6[n];  //gr["vz"];

      // confirm that these are reasonable values
      size_t nan_flag = 0;
      if(!std::isnan(xloc)) nan_flag++;
      if(!std::isnan(yloc)) nan_flag++;
      if(!std::isnan(zloc)) nan_flag++;
      if(!std::isnan(ux)  ) nan_flag++;
      if(!std::isnan(uy)  ) nan_flag++;
      if(!std::isnan(uz)  ) nan_flag++;

      // confirm that we are inside the tile
      size_t loc_flag = 0;
      if(D>= 1 && mins[0] <= xloc && xloc <= maxs[0]) loc_flag++;
      if(D>= 2 && mins[1] <= yloc && yloc <= maxs[1]) loc_flag++;
      if(D>= 3 && mins[2] <= zloc && zloc <= maxs[2]) loc_flag++;

      // confirm that velocities are reasonable
      //u = ux*ux + uy*uy + uz*uz;
      size_t vel_flag = 0;
      if( (ux > -max_vel) && (ux < max_vel) ) vel_flag++;
      if( (uy > -max_vel) && (uy < max_vel) ) vel_flag++;
      if( (uz > -max_vel) && (uz < max_vel) ) vel_flag++;

      //finally, check that all checks pass
      bool good_prtcl = true;
      if( nan_flag  < 6) good_prtcl = false;
      if( loc_flag  < D) good_prtcl = false;
      if( vel_flag  < 3) good_prtcl = false;

      if(good_prtcl) {

        // assumes old key
        container.add_identified_particle(
          {arr1[n], arr2[n], arr3[n]}, 
          {arr4[n], arr5[n], arr6[n]}, arr7[n],
          iarr1[n], iarr2[n]
          );

      } else {
        std::cerr << "skipping prtcl"
                  << " n=" << n
                  << " x=" << xloc
                  << " y=" << yloc
                  << " z=" << zloc
                  << " ux=" << ux
                  << " uy=" << uy
                  << " uz=" << uz
                  << " nan_flag=" << nan_flag
                  << " loc_flag=" << loc_flag
                  << " vel_flag=" << vel_flag
                  << std::endl;
        //assert(false);
      }
    }
  }

  return true;
}







