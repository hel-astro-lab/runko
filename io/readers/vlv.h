#pragma once

#include "io/namer.h"
#include "vlv/tile.h"

/// Read Vlasov tile content from hdf5 to Tile
template<size_t D>
bool 
h5io::Reader::read( 
  vlv::Tile<D>& tile,
  ezh5::File& file
  )
{
  //--------------------------------------------------
  // open file
  //int frank = locate_rank(tile); // find what rank had this tile

  //std::string fname = 
  //  folder + "/meshes-" + to_string(frank) + "_" + to_string(lap) + ".h5";
  ////std::cout << "reading file " << fname << std::endl;
  //  
  //ezh5::File file(fname, H5F_ACC_RDONLY);

  //--------------------------------------------------
  // open group

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // group for tile
  auto gr1 = file["tile-"+numbering];

  // get reference to the current time step 
  auto& step0 = tile.steps.get(0);

  // loop over different particle species 
  int ispc = 0; // ith particle species
  for(auto&& block0 : step0) {

    auto Nx = int(block0.Nx),
         Ny = int(block0.Ny),
         Nz = int(block0.Nz);

    // loop over different points in the tile
    for(int s=0; s<Nz; s++) {
      for(int r=0; r<Ny; r++) {
        for(int q=0; q<Nx; q++) {
          auto& mesh = block0.block(q,r,s);   // f_i
          mesh.clear();

          // group for location in tile
          auto gr2 = gr1["loc-" + create_numbering( q, r, s)];

          // group for species
          auto gr = gr2["sp-" + std::to_string(ispc)];

          // mesh metadata
          mesh.maximum_refinement_level << gr["maximum_refinement_level"];
          //mesh.error_cid              << gr["error_cid"];
          //mesh.error_index            << gr["error_index"];
          mesh.top_refinement_level     << gr["top_refinement_level"];

          std::vector<uint64_t> len_vec;
          len_vec << gr["length"];
          mesh.resize({{len_vec[0], len_vec[1], len_vec[2]}});

          //std::array<uint64_t, 3> length;
          //for (size_t ijk=0; ijk<3; ijk++) length[ijk] = len_vec[ijk];
          //mesh.resize(length);
          
          std::vector<float> mins_vec, maxs_vec;
          mins_vec << gr["mins"];
          maxs_vec << gr["maxs"];
          mesh.set_min({{ mins_vec[0], mins_vec[1], mins_vec[2]}} );
          mesh.set_max({{ maxs_vec[0], maxs_vec[1], maxs_vec[2]}} );

          //--------------------------------------------------
          // save actual mesh points
          //size_t len;
          //len << gr["number_of_cells"];

          std::vector<uint64_t> cids;
          std::vector<float>  vals;
          //cids.reserve( len );
          //vals.reserve( len );

          cids << gr["cids"];
          vals << gr["vals"];
          mesh.unserialize(cids, vals);

        } // q
      } // r
    } // s
    ispc++;
  } // end of species


  // overwrite also to other time steps
  // TODO: Does not take into account multi-stepping schemes.
  //       This only fixes boundary tile behavior.
  auto& step1 = tile.steps.get(1);
  step1 = step0;

  return true;
}



