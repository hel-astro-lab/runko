#pragma once

#include "namer.h"
#include "../vlasov/tile.h"
#include "../vlasov/amr/mesh.h"
#include "write_tags.h"


/// Write vlv::Tile 
template<size_t D>
bool
h5io::Writer::write(
  const vlv::Tile<D>& tile
  ) 
{

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // group for tile
  auto gr1 = file["tile-"+numbering];

  // get reference to the current time step 
  const auto& step0 = tile.steps.get(0);

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
          const auto& mesh = block0.block(q,r,s);   // f_i

          // group for location in tile
          auto gr2 = gr1["loc-" + create_numbering( q, r, s)];

          // group for species
          auto gr = gr2["sp-" + std::to_string(ispc)];

          //-------------------------------------------------- 
          // actual writing after all meta info

          // tile location inside node
          gr["i"]  = std::get<0>(my_ind);
          gr["j"]  = std::get<1>(my_ind);
          gr["k"]  = std::get<2>(my_ind);
          gr["sp"] = ispc;

          // mesh metadata
          gr["maximum_refinement_level"] = mesh.maximum_refinement_level;
          gr["error_cid"]   = mesh.error_cid;
          gr["error_index"] = mesh.error_index;
          gr["top_refinement_level"] = mesh.top_refinement_level;

          gr["length"] = mesh.length;
          gr["mins"] = mesh.mins;
          gr["maxs"] = mesh.maxs;

          size_t len = mesh.data.size();
          gr["number_of_cells"] = len;

          //--------------------------------------------------
          // save actual mesh points
          std::vector<uint64_t> cids;
          std::vector<Realf>    vals;
          cids.reserve( len );
          vals.reserve( len );

          // serialize into 1D containers
          for(uint64_t cid : mesh.get_cells(true) ) {
            cids.push_back(cid);
            vals.push_back(mesh.data.at(cid));
          }

          gr["cids"] = cids;
          gr["vals"] = vals;


        } // q
      } // r
    } // s
    ispc++;
  } // end of species

  return true;
}

