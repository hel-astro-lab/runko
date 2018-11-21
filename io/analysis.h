#pragma once

#include "namer.h"
#include "write_tags.h"
#include "../em-fields/tile.h"


/// Write PlasmaTile content into a hdf5 data group
template<size_t D>
bool 
h5io::Writer::write2( 
  const fields::Tile<D>& tile
    ) 
{

  // internal tile numbering 
  auto my_ind = expand_indices( &tile );
  string numbering = create_numbering(my_ind);

  // loop over all species
  int Nspecies = static_cast<int>(tile.analysis.size());
  for(int ispcs = 0; ispcs < Nspecies; ispcs++) {

    const auto& analysis = tile.get_const_analysis(ispcs);

    // open individual group for the data
    auto gr = file["analysis_"+numbering+"-"+to_string(ispcs) ];

    // tile location inside node
    gr["i"] = std::get<0>(my_ind);
    gr["j"] = std::get<1>(my_ind);
    gr["k"] = std::get<2>(my_ind);
    gr["ispcs"] = ispcs;

    // size
    gr["Nx"] = analysis.Nx;
    gr["Ny"] = analysis.Ny;
    gr["Nz"] = analysis.Nz;

    gr["rho"] = analysis.rho.serialize();
    gr["edens"] = analysis.edens.serialize();
    gr["temp"] = analysis.temp.serialize();

    gr["Vx"] = analysis.Vx.serialize();
    gr["Vy"] = analysis.Vy.serialize();
    gr["Vz"] = analysis.Vz.serialize();

    gr["momx"] = analysis.momx.serialize();
    gr["momy"] = analysis.momy.serialize();
    gr["momz"] = analysis.momz.serialize();

    gr["pressx"] = analysis.pressx.serialize();
    gr["pressy"] = analysis.pressy.serialize();
    gr["pressz"] = analysis.pressz.serialize();

    gr["shearxy"] = analysis.shearxy.serialize();
    gr["shearxz"] = analysis.shearxz.serialize();
    gr["shearyz"] = analysis.shearyz.serialize();

  }

  return true;
}








