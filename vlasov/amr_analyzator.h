#pragma once

#include <cmath> 

#include "cell.h"
#include "grid.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"



namespace vlasov {


/// General analyzator that computes moments for the vlasov meshes inside the cells
template<typename T>
class Analyzator {


  /// integrate Adaptive Mesh for number density
  // TODO: think about merging with amr_spatial_solver auxiliary functions
  T integrate(const toolbox::AdaptiveMesh<T,3>& m) const {
    T integ = T(0);

    // pre-create size of the elements
    std::vector<T> du;
    du.resize( m.maximum_refinement_level );

    int rfl; // current refinement level
    for(rfl=0; rfl<=m.maximum_refinement_level; rfl++) {
      auto lens = m.get_length(rfl);
      T len = lens[0]*lens[1]*lens[2];
      du[rfl] = len;
    }

    for(auto cid : m.get_cells(true) ) {
      if( !m.is_leaf(cid) ) continue;

      //auto index = m.get_indices(cid);
      rfl        = m.get_refinement_level(cid);
      //auto uvel  = m.get_center(index, rfl);
      //gam        = gamma<T,3>(uvel);

      integ += m.data.at(cid) * du[rfl];
    }

    return integ;
  }


  public:

  virtual void analyze( vlasov::VlasovCell& cell )
  {

    // Yee lattice reference
    auto& yee = cell.getYee();
    yee.rh.clear();


    // get reference to the Vlasov fluid that we are solving
    auto& step0 = cell.steps.get(0);

    // timestep
    // T dt = cell.dt;
    // T dx = cell.dx;


    // loop over different particle species 
    //int ispc = 0; // ith particle species
    for(auto&& block0 : step0) {
      //T qm = 1.0 / block0.qm;  // charge to mass ratio

      int Nx = int(block0.Nx),
          Ny = int(block0.Ny),
          Nz = int(block0.Nz);

      for (int s=0; s<Nz; s++) {
        for(int r=0; r<Ny; r++) {
          for(int q=0; q<Nx; q++) {
            const auto& M   = block0.block(q,r,s);   // f_i
            T dens = integrate(M);

            yee.rh(q,r,s) += dens; // 
          }
        }
      }

    }// end of loop over species


    return;
  }


};


}

