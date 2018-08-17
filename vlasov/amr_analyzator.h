#pragma once

#include <cmath> 

#include "cell.h"
#include "grid.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"

#include "../tools/signum.h"

#include <iostream>

using toolbox::sign;


namespace vlasov {

/// Relativistic gamma from velocity
template<typename T, int D>
inline T gamma(std::array<T,D>& uvel) 
{
  T gammasq = 1.0;
  for(size_t i=0; i<D; i++) gammasq += uvel[i]*uvel[i];
  return std::sqrt(gammasq);
}
//{
//  return 1.0;
//}



/// integrate Adaptive Mesh phase space with function Chi(u)
// TODO: think about merging with amr_spatial_solver auxiliary functions
template<typename T, typename Lambda>
T integrate_moment(
    const toolbox::AdaptiveMesh<T,3>& m,
    Lambda&& chi
    ) {
  auto integ = T(0);

  // pre-create size of the elements
  // TODO optimize this depending on top_refinement_level
  //std::vector<T> du;
  //du.resize( m.top_refinement_level+1 );
  int rfl; // current refinement level
  //for(rfl=0; rfl<=m.top_refinement_level; rfl++) {
  //  auto lens = m.get_length(rfl);
  //  T len = lens[0]*lens[1]*lens[2];
  //  du[rfl] = len;
  //}

  // simplified approach assuming no refinement.
  // TODO: rewrite this to accommade adaptivity
  auto lens = m.get_length(0);
  T du0 = lens[0]*lens[1]*lens[2];


  // convolve with chi(u) function
  for(auto cid : m.get_cells(false) ) {
    if( !m.is_leaf(cid) ) continue; //TODO: fixme

    auto index = m.get_indices(cid);
    rfl        = m.get_refinement_level(cid);
    auto uvel  = m.get_center(index, rfl);

    //std::cout << "cid:" << cid << " data:" << m.data.at(cid) << " du:" << du[rfl] << " chi:" << chi(uvel) << '\n';
    //std::cout << "cid:" << cid << " data:" << m.get(cid) << " du:" << du[rfl] << " chi:" << chi(uvel) << '\n';

    //integ += m.data.at(cid) * du[rfl] * chi(uvel);
    //integ += m.data.at(cid) * du[rfl] * chi(uvel);
    integ += m.data.at(cid) * du0 * chi(uvel);
  }

  return integ;
}




/// General analyzator that computes moments for the vlasov meshes inside the cells
template<typename T>
class Analyzator {


  public:

  virtual void analyze( vlasov::VlasovCell& cell )
  {

    // Yee lattice reference
    auto& yee = cell.getYee();
    yee.rho.clear();
    //yee.ekin.clear();
    //yee.jx1.clear();


    // get reference to the Vlasov fluid that we are solving
    auto& step0 = cell.steps.get(0);

    // timestep
    // T dt = cell.dt;
    // T dx = cell.dx;


    // loop over different particle species 
    int ispc = 0; // ith particle species
    for(auto&& block0 : step0) {

      // get reference to species-specific analysis mesh
      auto& analysis = cell.analysis[ispc];

      auto Nx = int(block0.Nx),
          Ny = int(block0.Ny),
          Nz = int(block0.Nz);

      for(int s=0; s<Nz; s++) {
        for(int r=0; r<Ny; r++) {
          for(int q=0; q<Nx; q++) {
            const auto& M   = block0.block(q,r,s);   // f_i

            //T qm = 1.0 / block0.qm;  // charge to mass ratio


            // TODO there is a possibility to optimize this by putting everything 
            // inside one loop at the expense of code readability... So it is not done.

            //-------------------------------------------------- 
            // Global quantities (sum over species)

            // number density; chi(u) = 1
            yee.rho(q,r,s) += 
              integrate_moment(
                  M,
                [](std::array<T,3>& ) -> T { return T(1);}
                );

            // Jx current; chi(u) = u/gamma = v
            // NOTE: Current is already computed in the time loop so omit it here
            //yee.jx1(q,r,s) += sign(qm)*
            //  integrate_moment(
            //      M,
            //    [](std::array<T,3> uvel) -> T 
            //    { return uvel[0]/gamma<T,3>(uvel); }
            //    );
              

            //-------------------------------------------------- 
            // Species-specific quantities
            // NOTE: different insert location

            // number density; chi(u) = 1
            T rho = 
              integrate_moment(
                  M,
                [](std::array<T,3>& ) -> T { return T(1); }
                );
            analysis.rho(q,r,s) = rho;


            // mean gamma (fluid energy density)
            // chi(u) = \gamma / rho
            analysis.mgamma(q,r,s) = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T { return gamma<T,3>(uvel); }
                );
            analysis.mgamma(q,r,s) /= rho; // normalize


            /// (mean) fluid/bulk velocity
            // chi(u) = v = u/gamma
            T Vx = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T 
                { return uvel[0]/gamma<T,3>(uvel); }
                );
            analysis.Vx(q,r,s) = Vx;

            // TODO add Vy and Vz components

              
            /// (non-relativistic) temperature
            // chi(u) = v - V
            analysis.Tx(q,r,s) = 
              integrate_moment(
                  M,
                [Vx](std::array<T,3> uvel) -> T 
                { return uvel[0]/gamma<T,3>(uvel) - Vx; }
                );

            // TODO add Ty and Tz components


            // Kinetic energy
            // chi(u) = 1/2 u.u
            analysis.ekin(q,r,s) = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T { 
                return 0.5*uvel[0]*uvel[0]/(gamma<T,3>(uvel)*gamma<T,3>(uvel)); }
                );


          }
        }
      }

      ispc++;
    }// end of loop over species


    return;
  }


};


}

