#pragma once

#include <cmath> 

#include "tile.h"
#include "grid.h"
#include "../em-fields/tile.h"
#include "amr/mesh.h"
#include "amr/integrate.h"

#include "../tools/signum.h"

#include <iostream>

using toolbox::sign;


namespace vlv {

/// Relativistic gamma from velocity
template<typename T, int V>
inline T gamma(std::array<T,V>& uvel) 
{
  T gammasq = 1.0;
  for(size_t i=0; i<V; i++) gammasq += uvel[i]*uvel[i];
  return std::sqrt(gammasq);
}

//{
//  return 1.0;
//}



/// General analyzator that computes moments for the vlasov meshes inside the tiles
template<typename T>
class Analyzator {

  public:

  Analyzator() {};

  virtual ~Analyzator() = default;


  virtual void analyze( vlv::Tile<1>& tile )
  {

    // Yee lattice reference
    auto& yee = tile.getYee();
    yee.rho.clear();
    //yee.ekin.clear();
    //yee.jx1.clear();


    // get reference to the Vlasov fluid that we are solving
    auto& step0 = tile.steps.get(0);

    // timestep
    // T dt = tile.dt;
    // T dx = tile.dx;


    // loop over different particle species 
    int ispc = 0; // ith particle species
    for(auto&& block0 : step0) {

      // get reference to species-specific analysis mesh
      auto& analysis = tile.analysis[ispc];

      auto Nx = static_cast<int>(block0.Nx),
           Ny = static_cast<int>(block0.Ny),
           Nz = static_cast<int>(block0.Nz);

        for(int s=0; s<Nz; s++) {
        for(int r=0; r<Ny; r++) {
        for(int q=0; q<Nx; q++) {
            const auto& M   = block0.block(q,r,s);   // f_i

            //T qm = 1.0 / block0.qm;  // charge to mass ratio


            // TODO there is a possibility to optimize this by putting everything 
            // inside one loop at the expense of code readability... 
            // ...So it is not done here atm.

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
            // NOTE: different insert location (i.e., analysis object 
            // instead of yee object)

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
            analysis.Vx(q,r,s) = Vx / rho;

            // TODO add Vy and Vz components

              
            /// (non-relativistic) temperature
            // chi(u) = v - V
            T Tx = 
              integrate_moment(
                  M,
                [Vx](std::array<T,3> uvel) -> T 
                { return std::pow(uvel[0]/gamma<T,3>(uvel) - Vx, 2); }
                );
            analysis.Tx(q,r,s) = std::sqrt(Tx);

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

