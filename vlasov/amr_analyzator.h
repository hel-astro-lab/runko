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

            /// (mean) non-relativistic fluid/bulk velocity
            // chi(v) = V = U/gamma 
            T Vx = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T 
                { return uvel[0]/gamma<T,3>(uvel); }
                );
            analysis.Vx(q,r,s) = Vx / rho;

            //--------------------------------------------------
            // stress energy tensor components
              
            // T^00: fluid energy density
            // chi(u) = gam*mass*c^2
            analysis.edens(q,r,s) = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T { return gamma<T,3>(uvel); }
                );
              
            /// T^01 = T^10: momentum density, chi(u) = u_x*mass*c = m u_x
            analysis.momx(q,r,s) = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T 
                { return uvel[0]; }
                );

            /// T^11: pressure, chi(u) = m u_x v_x *c^2
            analysis.pressx(q,r,s) = 
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T 
                { return uvel[0]*uvel[0]/gamma<T,3>(uvel); }
                );

            /// T^12 = T^21: shear, chi(u) = m u_x v_y c^2
            //analysis.shearxy(q,r,s) = 
            //  integrate_moment(
            //      M,
            //    [](std::array<T,3> uvel) -> T 
            //    { return uvel[0]*uvel[1]/gamma<T,3>(uvel); }
            //    );

            // + other components for u_y and u_z

            //--------------------------------------------------
              
            // mean gamma
            //analysis.mgamma(q,r,s) = analysis.edens(q,r,s)/rho;

            // TODO: correct?
            /// relativistic temperature
            // chi(u) = gam - 1
            analysis.temp(q,r,s) =  
              integrate_moment(
                  M,
                [](std::array<T,3> uvel) -> T 
                { return gamma<T,3>(uvel) - T(1); }
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

