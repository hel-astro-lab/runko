#include "amr_momentum_solver.h"

#include <cmath> 
//#include <Eigen/Dense>

#include "external/cppitertools/zip.hpp"
#include "definitions.h"

//using namespace Eigen;
using iter::zip;


/// Get snapshot current J_i^n+1 from momentum distribution
template< typename T, int D, int V>
void vlv::MomentumSolver<T,D,V>::update_future_current( vlv::Tile<D>& tile, T cfl)
{
  //auto& gs = tile.get_grids();
  tile.jx1.clear();

  auto& step0 = tile.steps.get(0);
  for(auto&& block0 : step0) {

    auto Nx = int(block0.Nx),
         Ny = int(block0.Ny),
         Nz = int(block0.Nz);

    for (int s=0; s<Nz; s++) {
      for(int r=0; r<Ny; r++) {
        for(int q=0; q<Nx; q++) {
          const auto& M   = block0.block(q,r,s);   // f_i

          T qm = 1.0 / block0.qm;  // charge to mass ratio

          // Jx current; chi(u) = u/gamma = v
          //
          // NOTE: needs to be given in units of grid speed 
          //       so we scale with dt/dx
          //gs.jx1(q,r,s) += qm*
          tile.jx1(q,r,s) += qm*cfl*
            integrate_moment(
                M,
                [](std::array<T,3> uvel) -> T 
                { return uvel[0]/gamma<T,3>(uvel); }
                );
        }
      }
    }

  }// end of loop over species

  }



/*! \brief Solve Vlasov tile contents
 *
 * Exposes the actual momentum mesh from the Vlasov Tile containers
 * and feeds those to the mesh solver.
 */
template< typename T, int D, int V>
void vlv::MomentumSolver<T,D,V>::solve( vlv::Tile<D>& tile, T step_size)
{

  // init lock/mutex for mesh.clear()

  // get reference to the Vlasov fluid that we are solving
  auto& step0 = tile.steps.get(0);
  auto& step1 = tile.steps.get(1);


  // get reference to the Yee grid 
  auto& gs = tile.get_grids();

  // timestep
  //auto dt   = (T) tile.dt;      
  //auto dx   = (T) tile.dx;      
  //T cfl  = step_size*dt/dx;
  auto cfl = step_size*tile.cfl;

  // block limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;


  /// Now get future current
  update_future_current(tile, cfl);

  // param object for solve_mesh
  vlv::tools::Params<T> params = {};
  params.cfl = cfl;


  // loop over different particle species (zips current [0] and new [1] solutions)
  for(auto&& blocks : zip(step0, step1) ) {
      
    // loop over the tile's internal grid
    auto& block0 = std::get<0>(blocks);
    auto& block1 = std::get<1>(blocks);

      

    for(int q=0; q<block0.Nx; q++) {
      for(int r=0; r<block0.Ny; r++) {
        for(int s=0; s<block0.Nz; s++) {
          T qm = 1.0 / block0.qm;  // charge to mass ratio


          // Get local field components
          vec 
            B = 
            {{
               (T) gs.bx(q,r,s),
               (T) gs.by(q,r,s),
               (T) gs.bz(q,r,s)
            }},              

            // E-field interpolated to the middle of the tile
            // E_i = (E_i+1/2 + E_i-1/2)
            // XXX
            E =                
            {{                 
               (T) (0.5*(gs.ex(q,r,s) + gs.ex(q-1,r,   s  ))),
               (T) (0.5*(gs.ey(q,r,s) + gs.ey(q,  r-1, s  ))),
               (T) (0.5*(gs.ez(q,r,s) + gs.ez(q,  r,   s-1)))
            }};

          // Now push E field to future temporarily
          //E[0] -= gs.jx1(q,r,s) * 0.5;


          // dig out velomeshes from blocks
          auto& mesh0 = block0.block(q,r,s);
          auto& mesh1 = block1.block(q,r,s);

          // fmt::print("solving for srq ({},{},{})\n",s,r,q);
          params.qm = qm;
          params.xloc = mins[0] + static_cast<T>(q);

          // then the final call to the actual mesh solver
          solve_mesh( mesh0, mesh1, E, B, params);
        }
      }
    }
  }

  // XXX update jx1 for debug
  //update_future_current(tile, cfl);
  //
  

  }

//--------------------------------------------------
// explicit template instantiation
template class vlv::MomentumSolver<float_m, 1, 1>;


