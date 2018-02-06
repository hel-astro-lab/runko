#pragma once

#include <cmath> 

#include "cell.h"
#include "grid.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"
#include "amr/numerics.h"
#include "amr/refiner.h"
#include "amr/operators.h"

#include "../tools/cppitertools/zip.hpp"
using iter::zip;


#include <Eigen/Dense>
using namespace Eigen;


namespace vlasov {


/*! General interface for AMR spatial solvers
 *
 */
template<typename T>
class SpatialSolver {

  public:
    typedef std::array<T, 3> vec;


    /// get neighboring cell from grid
    // TODO: think about dynamic casting; something might be wrong in the design
    vlasov::PlasmaBlock& get_external_data(
        int i, int j, int ispc,
        vlasov::VlasovCell& cell, 
        vlasov::Grid& grid)
    {
      auto neigh_index               = cell.neighs(i, j); 
      uint64_t neigh_cid             = grid.cellId( std::get<0>(neigh_index), std::get<1>(neigh_index) );
      vlasov::VlasovCell& cell_neigh = dynamic_cast<vlasov::VlasovCell&>( grid.getCell(neigh_cid) );

      auto& species = cell_neigh.steps.get();

      return species[ispc];
    }


    /// Actual solver implementation
    virtual void solve( vlasov::VlasovCell& cell, vlasov::Grid& grid) = 0;

};





/// Relativistic gamma from velocity
template<typename T, int D>
inline T gammafac(std::array<T,D>& uvel) 
{
  T gammasq = 1.0;
  for(size_t i=0; i<3; i++) gammasq += uvel[i]*uvel[i];
  return std::sqrt(gammasq);
}







/// Lagrangian conservative solver
template<typename T>
class AmrSpatialLagrangianSolver : public SpatialSolver<T> {

  public:

    using SpatialSolver<T>::get_external_data;



    /// velocity tensor product
    toolbox::AdaptiveMesh<T,3>& velocityTensorProduct(
        toolbox::AdaptiveMesh<T,3>&  lhs,
        size_t rank)
    {

      for(auto&& cid : lhs.get_cells(false)) {
        auto index = lhs.get_indices(cid);
        int rfl    = lhs.get_refinement_level(cid);
        auto uvel  = lhs.get_center(index, rfl);

        // compute (u_x/gamma)^R factor
        T vvel = T(1);
        T gam  = gamma(uvel);
        for(size_t i=0; i<rank; i++) {
          vvel *= (uvel[0]/gam);
        }

        // finally compute flux: (u_x dt/dx)^R f(x,v,t)
        lhs.data.at(cid) *= vvel;
      }


      return lhs;
    }

    

    /// Linear 1st order upwind-biased flux
    // inline toolbox::AdaptiveMesh<T,3> flux1st(
    //     const toolbox::AdaptiveMesh<T,3>& M0,
    //     const toolbox::AdaptiveMesh<T,3>& Mm1) {

    //   return M0.xx + Mm1.xx;
    // }
 

    /// Second order centralized flux; U_+
    inline toolbox::AdaptiveMesh<T,3> flux2nd(
        const toolbox::AdaptiveMesh<T,3>& M ,
        const toolbox::AdaptiveMesh<T,3>& Mp1,
        T dt,
        T dx) {

      return   0.5*(dt/dx)*        velocityTensorProduct(Mp1 + M, 1) 
             - 0.5*(dt/dx)*(dt/dx)*velocityTensorProduct(Mp1 - M, 2); 
    }



    /// sweep and solve internal blocks in X direction
    void xsweep(
        vlasov::PlasmaBlock& block0,
        vlasov::PlasmaBlock& block1,
        vlasov::PlasmaBlock& block0_left,
        vlasov::PlasmaBlock& block0_right,
        T qm,
        T dt,
        T dx)
    {

      size_t H=1; // width of halo region
      for (size_t s=0; s<block0.Nz; s++) {
        for(size_t r=0; r<block0.Ny; r++) {


          for(size_t q=0+H; q<block0.Nx-H; q++) {

            // dig out velomeshes from blocks (M's; constants)
            // const auto& Mm1 = block0.block(q-1,r,s); // f_i-1
            const auto& M   = block0.block(q,  r,s); // f_i
            const auto& Mp1 = block0.block(q+1,r,s); // f_i+1


            // new time step targets to update into (N's)
            auto& N   = block1.block(q,  r,s); // f_i  ^t+dt
            auto& Nm1 = block1.block(q-1,r,s); // f_i-1^t+dt


            // debugging test
            N = M + Mp1;


            // flux calculation, i.e., U_i+1/2
            // auto flux = flux2nd(M, Mp1, dt, dx);

            // N   -= flux; // - U_i+1/2
            // Nm1 += flux; // + U_i-1/2

            // calculate current



          }


          /*
          const auto& M   = block0.block(block0.Nx-1,  r,s); // f_i
          const auto& Mp1 = block0_right.block(0,r,s); // f_i-1
          auto& N   = block1.block(,  r,s); // f_i  ^t+dt

          auto flux = flux2nd(Mm1, M, dt, dx);
          N += flux; // + U_i-1/2
          */



        }
      }




    }


    // Strang splitted rotating (X/2 Y X/2) solver 
    virtual void solve( vlasov::VlasovCell& cell, vlasov::Grid& grid )
    {

      // Yee lattice reference
      auto& yee = cell.getYee();
      yee.jx.clear();
      yee.jy.clear();
      yee.jz.clear();


      // get reference to the Vlasov fluid that we are solving
      auto& step0 = cell.steps.get(0);
      auto& step1 = cell.steps.get(1);

      // timestep
      T dt = cell.dt;
      T dx = cell.dx;


      // loop over different particle species (zips current [0] and new [1] solutions)
      int ispc = 0; // ith particle species
      for(auto&& blocks : zip(step0, step1) ) {

        // internal blocks of vlasov meshes for current and future timestep
        auto& block0 = std::get<0>(blocks); // Solution^t
        auto& block1 = std::get<1>(blocks); // Solution^t+dt
        T qm = 1.0 / block0.qm;  // charge to mass ratio


        // external neighbors
        auto& block0_left   = get_external_data(-1, 0, ispc, cell, grid);
        auto& block0_right  = get_external_data(+1, 0, ispc, cell, grid);
        //auto& block0_bottom = get_external_data( 0,-1, ispc, cell, grid);
        //auto& block0_top    = get_external_data( 0,+1, ispc, cell, grid);

        // sweep in X
        xsweep(block0, block1, block0_left, block0_right, qm, dt, dx);


        // ysweep(block0, block1, block0_bottom, block0_top,   qm, dt, dx);
        // xsweep(block0, block1, block0_left,   block0_right, qm, dt/2, dx);

        ispc++;
      }

      // done

    }


};







} // end of vlasov namespace
