#pragma once

#include <cmath> 

#include "cell.h"
#include "grid.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"
#include "amr/numerics.h"
#include "amr/refiner.h"
#include "amr/operators.h"
#include "../tools/signum.h"
#include "amr_analyzator.h"

#include "../tools/cppitertools/zip.hpp"
using iter::zip;
using toolbox::sign;


#include <Eigen/Dense>
using namespace Eigen;


namespace vlasov {




/*! General interface for AMR spatial solvers
 *
 */
template<typename T>
class SpatialSolver {
    
  /// get neighboring cell from grid
  // TODO: think about dynamic casting; something might be wrong in the design
  // TODO: separate into own communication module/header
  //vlasov::PlasmaBlock& get_external_data(
  //    int i, int j, int ispc,
  //    vlasov::VlasovCell& cell, 
  //    vlasov::Grid& grid)
  //{
  //  auto neigh_index   = cell.neighs(i, j); 
  //  uint64_t neigh_cid = grid.cellId( std::get<0>(neigh_index), std::get<1>(neigh_index) );
  //  vlasov::VlasovCell& cell_neigh = dynamic_cast<vlasov::VlasovCell&>( grid.getCell(neigh_cid) );

  //  auto& species = cell_neigh.steps.get();

  //  return species[ispc];
  //}


  public:
    typedef std::array<T, 3> vec;

    /// Actual solver implementation
    virtual void solve( vlasov::VlasovCell& cell, vlasov::Grid& grid) = 0;
};





/// Relativistic gamma from velocity
//template<typename T, int D>
//inline T gamma(std::array<T,D>& uvel) 
//{
//  T gammasq = 1.0;
//  for(size_t i=0; i<D; i++) gammasq += uvel[i]*uvel[i];
//  return std::sqrt(gammasq);
//}


/// Simple (relativistic) box integration of a mesh flux
/*
template<typename T>
T integrate_current(
    toolbox::AdaptiveMesh<T,3>& m)
{

  int rfl; // current refinement level
  T gam; // differential element (box size)

  // pre-create size of the elements
  std::vector<T> du;
  du.resize( m.maximum_refinement_level );

  for(rfl=0; rfl<=m.maximum_refinement_level; rfl++) {
    auto lens = m.get_length(rfl);
    T len = lens[0]*lens[1]*lens[2];
    du[rfl] = len;
  }

  // integrate leafs; j = int{ U du/gamma}
  T integ = T(0);
  for(auto cid : m.get_cells(false) ) {
    if( !m.is_leaf(cid) ) continue;

    auto index = m.get_indices(cid);
    rfl        = m.get_refinement_level(cid);
    auto uvel  = m.get_center(index, rfl);
    gam        = gamma<T,3>(uvel);

    integ += m.data.at(cid)*du[rfl]/gam;
  }

  return integ;
}
*/



/// Lagrangian conservative solver
template<typename T>
class AmrSpatialLagrangianSolver : public SpatialSolver<T> {

  public:
    //using SpatialSolver<T>::get_external_data;

    /// velocity tensor product
    toolbox::AdaptiveMesh<T,3>& velocityTensorProduct(
        toolbox::AdaptiveMesh<T,3>&  lhs,
        size_t rank,
        T Const)
    {

      for(auto&& cid : lhs.get_cells(false)) {
        auto index = lhs.get_indices(cid);
        int rfl    = lhs.get_refinement_level(cid);
        auto uvel  = lhs.get_center(index, rfl);

        // compute (u_x/gamma)^R factor
        T vvel = T(1);
        T gam  = gamma<T,3>(uvel);
        for(size_t i=0; i<rank; i++) {
          vvel *= (uvel[0]/gam);
        }

        // finally compute flux: C*(u_x dt/dx)^R f(x,v,t)
        lhs.data.at(cid) *= vvel*Const;
      }

      return lhs;
    }

    

    /// Linear 1st order upwind-biased flux
    // inline toolbox::AdaptiveMesh<T,3> flux1st(
    //     const toolbox::AdaptiveMesh<T,3>& M0,
    //     const toolbox::AdaptiveMesh<T,3>& Mm1) {

    //   return M0.xx + Mm1.xx;
    // }
 

    /// Second order flux; U_+
    inline toolbox::AdaptiveMesh<T,3> flux2nd(
        const toolbox::AdaptiveMesh<T,3>& M ,
        const toolbox::AdaptiveMesh<T,3>& Mp1,
        T cfl) {

      auto Mp = Mp1 + M; // explicitly allocate tmp variables
      auto Mn = Mp1 - M;    

      return   velocityTensorProduct(Mp, 1, 0.5*cfl ) 
             - velocityTensorProduct(Mn, 2, 0.5*cfl*cfl ); 
    }




    /// sweep and solve internal blocks in X direction
    void xsweep(
        vlasov::PlasmaBlock& block0,
        vlasov::PlasmaBlock& block1,
        vlasov::PlasmaBlock& block0_left,
        vlasov::PlasmaBlock& block0_right,
        T qm,
        T cfl,
        fields::YeeLattice& yee)
    {


      // initialize new step
      for (size_t s=0; s<block0.Nz; s++) {
        for(size_t r=0; r<block0.Ny; r++) {
          for(size_t q=0; q<block0.Nx; q++) {
            const auto& M = block0.block(q,r,s); // f_i
            auto& N       = block1.block(q,r,s); // f_i^t+dt
            N = M;
          }
        }
      }


      // local flows
      toolbox::AdaptiveMesh<T,3> flux;
      int Nx = int(block0.Nx),
          Ny = int(block0.Ny),
          Nz = int(block0.Nz);

      for (int s=0; s<Nz; s++) {
        for(int r=0; r<Ny; r++) {
          for(int q=-1; q<Nx; q++) {

            // dig out velomeshes from blocks (M's; constants)
            // and calculate flux, i.e., U_i+1/2
            //
            // NOTE: flux has to be calculated inside every if-clause
            //       because we can not pre-initialize references 
            //       (these are there to save memory by avoiding 
            //       copy of consts)
              
            if (q < 0) { //left halo
              const auto& M   = block0_left.block(block0_left.Nx-1,r,s); // f_i
              const auto& Mp1 = block0.block(0,r,s);                     // f_i+1

              flux = flux2nd(M, Mp1, cfl);
            } else if ( (q >= 0) && (q <= Nx-2) ) { // inside
              const auto& M   = block0.block(q,r,s);   // f_i
              const auto& Mp1 = block0.block(q+1,r,s); // f_i+1

              flux = flux2nd(M, Mp1, cfl);
            } else if (q >= Nx-1) { //right halo
              const auto& M   = block0.block(q,r,s);       // f_i
              const auto& Mp1 = block0_right.block(0,r,s); // f_i+1

              flux = flux2nd(M, Mp1, cfl);
            }

            // new local time step targets to update into (N's)
            auto& N   = block1.block(q,  r,s); // f_i  ^t+dt
            auto& Np1 = block1.block(q+1,r,s); // f_i+1^t+dt

            // now flow to neighbors; only local flows are allowed
            if(q >= 0)    N   -= flux; // - U_i+1/2 (outflowing from cell)
            if(q <= Nx-2) Np1 += flux; // + U_i-1/2 (inflowing to neighbor)

            // calculate current
            T jx = qm*integrate_moment( 
                flux,
                [](std::array<T,3>& uvel) -> T { return 1.0;}
                );
              
            if(q >= 0)    yee.jx(q,r,s)   += jx; //U_i+1/2
            //if(q <= Nx-2) yee.jx(q+1,r,s) += jx; //U_i-1/2
          }
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
      T cfl = cell.dt/cell.dx;


      // loop over different particle species (zips current [0] and new [1] solutions)
      int ispc = 0; // ith particle species
      for(auto&& blocks : zip(step0, step1) ) {

        // internal blocks of vlasov meshes for current and future timestep
        auto& block0 = std::get<0>(blocks); // Solution^t
        auto& block1 = std::get<1>(blocks); // Solution^t+dt
        T qm = 1.0 / block0.qm;  // charge to mass ratio


        // external neighbors
        auto& block0_left   = cell.get_external_data(-1, 0, ispc, grid);
        auto& block0_right  = cell.get_external_data(+1, 0, ispc, grid);
        //auto& block0_bottom = get_external_data( 0,-1, ispc, cell, grid);
        //auto& block0_top    = get_external_data( 0,+1, ispc, cell, grid);

        // sweep in X
        xsweep(block0, block1, block0_left, block0_right, qm, cfl, yee);
        // ysweep(block0, block1, block0_bottom, block0_top,   qm, dt, dx);
        // xsweep(block0, block1, block0_left,   block0_right, qm, dt/2, dx);


        ispc++;
      }

      // done
    }


};





} // end of vlasov namespace
