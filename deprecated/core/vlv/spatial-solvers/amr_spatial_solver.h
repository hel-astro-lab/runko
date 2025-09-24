#pragma once

#include <cmath> 

#include "core/emf/tile.h"

#include "core/vlv/tile.h"
#include "core/vlv/grid.h"
#include "core/vlv/amr/mesh.h"
#include "core/vlv/amr/numerics.h"
#include "core/vlv/amr/refiner.h"
#include "core/vlv/amr/operators.h"

#include "tools/signum.h"
#include "external/cppitertools/zip.hpp"

using iter::zip;
using toolbox::sign;


//#include <Eigen/Dense>
//using namespace Eigen;
using std::min;
using std::max;


namespace vlv {



/*! General interface for AMR spatial solvers
 *
 */
template<typename T>
class SpatialSolver {
    
  /// get neighboring tile from grid
  // TODO: think about dynamic casting; something might be wrong in the design
  // TODO: separate into own communication module/header
  //vlv::PlasmaBlock& get_external_data(
  //    int i, int j, int ispc,
  //    vlv::VlasovTile& tile, 
  //    vlv::Grid& grid)
  //{
  //  auto neigh_index   = tile.neighs(i, j); 
  //  uint64_t neigh_cid = grid.tile_id( std::get<0>(neigh_index), std::get<1>(neigh_index) );
  //  vlv::VlasovTile& tile_neigh = dynamic_cast<vlv::VlasovTile&>( grid.get_tile(neigh_cid) );

  //  auto& species = tile_neigh.steps.get();

  //  return species[ispc];
  //}


  public:

    SpatialSolver() = default;

    virtual ~SpatialSolver() = default;

    using vec = std::array<T, 3>;

    /// Actual solver implementation
    virtual void solve( vlv::Tile<1>& tile, corgi::Grid<1>& grid) = 0;
};




/// Lagrangian conservative solver
template<typename T>
class AmrSpatialLagrangianSolver : 
  virtual public SpatialSolver<T> 
{

  public:

    AmrSpatialLagrangianSolver() = default;

    ~AmrSpatialLagrangianSolver() override = default;
      

    //using SpatialSolver<T>::get_external_data;

    /// velocity tensor product
    toolbox::AdaptiveMesh<T,3>& velocity_tensor_product(
        toolbox::AdaptiveMesh<T,3>&  lhs,
        size_t rank,
        T Const)
    {

      for(auto& it : lhs.data) {
        auto index = lhs.get_indices(it.first);
        int rfl    = lhs.get_refinement_level(it.first);
        auto uvel  = lhs.get_center(index, rfl);

        // compute (u_x/gamma)^R factor
        T vvel = static_cast<T>(1);
        T gam  = gamma<T,3>(uvel);
        for(size_t i=0; i<rank; i++) {
          vvel *= (uvel[0]/gam);
        }

        // finally compute flux: C*(u_x dt/dx)^R f(x,v,t)
        it.second *= vvel*Const;
      }

      return lhs;
    }

    

    /// Linear 1st order upwind-biased flux
    inline toolbox::AdaptiveMesh<T,3> flux1st(
        const toolbox::AdaptiveMesh<T,3>& M0,
        const toolbox::AdaptiveMesh<T,3>& Mp1,
        T cfl) {

      // make a new fresh mesh for updating (based on M0)
      // XXX shallow or deep copy?
      toolbox::AdaptiveMesh<T,3> flux(M0);
      flux.data.clear();

      // Here we compute (u_x/gamma) and then
      // switch to units of grid speed by multiplying with Cfl

      //left side
      for(const auto& it : M0.data) {
        auto& cid  = it.first;
        auto index = M0.get_indices(cid);
        int rfl    = M0.get_refinement_level(cid);
        auto uvel  = M0.get_center(index, rfl);

        if (-uvel[0] < 0.0) {
          T gam  = gamma<T,3>(uvel);
          //T gam = 1.0;
          flux.data[cid] = cfl*(uvel[0]/gam)*it.second;

          //if (M0.data.at(cid) > 0.01) std::cout << "LL:" << flux.data[cid] << " / " << M0.data.at(cid) << std::endl;
        } 

      }

      //right side
      for(const auto& it : Mp1.data) {
        auto& cid  = it.first;
        auto index = Mp1.get_indices(cid);
        int rfl    = Mp1.get_refinement_level(cid);
        auto uvel  = Mp1.get_center(index, rfl);

        if (-uvel[0] > 0.0) {
          T gam  = gamma<T,3>(uvel);
          //T gam = 1.0;
          flux.data[cid] = cfl*(uvel[0]/gam)*it.second;

          //if (Mp1.data.at(cid) > 0.01) std::cout << "RR:" << uvel[0] << std::endl;
        }
      }

      return flux; 
    }
 


    //--------------------------------------------------
    //inline T fmaxf(T fm2, T fm1, T f0, T fp1, T fp2) { 
    //  return T(100);
    //}
    //inline T fminf(T fm2, T fm1, T f0, T fp1, T fp2) { 
    //  return T(0);
    //}
      
    inline T fmaxf(T fm2, T fm1, T f0, T fp1, T fp2) { 
      T fmax1 = max<T>( max<T>(fm1, f0), min<T>( T(2)*fm1 - fm2, T(2)*f0 - fp1) );
      T fmax2 = max<T>( max<T>(fp1, f0), min<T>( T(2)*fp1 - fp2, T(2)*f0 - fm1) );
      return max<T>(fmax1, fmax2);
    }

    inline T fminf(T fm2, T fm1, T f0, T fp1, T fp2) { 
      T fmin1 = min<T>( min<T>(fm1, f0), max<T>(T(2)*fm1 - fm2, T(2)*f0 - fp1) );
      T fmin2 = min<T>( min<T>(fp1, f0), max<T>(T(2)*fp1 - fp2, T(2)*f0 - fm1) );
      return max<T>( T(0), min<T>(fmin1, fmin2) );
    }

    inline T Lpf(T fm2, T fm1, T f0, T fp1, T fp2) {
      if (fp1 >= f0){ 
        T fmin = fminf(fm2, fm1, f0, fp1, fp2);
        return min<T>( T(2)*(f0 - fmin), fp1 - f0);
      } 
        T fmax = fmaxf(fm2, fm1, f0, fp1, fp2);
        return max<T>( T(2)*(f0 - fmax), fp1 - f0);
      
    }

    inline T Lmf(T fm2, T fm1, T f0, T fp1, T fp2) {
      if (f0 >= fm1){ 
        T fmax = fmaxf(fm2, fm1, f0, fp1, fp2);
        return min<T>( T(2)*(fmax - f0), f0 - fm1);
      } 
        T fmin = fminf(fm2, fm1, f0, fp1, fp2);
        return max<T>( T(2)*(fmin - f0), f0 - fm1);
      
    }



    /// 3rd order upwind-biased flux with flux limiters
    inline toolbox::AdaptiveMesh<T,3> flux3rdU(
        const toolbox::AdaptiveMesh<T,3>& Mm2,
        const toolbox::AdaptiveMesh<T,3>& Mm1,
        const toolbox::AdaptiveMesh<T,3>& M0,
        const toolbox::AdaptiveMesh<T,3>& Mp1,
        const toolbox::AdaptiveMesh<T,3>& Mp2,
        const toolbox::AdaptiveMesh<T,3>& Mp3,
        T cfl) {

      // make a new fresh mesh for updating (based on M0)
      // XXX shallow or deep copy?
      toolbox::AdaptiveMesh<T,3> flux(M0);
      flux.data.clear();

      // Here we compute (u_x/gamma) and then
      // switch to units of grid speed by multiplying with Cfl

      //left side
      for(const auto& it : Mm1.data) {
        auto& cid  = it.first;
        auto index = Mm1.get_indices(cid);
        int rfl    = Mm1.get_refinement_level(cid);
        auto uvel  = Mm1.get_center(index, rfl);

        if (uvel[0] >= 0.0) {
          auto gam  = gamma<T,3>(uvel);
          //T gam = 1.0;
          T v   = cfl*uvel[0]/gam; // CFL

          T fm2 = Mm2.get(cid);
          T fm1 = it.second;
          T f0  =  M0.get(cid);
          T fp1 = Mp1.get(cid);
          T fp2 = Mp2.get(cid);

          //T Lp = fp1 - f0;
          //T Lm = f0  - fm1;
            
          //minmax limited values
          T Lp = Lpf(fm2, fm1, f0, fp1, fp2);
          T Lm = Lmf(fm2, fm1, f0, fp1, fp2);

          flux.data[cid] = 
            v*f0 + 
            v*(1.0-v)*(2.0-v)*Lp/6.0 +
            v*(1.0-v)*(1.0+v)*Lm/6.0;
        }
      }

      //right side
      for(const auto& it : Mp1.data) {
        auto& cid  = it.first;
        auto index = M0.get_indices(cid);
        int rfl    = M0.get_refinement_level(cid);
        auto uvel  = M0.get_center(index, rfl);

        if (uvel[0] < 0.0) {
          auto gam  = gamma<T,3>(uvel);
          //T gam = 1.0;
          // v = i - x
          T v = cfl*uvel[0]/gam; // CFL

          //T fm2 = Mm2.get(cid);
          T fm1 = Mm1.get(cid);
          T f0  = M0.get(cid);
          T fp1 = it.second;
          T fp2 = Mp2.get(cid);
          T fp3 = Mp3.get(cid);
            
          //T Lp = f0  - fp1;
          //T Lm = fp1 - fp2;
            
          //minmax limited values
          T Lp = Lpf(fp3, fp2, fp1, f0, fm1);
          T Lm = Lmf(fp3, fp2, fp1, f0, fm1);

          flux.data[cid] = 
            v*fp1 + 
            v*(1.0+v)*(2.0+v)*Lp/6.0 +
            v*(1.0-v)*(1.0+v)*Lm/6.0;

        }
      }

      return flux; 
    }
 


    /// Second order flux; U_+
    inline toolbox::AdaptiveMesh<T,3> flux2nd(
        const toolbox::AdaptiveMesh<T,3>& M ,
        const toolbox::AdaptiveMesh<T,3>& Mp1,
        T cfl) {

      // explicitly allocate tmp variables
      auto Mp = Mp1 + M; 
      auto Mn = Mp1 - M;    

      return   velocity_tensor_product(Mp, 1, 0.5*cfl ) 
             - velocity_tensor_product(Mn, 2, 0.5*cfl*cfl ); 
    }




    /// sweep and solve internal blocks in X direction
    void xsweep(
        vlv::PlasmaBlock& block0,
        vlv::PlasmaBlock& block1,
        vlv::PlasmaBlock& block0_left,
        vlv::PlasmaBlock& block0_right,
        T qm,
        T cfl,
        emf::Grids& gs)
    {


      // initialize new step
      for (int s=0; s<block0.Nz; s++) {
        for(int r=0; r<block0.Ny; r++) {
          for(int q=0; q<block0.Nx; q++) {
            const auto& M = block0.block(q,r,s); // f_i
            auto& N       = block1.block(q,r,s); // f_i^t+dt
            N = M;
            N.data = M.data;
          }
        }
      }


      // local flows
      toolbox::AdaptiveMesh<T,3> flux;
      auto Nx = int(block0.Nx),
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
              
            // 2nd order central scheme
            //if (q < 0) { //left halo
            //  const auto& M   = block0_left.block(block0_left.Nx-1,r,s); // f_i
            //  const auto& Mp1 = block0.block(0,r,s);                     // f_i+1

            //  flux = flux2nd(M, Mp1, cfl);
            //} else if ( (q >= 0) && (q <= Nx-2) ) { // inside
            //  const auto& M   = block0.block(q,r,s);   // f_i
            //  const auto& Mp1 = block0.block(q+1,r,s); // f_i+1

            //  flux = flux2nd(M, Mp1, cfl);
            //} else if (q >= Nx-1) { //right halo
            //  const auto& M   = block0.block(q,r,s);       // f_i
            //  const auto& Mp1 = block0_right.block(0,r,s); // f_i+1

            //  flux = flux2nd(M, Mp1, cfl);
            //}

            // 1st order upwind scheme
            /*
            if (q < 0) { //left halo
              const auto& M   = block0_left.block(block0_left.Nx-1,r,s); // f_i
              const auto& Mp1 = block0.block(0,r,s);                     // f_i+1

              flux = flux1st(M, Mp1, cfl);
            } else if ( (q >= 0) && (q <= Nx-2) ) { // inside
              const auto& M   = block0.block(q,r,s);   // f_i
              const auto& Mp1 = block0.block(q+1,r,s); // f_i+1

              flux = flux1st(M, Mp1, cfl);
            } else if (q >= Nx-1) { //right halo
              const auto& M   = block0.block(q,r,s);       // f_i
              const auto& Mp1 = block0_right.block(0,r,s); // f_i+1

              flux = flux1st(M, Mp1, cfl);
            }
            */

            // 3rd order upwind scheme
            if (q == -1) { //left halo
              const auto& Mm2 = block0_left.block(block0_left.Nx-3,r,s);   // f_i-2
              const auto& Mm1 = block0_left.block(block0_left.Nx-2,r,s);   // f_i-1
              const auto& M   = block0_left.block(block0_left.Nx-1,r,s);   // f_i
              const auto& Mp1 = block0.block(q+1,r,s);                     // f_i+1
              const auto& Mp2 = block0.block(q+2,r,s);                     // f_i+2
              const auto& Mp3 = block0.block(q+3,r,s);                     // f_i+3

              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            } else if (q == 0) { // left halo / inside
              const auto& Mm2 = block0_left.block(block0_left.Nx-2,r,s); // f_i-2
              const auto& Mm1 = block0_left.block(block0_left.Nx-1,r,s); // f_i-1
              const auto& M   = block0.block(q,r,s);                     // f_i
              const auto& Mp1 = block0.block(q+1,r,s);                   // f_i+1
              const auto& Mp2 = block0.block(q+2,r,s);                   // f_i+2
              const auto& Mp3 = block0.block(q+3,r,s);                   // f_i+3

              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            } else if (q == 1) { // left halo / inside
              const auto& Mm2 = block0_left.block(block0_left.Nx-1,r,s); // f_i-2
              const auto& Mm1 = block0.block(q-1,r,s);                   // f_i-1
              const auto& M   = block0.block(q,r,s);                     // f_i
              const auto& Mp1 = block0.block(q+1,r,s);                   // f_i+1
              const auto& Mp2 = block0.block(q+2,r,s);                   // f_i+2
              const auto& Mp3 = block0.block(q+3,r,s);                   // f_i+3

              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            } else if ( (q >= 2) && (q <= Nx-4) ) { // inside
              const auto& Mm2 = block0.block(q-2,r,s); // f_i-2
              const auto& Mm1 = block0.block(q-1,r,s); // f_i-1
              const auto& M   = block0.block(q,r,s);   // f_i
              const auto& Mp1 = block0.block(q+1,r,s); // f_i+1
              const auto& Mp2 = block0.block(q+2,r,s); // f_i+2
              const auto& Mp3 = block0.block(q+3,r,s); // f_i+3

              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            } else if (q == Nx-3) { // right halo / inside
              const auto& Mm2 = block0.block(q-2,r,s);     // f_i-2
              const auto& Mm1 = block0.block(q-1,r,s);     // f_i-1
              const auto& M   = block0.block(q,r,s);       // f_i
              const auto& Mp1 = block0.block(q+1,r,s);     // f_i+1
              const auto& Mp2 = block0.block(q+2,r,s);     // f_i+2
              const auto& Mp3 = block0_right.block(0,r,s); // f_i+3
            
              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            } else if (q == Nx-2) { // right halo / inside
              const auto& Mm2 = block0.block(q-2,r,s);     // f_i-2
              const auto& Mm1 = block0.block(q-1,r,s);     // f_i-1
              const auto& M   = block0.block(q,r,s);       // f_i
              const auto& Mp1 = block0.block(q+1,r,s);     // f_i+1
              const auto& Mp2 = block0_right.block(0,r,s); // f_i+2
              const auto& Mp3 = block0_right.block(1,r,s); // f_i+3

              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            } else if (q >= Nx-1) { //right halo
              const auto& Mm2 = block0.block(q-2,r,s);     // f_i-2
              const auto& Mm1 = block0.block(q-1,r,s);     // f_i-1
              const auto& M   = block0.block(q,r,s);       // f_i
              const auto& Mp1 = block0_right.block(0,r,s); // f_i+1
              const auto& Mp2 = block0_right.block(1,r,s); // f_i+2
              const auto& Mp3 = block0_right.block(2,r,s); // f_i+3

              flux = flux3rdU(Mm2, Mm1, M, Mp1, Mp2, Mp3, cfl);
            }

            //const auto& M   = block0.block(q,r,s);       // f_i
            //const auto& Mp1 = block0.block(q+1,r,s);     // f_i+1


            // new local time step targets to update into (N's)
            auto& N   = block1.block(q,  r,s); // f_i  ^t+dt
            auto& Np1 = block1.block(q+1,r,s); // f_i+1^t+dt


            // now flow to neighbors; only local flows are allowed
            if(q >= 0)    N   -= flux; // - (dt/dx)U_i+1/2 (outflowing from cell)
            if(q <= Nx-2) Np1 += flux; // + (dt/dx)U_i-1/2 (inflowing to neighbor)


            // calculate current
            // NOTE: Flux (dt/dx)U is given in units of grid velocity.
            //       Since Maxwell's fields are also in the same 'units'
            //       we do not scale back to speed of light with (qm/cfl)
            T jx = qm*integrate_moment( 
                flux,
                [](std::array<T,3>&  /*unused*/) -> T { return 1.0;}
                );
              

            // vertex centered
            if(q >= 0)    gs.jx(q,r,s)   += jx; //U_i+1/2
            
            // cell centered
            //if(q >= 0)    gs.jx(q,r,s)   += jx/2.0; //U_i+1/2
            //if(q <= Nx-2) gs.jx(q+1,r,s) += jx/2.0; //U_i-1/2
            
          }
        }
      }
    }



    // Strang splitted rotating (X/2 Y X/2) solver 
    void solve( vlv::Tile<1>& tile, corgi::Grid<1>& grid ) override
    {

      // Yee lattice reference
      auto& gs = tile.get_grids();
      gs.jx.clear();
      gs.jy.clear();
      gs.jz.clear();


      // get reference to the Vlasov fluid that we are solving
      auto& step0 = tile.steps.get(0);
      auto& step1 = tile.steps.get(1);

      // timestep
      //T cfl = tile.dt/tile.dx;
      auto cfl = tile.cfl;


      // loop over different particle species (zips current [0] and new [1] solutions)
      int ispc = 0; // ith particle species
      for(auto&& blocks : zip(step0, step1) ) {

        // internal blocks of vlv meshes for current and future timestep
        auto& block0 = std::get<0>(blocks); // Solution^t
        auto& block1 = std::get<1>(blocks); // Solution^t+dt
        T qm = 1.0 / block0.qm;  // charge to mass ratio


        // external neighbors
        auto& block0_left   = tile.get_external_data(grid, ispc, -1);
        auto& block0_right  = tile.get_external_data(grid, ispc, +1);
        //auto& block0_bottom = get_external_data( 0,-1, ispc, tile, grid);
        //auto& block0_top    = get_external_data( 0,+1, ispc, tile, grid);

        // sweep in X
        xsweep(block0, block1, block0_left, block0_right, qm, cfl, gs);
        // ysweep(block0, block1, block0_bottom, block0_top,   qm, dt, dx);
        // xsweep(block0, block1, block0_left,   block0_right, qm, dt/2, dx);


        ispc++;
      }

      // done
    }


};



} // end of vlv namespace
