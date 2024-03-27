#pragma once

#include "vlv/momentum-solvers/amr_momentum_solver.h"


namespace vlv {

/// \brief Forward semi-Lagrangian adaptive advection solver
template<typename T, size_t D>
class AmrMomentumFwdLagrangianSolver : public MomentumSolver<T, D> {

  public:

    typedef std::array<T, 3> vec;

    virtual void solve_mesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& Einc,
        vec& Binc,
        T qm)
    {

      // empty the target mesh
      // TODO: is this efficient; should we recycle instead?
      mesh1.data.clear();

      // create vectors
      Vec3E B(Binc);  
      Vec3E E(Einc);  

      Vec3E Fhalf = E/2.0; // construct (1/2) force

      // fmt::print("F: {} {} {}\n", Fhalf(0), Fhalf(1), Fhalf(2));

      T val = T(0);
      for(auto&& cid : mesh0.get_tiles(false) ) {
        if(! mesh0.is_leaf(cid)) continue;

        auto index = mesh0.get_indices(cid);
        int rfl    = mesh0.get_refinement_level(cid);
        auto len   = mesh0.get_size(rfl);

        vec uvel   = mesh0.get_center(index, rfl);  // velocity
        vec du     = mesh0.get_length(rfl);         // box size, i.e., \Delta u


        // get shift of the characteristic solution
        Vec3E P = qm*(Fhalf + Fhalf);


        // shift in units of tile (i.e., CFL number)
        int CFLr = (int) floor(P(0)*(dt/du[0])),
            CFLs = (int) floor(P(1)*(dt/du[1])),
            CFLt = (int) floor(P(2)*(dt/du[2]));

        // compute how many grid indices did we advect
        // int r1 = index[0] - CFLr,
        //     s1 = index[1] - CFLs,
        //     t1 = index[2] - CFLt;

        std::array<uint64_t, 3> index1 = 
        {{
          (uint64_t) index[0] - CFLr,
          (uint64_t) index[1] - CFLs,
          (uint64_t) index[2] - CFLt
        }};


        // set boundary conditions (zero outside the grid limits)
        if(    index1[0] < 2  
            || index1[1] < 2       
            || index1[2] < 2 
            || index1[0] > len[0]-1 
            || index1[1] > len[1]-1 
            || index1[2] > len[2]-1 ) 
        { 
          val = T(0); 
        } else {
          // interpolation branch
            

          // internal shift in units of tile's length
          std::array<T, 3> deltau = {{
            P(0)*(dt/du[0]) - (T) CFLr,
            P(1)*(dt/du[1]) - (T) CFLs,
            P(2)*(dt/du[2]) - (T) CFLt
          }};

          /*
          fmt::print("index: ({},{},{}); shift: ({},{},{}); newind ({},{},{}) dv: ({},{},{}) \n",
            index[0], index[1], index[2],
            CFLr, CFLs, CFLt,
            r1, s1, t1,
            deltau[0], deltau[1], deltau[2]);
           */


          //val = trilinear_interp(mesh0, index, deltau, rfl);
          val = tricubic_interp(mesh0, index, deltau, rfl);
        }

        uint64_t cid1 = mesh0.get_tile_from_indices(index1, rfl);

        mesh1.set_recursively(cid1, val);
      }

      return;
    }

};

} // end of namespace vlv

