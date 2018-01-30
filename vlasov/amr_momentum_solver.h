#pragma once

#include <cmath> 

#include <Eigen/Dense>

#include "cell.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"
#include "amr/numerics.h"
#include "../tools/cppitertools/zip.hpp"
#include "../tools/cppitertools/enumerate.hpp"


using namespace Eigen;

using std::floor;
using iter::zip;


namespace vlasov {



/*! General interface for AMR Momentum space solvers
 *
 * Solve is the general function that feeds stuff to the actual solveMesh()
 * that does all the dirty work. Hence, new solvers should only introduce
 * solveMesh() functions and let solve do what ever it does to loop the meshes
 * correctly from the container.
 *
 */
template<typename T>
class MomentumSolver {


  public:
    typedef std::array<T, 3> vec;


    /*! \brief Solve Vlasov cell contents
     *
     * Exposes the actual momentum mesh from the Vlasov Cell containers
     * and feeds those to the mesh solver.
     */
    void solve( vlasov::VlasovCell& cell ) 
    {
        
      // get reference to the Vlasov fluid that we are solving
      auto& step0 = cell.steps.get(0);
      auto& step1 = cell.steps.get(1);


      // get reference to the Yee grid 
      auto& yee = cell.getYee();

      // timestep
      T dt = (T) cell.dt;


      // loop over different particle species (zips current [0] and new [1] solutions)
      for(auto&& blocks : zip(step0, step1) ) {
          
        // loop over the cell's internal grid
        auto& block0 = std::get<0>(blocks);
        auto& block1 = std::get<1>(blocks);
          
        for(size_t q=0; q<block0.Nx; q++) {
          for(size_t r=0; r<block0.Ny; r++) {
            for (size_t s=0; s<block0.Nz; s++) {
              T qm = 1.0 / block0.qm;  // charge to mass ratio

              // Get local field components
              vec 
                B = 
                {{
                   (T) yee.bx(q,r,s),
                   (T) yee.by(q,r,s),
                   (T) yee.bz(q,r,s)
                 }},               
                                   
                E =                
                {{                 
                   (T) yee.ex(q,r,s),
                   (T) yee.ey(q,r,s),
                   (T) yee.ez(q,r,s) 
                 }};

              // dig out velomeshes from blocks
              auto& mesh0 = block0.block(q,r,s);
              auto& mesh1 = block1.block(q,r,s);

              // then the final call to the actual mesh solver
              solveMesh( mesh0, mesh1, E, B, qm, dt);
            }
          }
        }
      }

      return;
    }


    virtual void solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& E,
        vec& B,
        T qm,
        T dt) = 0;

};






/// \brief Trilinear semi-Lagrangian adaptive advection solver
template<typename T>
class AmrMomentumLagrangianSolver : public MomentumSolver<T> {

  public:

    typedef std::array<T, 3> vec;

    virtual void solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& Einc,
        vec& Binc,
        T qm,
        T dt)
    {

      // empty the target mesh
      // TODO: is this efficient; should we recycle instead?
      mesh1.data.clear();

      // create vectors
      Vector3f B(Binc.data());  
      Vector3f E(Einc.data());  

      Vector3f Fhalf = E/2.0; // construct (1/2) force

      fmt::print("F: {} {} {}\n", Fhalf(0), Fhalf(1), Fhalf(2));

      T val = T(0);
      for(auto&& cid : mesh0.get_cells(false) ) {
        if(! mesh0.is_leaf(cid)) continue;

        auto index = mesh0.get_indices(cid);
        int rfl    = mesh0.get_refinement_level(cid);
        auto len   = mesh0.get_size(rfl);

        vec uvel   = mesh0.get_center(index, rfl);  // velocity
        vec du     = mesh0.get_length(rfl);         // box size, i.e., \Delta u


        // get shift of the characteristic solution
        Vector3f P = qm*(Fhalf + Fhalf);


        // shift in units of cell (i.e., CFL number)
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
            

          // internal shift in units of cell's length
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

        uint64_t cid1 = mesh0.get_cell_from_indices(index1, rfl);

        mesh1.set_recursively(cid1, val);
      }



      return;
    }

};




} // end of namespace vlasov
