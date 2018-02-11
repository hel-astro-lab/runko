#pragma once

#include <cmath> 

#include <Eigen/Dense>

#include "cell.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"
#include "amr/numerics.h"
#include "amr/refiner.h"
#include "../tools/cppitertools/zip.hpp"
#include "../tools/cppitertools/enumerate.hpp"
#include "../tools/signum.h"


using namespace Eigen;

using std::floor;
using iter::zip;
using toolbox::sign;

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

              // fmt::print("solving for srq ({},{},{})\n",s,r,q);

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




/// \brief Forward semi-Lagrangian adaptive advection solver
template<typename T>
class AmrMomentumFwdLagrangianSolver : public MomentumSolver<T> {

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

      // fmt::print("F: {} {} {}\n", Fhalf(0), Fhalf(1), Fhalf(2));

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



/// \brief Backward semi-Lagrangian adaptive advection solver
template<typename T>
class AmrMomentumLagrangianSolver : public MomentumSolver<T> {

  public:

    typedef std::array<T, 3> vec;


    /// Relativistic Lorentz force
    inline Vector3f lorentz_force(
        Vector3f& uvel,
        Vector3f& E,
        Vector3f& B,
        T qm,
        T dt)
    {

      Vector3f Ehalf = dt*E/2; // half step in E
      Vector3f us = uvel - Ehalf; // P^*, i.e., velocity in the middle of the step

      // B-field rotation matrix
      Vector3f b = B.normalized();
      T gamma = sqrt(1.0 + us.transpose()*us );
      T wt = dt*(qm*B.norm() / gamma); // relativistic cyclotron frequency

      Matrix3f Rot;
      Rot <<
        b(0)*b(0)+(1-b(0)*b(0))*cos(wt),    b(0)*b(1)*(1-cos(wt))-b(2)*sin(wt), b(0)*b(2)*(1-cos(wt))+b(1)*sin(wt),
        b(0)*b(1)*(1-cos(wt))+b(2)*sin(wt), b(1)*b(1)+(1-b(1)*b(1))*cos(wt),    b(1)*b(2)*(1-cos(wt))-b(0)*sin(wt),
        b(0)*b(2)*(1-cos(wt))-b(1)*sin(wt), b(1)*b(2)*(1-cos(wt))+b(0)*sin(wt), b(2)*b(2)+(1-b(2)*b(2))*cos(wt);

      return qm*( -Ehalf + Rot*us - uvel )*dt;
    }



    // backward advected Lorentz force
    T backward_advect(
        std::array<uint64_t, 3>& index,
        int rfl,
        const toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        Vector3f& E,
        Vector3f& B,
        T qm,
        T dt) 
    {
      T val; // return value
      vec u    = mesh1.get_center(index, rfl);  // velocity
      vec du   = mesh1.get_length(rfl);         // box size, i.e., \Delta u
      auto len = mesh1.get_size(rfl);


      // get shift of the characteristic solution from Lorentz force
      Vector3f uvel( u.data() );
      Vector3f F = lorentz_force(uvel, E, B, qm, dt);


      // advection in units of cells
      std::array<T,3> shift, cell_shift;
      for(size_t i=0; i<3; i++) shift[i] = F(i)/du[i];

      // advected cells
      std::array<int,3> index_shift;
      for(size_t i=0; i<3; i++) index_shift[i] = static_cast<int>( trunc(shift[i]) );

      // advection inside cell
      T tmp;
      for(size_t i=0; i<3; i++) cell_shift[i] = modf(shift[i], &tmp);

      // new grid indices
      std::array<uint64_t, 3> index_new;
      for(size_t i=0; i<3; i++) index_new[i] = index[i] + index_shift[i];


      // set boundary conditions (zero outside the grid limits)
      if(    index_new[0] < 2
          || index_new[1] < 2
          || index_new[2] < 2
          || index_new[0] >= len[0]-2 
          || index_new[1] >= len[1]-2 
          || index_new[2] >= len[2]-2 ) 
      { 
        val = T(0); 
      } else {
        // interpolation branch
        val = tricubic_interp(mesh0, index_new, cell_shift, rfl);
      }

      return val;
    }


    virtual void solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& Einc,
        vec& Binc,
        T qm,
        T dt)
    {

      toolbox::Adapter<T,3> adapter;
      adapter.cells_to_refine.clear();
      adapter.cells_to_unrefine.clear();



      // empty the target mesh
      // TODO: is this efficient; should we recycle instead?
      mesh1.data.clear();


      std::array<uint64_t, 3> index;
      // std::array<T, 3> grad;

      Vector3f B(Binc.data());  
      Vector3f E(Einc.data());  


      // level zero fill
      T val = T(0);
      T refine_indicator, unrefine_indicator;
      auto len = mesh1.get_size(0);
      for(uint64_t r=0; r<len[0]; r++) {
        index[0] = r;
        for(uint64_t s=0; s<len[1]; s++) {
          index[1] = s;
          for(uint64_t t=0; t<len[2]; t++) {
            index[2] = t;

            uint64_t cid = mesh1.get_cell_from_indices(index, 0);
            val = backward_advect(index, 0, mesh0, mesh1, E, B, qm, dt);

            // refinement
            // TODO specialize to value & gradient instead of just value
            refine_indicator   = val;
            unrefine_indicator = val;
            adapter.checkCell(mesh1, cid, refine_indicator, unrefine_indicator);

            mesh1.set(cid, val);
          }
        }
      }


      // create new leafs
      // TODO fixme
      // for(size_t sweep=1; sweep<=mesh1.maximum_refinement_level; sweep++){
      for(size_t sweep=1; sweep<=0; sweep++){
      
        adapter.refine(mesh1);

        // fmt::print(" sol: cells created {}\n", adapter.cells_to_refine.size());
        // fmt::print(" sol: cells removed {}\n", adapter.cells_removed.size());

        adapter.cells_to_refine.clear();
        adapter.cells_to_unrefine.clear();

        // next we refine
        for(auto&& cid : adapter.cells_created) {
          int rfl = mesh1.get_refinement_level(cid);
          auto index2 = mesh1.get_indices(cid);

          // fmt::print("creating {} at {}\n", cid, rfl);

          val = backward_advect(index2, rfl, mesh0, mesh1, E, B, qm, dt);
          mesh1.set(cid, val);

          refine_indicator   = val;
          unrefine_indicator = val;
          adapter.checkCell(mesh1, cid, refine_indicator, unrefine_indicator);
        }

        // now unrefine 
        adapter.unrefine(mesh1);

      }

      return;
    }

};


} // end of namespace vlasov
