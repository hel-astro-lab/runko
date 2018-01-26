#pragma once

#include <cmath> 


#include "cell.h"
#include "../em-fields/fields.h"
#include "amr/mesh.h"
#include "amr/numerics.h"
#include "../tools/cppitertools/zip.hpp"
#include "../tools/cppitertools/enumerate.hpp"


using iter::zip;


namespace vlasov {



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
      std::vector<PlasmaBlock>& step0 = cell.steps.get(0);
      auto& step1 = cell.steps.get(1);


      // get reference to the Yee grid 
      fields::YeeLattice& yee = cell.getYee();

      // loop over different particle species
      for(auto&& blocks : zip(step0, step1) ) {
          
        // loop over the cell's internal grid
        auto& block0 = std::get<0>(blocks);
        auto& block1 = std::get<1>(blocks);
          
        for(size_t q=0; q<block0.Nx; q++) {
          for(size_t r=0; r<block0.Ny; r++) {
            for (size_t s=0; s<block0.Nz; s++) {
              T mq = 1.0 / block0.qm;  // mass to charge ratio

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

              toolbox::AdaptiveMesh<Realf, 3>& mesh0 = block0.block(q,r,s);
              toolbox::AdaptiveMesh<Realf, 3>& mesh1 = block1.block(q,r,s);

              // then the final call to the actual mesh solver
              solveMesh( mesh0, mesh1, E, B, mq);
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
        T qm) = 0;

};



/// \brief Trilinear semi-Lagrangian adaptive advection solver
template<typename T>
class AmrMomentumLagrangianSolver : public MomentumSolver<T> {

  public:

    typedef std::array<T, 3> vec;

    virtual void solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& E,
        vec& B,
        T mq)
    {

      // empty the target mesh
      // TODO: is this efficient; should we recycle instead?
      mesh1.data.clear();

      for(auto&& cid : mesh0.get_cells(false) ) {
        if(! mesh0.is_leaf(cid)) continue;

        T val = mesh0.get(cid);
        val = 1.0;

        mesh1.set(cid, val);
      }

      return;
    }
  

};




} // end of namespace vlasov
