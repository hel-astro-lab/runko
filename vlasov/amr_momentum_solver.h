#pragma once

#include <cmath> 


#include "cell.h"
#include "../em-fields/fields.h"

#include "amr_mesh.h"
#include "amr_numerics.h"


namespace vlasov {




template<typename T>
class MomentumSolver {


  public:
    typedef std::array<T, 3> vec;

    virtual void solveMesh( toolbox::AdaptiveMesh<T, 3>& mesh ) = 0;



    /*! \brief Solve Vlasov cell contents
     *
     * Exposes the actual momentum mesh from the Vlasov Cell containers
     * and feeds those to the actual solver.
     */
    virtual void solve( vlasov::VlasovCell& cell ) 
    {
        
      // get reference to the Vlasov fluid that we are solving
      VlasovFluid& gr          = cell.getPlasmaGrid();

      // get reference to the Yee grid 
      fields::YeeLattice& yee = cell.getYee();

      // loop over the cell's internal grid
      for(size_t q=0; q<gr.Nx; q++) {
        for(size_t r=0; r<gr.Ny; r++) {

          // loop over species
          size_t ispcs = 0;
          for(auto&& spcs : gr.species() ) {

            // mass to charge ratio
            T qm = 1.0 / gr.getQ(ispcs);

            // Get local field components
            vec 
              B = 
            {{
              yee.bx(q,r,0),
              yee.by(q,r,0),
              yee.bz(q,r,0)
            }},
              
              E = 
            {{
              yee.ex(q,r,0),
              yee.ey(q,r,0),
              yee.ez(q,r,0) 
            }};


            // and finally the mesh we are working on
            toolbox::AdaptiveMesh<T, 3>& mesh = spcs(q,r,0);

            // and then the final call to the actual mesh solver
            solveMesh( mesh, E, B, qm);

            ispcs++;
          }
        }
      }

      return;
    }


};



/// \brief Trilinear semi-Lagrangian adaptive advection solver
template<typename T>
class MomentumLagrangianSolver : public MomentumSolver<T> {

  public:

    typedef std::array<T, 3> vec;

    void solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh,
        vec& E,
        vec& B,
        T qm)
    {


      for(auto&& cell : mesh.get_cells() ) {



      }

      return;
    }
  

};




} // end of namespace vlasov
