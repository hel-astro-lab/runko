#pragma once

#include <string>
#include <memory>

#include "cell.h"
#include "../corgi/corgi.h"
#include "../tools/rotator.h"

#include "amr_momentum_solver.h"


namespace vlasov {

/*! \brief Plasma grid
 *
 * Grid infrastructure methods are inherited from corgi::Node
 */
class Grid : public corgi::Node {
  public:

    typedef std::shared_ptr<vlasov::VlasovCell> CellPtr;


    // copy constructor
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }
  
    // default destructor
    ~Grid() { };

    /// simple method class extension
    std::string howl() { return std::string("Auuu!"); };

    /// Cycle data containers of each cell forward
    /*
    void cycle() {
      for (auto& it: cells) {
        CellPtr cellptr = std::dynamic_pointer_cast<vlasov::VlasovCell>( it.second );
        cellptr->vmeshes.cycle();
        cellptr->yee.cycle();
      }
    }
    */


    /// Update Yee lattice boundaries
    void updateBoundaries()
    {

      for(auto cid : getCellIds() ){
        vlasov::VlasovCell& cell = dynamic_cast<vlasov::VlasovCell& >(getCell( cid ));
        cell.updateBoundaries( *this );
      }

    }


    void stepVelocity()
    {

      #pragma omp parallel
      {
        #pragma omp single
        {


          for(auto cid : getCellIds() ){
            #pragma omp task
            {
              vlasov::AmrMomentumLagrangianSolver<Realf> vsol;
              vlasov::VlasovCell& cell = dynamic_cast<vlasov::VlasovCell& >(getCell( cid ));
              vsol.solve(cell);
            }// end of omp task
          }


        }// end of omp single
      }// end of omp parallel

    }






};






} // end of namespace vlasov

