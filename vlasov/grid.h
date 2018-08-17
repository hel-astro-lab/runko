#pragma once

#include <string>
#include <memory>

#include "../corgi/corgi.h"
#include "../tools/rotator.h"



namespace vlasov {

/*! \brief Plasma grid
 *
 * Grid infrastructure methods are inherited from corgi::Node
 */
class Grid : public corgi::Node {
  public:

    //typedef std::shared_ptr<vlasov::VlasovCell> CellPtr;


    // copy constructor
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }
  
    // default destructor
    ~Grid() = default;

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



};






} // end of namespace vlasov

