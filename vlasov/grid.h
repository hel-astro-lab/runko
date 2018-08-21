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

    //typedef std::shared_ptr<vlasov::VlasovTile> TilePtr;


    // copy constructor
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }
  
    // default destructor
    ~Grid() = default;

    /// simple method class extension
    std::string howl() { return std::string("Auuu!"); };

    /// Cycle data containers of each tile forward
    /*
    void cycle() {
      for (auto& it: tiles) {
        TilePtr tileptr = std::dynamic_pointer_cast<vlasov::VlasovTile>( it.second );
        tileptr->vmeshes.cycle();
        tileptr->yee.cycle();
      }
    }
    */



};






} // end of namespace vlasov

