#pragma once

#include <string>
#include <memory>

#include "../corgi/corgi.h"
#include "../tools/rotator.h"



namespace vlv {

/*! \brief Plasma grid
 *
 * Grid infrastructure methods are inherited from corgi::Node
 */

template<std::size_t D>
class Grid : 
  virtual public corgi::Node<D> 
{
  public:


  // Standard construction with grid size
  template< typename... Dims,
    typename = corgi::internals::enable_if_t< (sizeof...(Dims) == D) && 
               corgi::internals::are_integral<Dims...>::value, void >
  > 
  Grid(Dims... dims) :
    corgi::Node<D>(dims...)
  {}


  // default destructor
  ~Grid() = default;

  /// simple method class extension test
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



} // end of namespace vlv
