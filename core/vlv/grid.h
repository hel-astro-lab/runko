#pragma once

#include <string>
#include <memory>
#include <concepts>

#include "corgi/corgi.h"
#include "tools/deprecated/rotator.h"



namespace vlv {

/*! \brief Plasma grid
 *
 * Grid infrastructure methods are inherited from corgi::Grid
 */

template<std::size_t D>
class Grid : 
  public corgi::Grid<D> 
{
  public:


  // Standard construction with grid size
  template<typename... Dims>
    requires (sizeof...(Dims) == D) and (std::integral<Dims> and ...)
  Grid(Dims... dims) :
    corgi::Grid<D>(dims...)
  {}


  /// simple method class extension test
  std::string howl() { return std::string("Auuu!"); };


  /// Cycle data containers of each tile forward
  /*
     void cycle() {
     for (auto& it: tiles) {
     Tileptr tileptr = std::dynamic_pointer_cast<vlv::VlasovTile>( it.second );
     tileptr->vmeshes.cycle();
     tileptr->yee.cycle();
     }
     }
     */

  /*
  void add_tile2(
    std::shared_ptr< vlv::Tile<D> > tileptr,
    corgi::internals::tuple_of<D, size_t> indices
    )
  {
    std::cout << "adding tile (NEW)\n";
    uint64_t cid = id( indices );
    std::cout << "1 \n";
    corgi::Grid<D>::tiles.erase(cid);
    std::cout << "2 \n";
    tileptr->index               = indices;
    tileptr->cid                 = cid;
    tileptr->communication.owner = corgi::Grid<D>::rank;
    tileptr->communication.local = true; //TODO Catch error if tile is not already mine?
    std::cout << "3 \n";
    corgi::Grid<D>::tiles[cid] = tileptr;
    std::cout << "4 \n";
  }
  */



};



} // end of namespace vlv
