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
  public corgi::Node<D> 
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
     Tileptr tileptr = std::dynamic_pointer_cast<vlasov::VlasovTile>( it.second );
     tileptr->vmeshes.cycle();
     tileptr->yee.cycle();
     }
     }
     */

  /*
  void addTile2(
    std::shared_ptr< vlv::Tile<D> > tileptr,
    corgi::internals::tuple_of<D, size_t> indices
    )
  {
    std::cout << "adding tile (NEW)\n";
    uint64_t cid = id( indices );
    std::cout << "1 \n";
    corgi::Node<D>::tiles.erase(cid);
    std::cout << "2 \n";
    tileptr->index               = indices;
    tileptr->cid                 = cid;
    tileptr->communication.owner = corgi::Node<D>::rank;
    tileptr->communication.local = true; //TODO Catch error if tile is not already mine?
    std::cout << "3 \n";
    corgi::Node<D>::tiles[cid] = tileptr;
    std::cout << "4 \n";
  }
  */



};



} // end of namespace vlv
