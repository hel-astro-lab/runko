#include <iostream>

#include "tile.h"


namespace ffe {
  using namespace mpi4cpp;


/// Get current integration time snapshot of Yee lattice
template<std::size_t D>
fields::YeeLattice& Tile<D>::get_yee(size_t i) 
{
  return this->yee;
}




//--------------------------------------------------
// explicit template instantiation
template class Tile<2>;
//template class Tile<3>;

} // end of ns ffe
