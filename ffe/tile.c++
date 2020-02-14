#include <iostream>

#include "tile.h"


namespace ffe {
  using namespace mpi4cpp;



//--------------------------------------------------
// explicit template instantiation
template class Tile<2>;
template class Tile<3>;

} // end of ns ffe
