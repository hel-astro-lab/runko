#include "core/emf2/tile.h"

#include <iostream>

namespace emf2 {

template<std::size_t D>
Tile<D>::Tile(const toolbox::ConfigParser& p) : corgi::Tile<D>()
{
}

}  // namespace emf2

template class emf2::Tile<3>;
