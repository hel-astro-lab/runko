#include "core/pic2/tile.h"

namespace pic2 {

template<std::size_t D>
Tile<D>::Tile(const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf2::Tile<D>(conf)
{
}

}  // namespace pic2

template class pic2::Tile<3>;
