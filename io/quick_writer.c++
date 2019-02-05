#include "quick_writer.h"

#include "../tools/mesh.h"
#include "../em-fields/tile.h"
#include "../tools/ezh5/src/ezh5.hpp"


using ezh5::File;

template<size_t D>
inline bool h5io::QuickWriter<D>::write(
    corgi::Node<D>& grid)
{

  for(auto cid : grid.get_tile_ids() ){
    auto& tile = dynamic_cast<fields::Tile<D>&>(grid.get_tile( cid ));

  }


  return true;
}



//--------------------------------------------------
// explicit template member instantiation
//template bool QuickWriter::write(const corgi::Node<1>& grid);
//template bool QuickWriter::write(const corgi::Node<2>& grid);
//template bool QuickWriter::write(const corgi::Node<3>& grid);


template class h5io::QuickWriter<1>;
template class h5io::QuickWriter<2>;
