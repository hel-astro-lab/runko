#pragma once

#include <string>
#include <omp.h>

#include "../io/writer.h"
#include "../io/reader.h"
#include "../io/write_tags.h"



namespace pic {


template<size_t D>
void write_particles( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir 
    )
{

  std::string prefix = dir + "particles-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  for(auto cid : grid.get_local_tiles() ){
    const auto& tile 
      = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    writer.write(tile);
  }

}


template<size_t D>
inline void read_particles( 
    corgi::Node<D>& grid, 
    int lap,
    const std::string& dir 
    )
{
  h5io::Reader reader(dir, lap, grid.comm.rank());

  for(auto cid : grid.get_tile_ids() ){
    auto& tile 
      = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    reader.read(tile);
  }
}


}
