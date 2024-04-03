#pragma once

#include <string>

#include "io/writers/writer.h"
#include "io/readers/reader.h"


namespace vlv{

template<size_t D>
inline void write_mesh( 
    corgi::Grid<D>& grid, 
    int lap,
    std::string dir 
    )
{
  if(dir.back() != '/') dir += '/';

  std::string prefix = dir + "meshes-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  ezh5::File file(writer.fname.name, H5F_ACC_TRUNC);

  for(auto cid : grid.get_local_tiles() ){
    const auto& tile 
      = dynamic_cast<vlv::Tile<D>&>(grid.get_tile( cid ));
    writer.write(tile, file);
  }
}


template<size_t D>
inline void read_mesh( 
    corgi::Grid<D>& grid, 
    int lap,
    std::string dir 
    )
{
  if(dir.back() != '/') dir += '/';

  std::string prefix = dir + "meshes-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Reader reader(prefix, lap);

  ezh5::File file(reader.fname.name, H5F_ACC_RDONLY);

  for(auto cid : grid.get_tile_ids() ){
    auto& tile 
      = dynamic_cast<vlv::Tile<D>&>(grid.get_tile( cid ));
    reader.read(tile, file);
  }
}


}// end of namespace vlv


//--------------------------------------------------

namespace emf {

template<size_t D>
inline void write_grids( 
    corgi::Grid<D>& grid, 
    int lap,
    std::string dir
    )
{
  if(dir.back() != '/') dir += '/';

  std::string prefix = dir + "fields-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  ezh5::File file(writer.fname.name, H5F_ACC_TRUNC);

  for(auto cid : grid.get_local_tiles() ){
    const auto& tile 
      = dynamic_cast<emf::Tile<D>&>(grid.get_tile( cid ));
    writer.write(tile, file);
  }
}


//template<size_t D>
//inline void write_analysis( 
//    corgi::Grid<D>& grid, 
//    int lap,
//    const std::string& dir
//    )
//{
//
//  std::string prefix = dir + "analysis-"; 
//  prefix += std::to_string(grid.comm.rank());
//  h5io::Writer writer(prefix, lap);
//
//  for(auto cid : grid.get_local_tiles() ){
//    const auto& tile 
//      = dynamic_cast<emf::Tile<D>&>(grid.get_tile( cid ));
//    writer.write2(tile);
//  }
//}

template<size_t D>
inline void read_grids( 
    corgi::Grid<D>& grid, 
    int lap,
    std::string dir 
    )
{
  if(dir.back() != '/') dir += '/';

  std::string prefix = dir + "fields-"; 
  prefix += std::to_string(grid.comm.rank());

  h5io::Reader reader(prefix, lap);
  ezh5::File file(reader.fname.name, H5F_ACC_RDONLY);

  for(auto cid : grid.get_tile_ids() ){
    auto& tile 
      = dynamic_cast<emf::Tile<D>&>(grid.get_tile( cid ));
    reader.read(tile, file);
  }

}

} // end of ns emf


//--------------------------------------------------

namespace pic {


template<size_t D>
void write_particles( 
    corgi::Grid<D>& grid, 
    int lap,
    std::string dir 
    )
{
  if(dir.back() != '/') dir += '/';

  std::string prefix = dir + "particles-"; 
  prefix += std::to_string(grid.comm.rank());
  h5io::Writer writer(prefix, lap);

  ezh5::File file(writer.fname.name, H5F_ACC_TRUNC);

  for(auto cid : grid.get_local_tiles() ){
    const auto& tile 
      = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    writer.write(tile, file);
  }
}


template<size_t D>
inline void read_particles( 
    corgi::Grid<D>& grid, 
    int lap,
    std::string dir 
    )
{
  if(dir.back() != '/') dir += '/';

  std::string prefix = dir + "particles-"; 
  prefix += std::to_string(grid.comm.rank());

  h5io::Reader reader(prefix, lap);
  ezh5::File file(reader.fname.name, H5F_ACC_RDONLY);

  for(auto cid : grid.get_tile_ids() ){
    auto& tile 
      = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    reader.read(tile, file);
  }
}


}
