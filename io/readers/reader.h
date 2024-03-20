#pragma once

#include <string>
#include <utility>
#include <vector>

#include "../../emf/tile.h"
#include "../../pic/tile.h"
#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../vlasov/tile.h"

#include "../namer.h"


namespace h5io {

//using ezh5::File;

class Reader 
{

  public:

    /// actual iostream of hdf5 file
    //File file;

  public:

  /// Object to handle file names and extensions
  Namer fname;

  // folder to read
  //std::string folder;

  // file name
  //std::string fname;

  // lap to consider
  //int lap;

  // my rank
  //int my_rank;

  public:

  //Reader(std::string  folder, int lap, int rank) : 
  //  folder{std::move(folder)},
  //  lap(lap),
  //  my_rank(rank)
  //{};

  //Reader(string& prefix, int lap) : 
  //    fname(prefix, lap),
  //    file(fname.name, H5F_ACC_RDONLY) 
  //{};

  Reader(string& prefix, int lap) : 
      fname(prefix, lap)
  {};

  /// read directory and locate what rank holds the given tile
  // FIXME: implement this
  template<size_t D>
  int locate_rank(const corgi::Tile<D>& /*tile*/) const
  {
    // c++-17 solution
    //#include <filesystem>
    //namespace fs = std::filesystem;//std::string path = "/path/to/directory";
    //for (auto & p : fs::directory_iterator(path))
    //    std::cout << p << std::endl;

    // assuming no tile migration
    //return my_rank;
    return 0;
  }


  template<size_t D>
  bool read(emf::Tile<D>& tile, ezh5::File& file);

  template<size_t D>
  bool read(   vlv::Tile<D>& tile, ezh5::File& file);

  template<size_t D>
  bool read(   pic::Tile<D>& tile, ezh5::File& file);

};


} // ns h5io


//--------------------------------------------------
// template implementations

#include "fields.h"
#include "vlasov.h"
#include "pic.h"


