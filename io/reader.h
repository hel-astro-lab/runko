#pragma once

#include <vector>
#include <string>

#include "../tools/ezh5/src/ezh5.hpp"

#include "../em-fields/tile.h"
#include "../vlasov/tile.h"
#include "../pic/tile.h"

namespace h5io {

using ezh5::File;


class Reader 
{

  private:

  // folder to read
  std::string folder;

  // lap to consider
  int lap;

  // my rank
  int my_rank;

  public:

  Reader(const std::string& folder, int lap, int rank) : 
    folder(folder),
    lap(lap),
    my_rank(rank)
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
    return my_rank;
  }


  template<size_t D>
  bool read(fields::Tile<D>& tile);

  template<size_t D>
  bool read(vlv::Tile<D>& tile);

  template<size_t D>
  bool read(pic::Tile<D>& tile);

};


} // ns h5io


//--------------------------------------------------
// template implementations

#include "input_fields.h"
#include "input_vlasov.h"
#include "input_pic.h"


