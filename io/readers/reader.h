#pragma once

#include "../../em-fields/tile.h"
#include "../../pic/tile.h"
#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../vlasov/tile.h"

#include <string>
#include <utility>
#include <vector>

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

  Reader(std::string  folder, int lap, int rank) : 
    folder{std::move(folder)},
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

#include "fields.h"
#include "vlasov.h"
#include "pic.h"


