#pragma once
#include <tuple>
#include <string>

#include "corgi/tile.h"


namespace h5io {

using std::string;
using std::to_string;


// helper function to get variable length indices
inline std::tuple<size_t, size_t, size_t> expand_indices( 
    const corgi::Tile<1>* tile )
{
  int my_i = static_cast<int>( std::get<0>(tile->index) );
  return std::make_tuple(my_i, 0, 0);
}

// helper function to get variable length indices
inline std::tuple<size_t, size_t, size_t> expand_indices( 
    const corgi::Tile<2>* tile )
{
  int my_i = static_cast<int>( std::get<0>(tile->index) );
  int my_j = static_cast<int>( std::get<1>(tile->index) );
  return std::make_tuple(my_i, my_j, 0);
}

// helper function to get variable length indices
inline std::tuple<size_t, size_t, size_t> expand_indices( 
    const corgi::Tile<3>* tile )
{
  int my_i = static_cast<int>( std::get<0>(tile->index) );
  int my_j = static_cast<int>( std::get<1>(tile->index) );
  int my_k = static_cast<int>( std::get<2>(tile->index) );
  return std::make_tuple(my_i, my_j,my_k );
}


//standard numbering scheme
// TODO generalize to variable argument
inline string create_numbering(size_t i, size_t j, size_t k) 
{
  return to_string(i) + "_" + to_string(j) + "_" + to_string(k);
}

inline string create_numbering(std::tuple<size_t, size_t, size_t>& ind)
{
  return create_numbering(
      std::get<0>(ind),
      std::get<1>(ind),
      std::get<2>(ind)
      );
}



/// General class to do all the text and file name mangling 
class Namer {

  private:
    const string extension = ".h5";

  public:
    string name;

    Namer(string& prefix) 
    {
      name = prefix + extension;
    }

    Namer(string& prefix, int lap) 
    {
      name = prefix + "_" + to_string(lap) + extension;
    }
};


} // end of namespace io
