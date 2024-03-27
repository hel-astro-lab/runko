#pragma once

#include <vector>
#include <string>

#include "io/namer.h"

#include "tools/mesh.h"
#include "emf/tile.h"
#include "vlv/tile.h"
#include "vlv/amr/mesh.h"
#include "pic/tile.h"
#include "external/ezh5/src/ezh5.hpp"


namespace h5io {

//using ezh5::File;

/// data struct for location information about tiles
class TileInfo {
  public:
    int i, j, k;
    int q, r, s;
    int sp;

    string prefix;
};



class Writer {

  private:
    /// actual iostream of hdf5 file
    //File file;

  public:

    /// Object to handle file names and extensions
    Namer fname;

    /// constructor that creates a name and opens the file handle
    //Writer(string& prefix) : 
    //  fname(prefix),
    //  file(fname.name, H5F_ACC_TRUNC) 
    //{};

    //Writer(string& prefix, int lap) : 
    //  fname(prefix, lap),
    //  file(fname.name, H5F_ACC_TRUNC) 
    //{};

    Writer(string& prefix, int lap) : 
      fname(prefix, lap)
    {};


    /// Destructor that explicitly closes the file handle
    //~Writer() {
    //  file.~File(); // call destructor explicitly
    //}
    
    //void close() {
    //    file.~File();
    //}


    // FIXME: use SFINAE to pick the right specialization without
    //        explicitly duplicating function signatures here.
    //
    // try1:
    //
    // template< template<typename> class T>
    // bool write(const T<size_t>& tile) = delete;
    //
    // try 2:
    //
    // template<typename T, typename W>
    // bool write( const T& tile, W iomode = Write_mode::standard )
    // {
    //   return write(tile, iomode);
    // }

    template<size_t D>
    bool write(const emf::Tile<D>& tile, ezh5::File& file);

    //template<size_t D>
    //bool write2(const emf::Tile<D>& tile);

    template<size_t D>
    bool write(const vlv::Tile<D>& tile, ezh5::File& file);

    template<size_t D>
    bool write(const pic::Tile<D>& tile, ezh5::File& file);


};

} // end of namespace io

//--------------------------------------------------
// template implementations

#include "fields.h"
#include "vlv.h"
#include "pic.h"

