#pragma once

#include <vector>
#include <string>

#include "namer.h"
#include "write_tags.h"

#include "../tools/mesh.h"
#include "../em-fields/tile.h"

#include "../vlasov/tile.h"
#include "../vlasov/amr/mesh.h"

//#include "../tools/HighFive/include/highfive/H5File.hpp"
//#include "../tools/HighFive/include/highfive/H5DataSet.hpp"
//#include "../tools/HighFive/include/highfive/H5DataSpace.hpp"
#include "../tools/ezh5/src/ezh5.hpp"



namespace h5io {

using ezh5::File;

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
    
    /// Object to handle file names and extensions
    Namer fname;

    /// actual iostream of hdf5 file
    File file;

  public:

    /// constructor that creates a name and opens the file handle
    Writer(string& prefix) : 
      fname(prefix),
      file(fname.name, H5F_ACC_TRUNC) 
    {};

    Writer(string& prefix, int lap) : 
      fname(prefix, lap),
      file(fname.name, H5F_ACC_TRUNC) 
    {};

    /// Destructor that explicitly closes the file handle
    ~Writer() {
      file.~File(); // call destructor explicitly
    }


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
    // bool write( const T& tile, W iomode = WriteMode::standard )
    // {
    //   return write(tile, iomode);
    // }

    template<size_t D>
    bool write(const fields::Tile<D>& tile, WriteMode::Standard = WriteMode::standard);

    template<size_t D>
    bool write(const fields::Tile<D>& tile, WriteMode::Analysis);

    template<size_t D>
    bool write(const vlv::Tile<D>& tile,    WriteMode::Standard = WriteMode::standard);



};

} // end of namespace io

//--------------------------------------------------
// template implementations

#include "fields.h"
#include "analysis.h"
#include "vlasov.h"





