#pragma once

#include <vector>
#include <string>


#include "../em-fields/tile.h"
#include "../tools/mesh.h"
#include "../vlasov/amr/mesh.h"

//#include "../tools/HighFive/include/highfive/H5File.hpp"
//#include "../tools/HighFive/include/highfive/H5DataSet.hpp"
//#include "../tools/HighFive/include/highfive/H5DataSpace.hpp"
#include "../tools/ezh5/src/ezh5.hpp"



namespace h5io {

using ezh5::File;

using std::string;
using std::to_string;


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
      
    //standard numbering scheme
    // TODO generalize to variable argument
    string numbering(size_t i, size_t j, size_t k)
    {
      return to_string(i) + "_" + to_string(j) + "_" + to_string(k);
    }

    string numbering(std::tuple<size_t, size_t, size_t>& ind)
    {
      return numbering(
          std::get<0>(ind),
          std::get<1>(ind),
          std::get<2>(ind)
          );
    }

};


/// data struct for location information about tiles
class TileInfo {

  public:
    int i, j, k;
    int q, r, s;
    int sp;

    string prefix;

};


// helper function to get variable length indices
std::tuple<size_t, size_t, size_t> expand_indices( 
    corgi::Tile<1>* tile )
{
  int my_i = static_cast<int>( std::get<0>(tile->index) );
  return std::make_tuple(my_i, 0,0 );
};

// helper function to get variable length indices
std::tuple<size_t, size_t, size_t> expand_indices( 
    corgi::Tile<2>* tile )
{
  int my_i = static_cast<int>( std::get<0>(tile->index) );
  int my_j = static_cast<int>( std::get<1>(tile->index) );
  return std::make_tuple(my_i, my_j,0 );
};

// helper function to get variable length indices
std::tuple<size_t, size_t, size_t> expand_indices( 
    corgi::Tile<3>* tile )
{
  int my_i = static_cast<int>( std::get<0>(tile->index) );
  int my_j = static_cast<int>( std::get<1>(tile->index) );
  int my_k = static_cast<int>( std::get<2>(tile->index) );
  return std::make_tuple(my_i, my_j,my_k );
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



    /// Write PlasmaTile content into a hdf5 data group
    template<size_t D>
    bool writeYee( fields::Tile<D>& tile )
    {
      auto& yee = tile.getYee();

      // internal tile numbering 
      auto my_ind = expand_indices( &tile );
      string numbering = fname.numbering(my_ind);

      // open individual group for the data
      auto gr = file["yee_"+numbering];

      // tile location inside node
      gr["i"] = std::get<0>(my_ind);
      gr["j"] = std::get<1>(my_ind);
      gr["k"] = std::get<2>(my_ind);

      // size
      gr["Nx"] = yee.Nx;
      gr["Ny"] = yee.Ny;
      gr["Nz"] = yee.Nz;

      //--------------------------------------------------
      // Yee lattice quantities

      gr["jx"] = yee.jx.serialize();
      gr["jy"] = yee.jy.serialize();
      gr["jz"] = yee.jz.serialize();

      gr["ex"] = yee.ex.serialize();
      gr["ey"] = yee.ey.serialize();
      gr["ez"] = yee.ez.serialize();

      gr["bx"] = yee.bx.serialize();
      gr["by"] = yee.by.serialize();
      gr["bz"] = yee.bz.serialize();

      gr["rho"] = yee.rho.serialize();


      return true;
    }



    /// Write PlasmaTile content into a hdf5 data group
    template<size_t D>
    bool writeAnalysis( fields::Tile<D>& tile )
    {

      int Nspecies = static_cast<int>(tile.analysis.size());

      for(int ispcs = 0; ispcs < Nspecies; ispcs++) {

        auto& analysis = tile.getAnalysis(ispcs);

        // internal tile numbering 
        auto my_ind = expand_indices( &tile );
        string numbering = fname.numbering(my_ind);

        // open individual group for the data
        auto gr = file["analysis_"+numbering+"-"+to_string(ispcs) ];

        // tile location inside node
        gr["i"] = std::get<0>(my_ind);
        gr["j"] = std::get<1>(my_ind);
        gr["k"] = std::get<2>(my_ind);
        gr["ispcs"] = ispcs;

        // size
        gr["Nx"] = analysis.Nx;
        gr["Ny"] = analysis.Ny;
        gr["Nz"] = analysis.Nz;

        gr["rho"] = analysis.rho.serialize();

        gr["mgamma"] = analysis.mgamma.serialize();

        gr["Vx"] = analysis.Vx.serialize();

        gr["Tx"] = analysis.Tx.serialize();

        gr["ekin"] = analysis.ekin.serialize();


      }

      return true;
    }


    /// Write 3D adaptive mesh into file
    template<typename T>
    bool write( 
        const toolbox::AdaptiveMesh<T, 3>& mesh,
        TileInfo tinfo
        )
    {

      //--------------------------------------------------
      // write meta info about grid and create internal directory structure
      size_t len = mesh.data.size();

      std::string name1 = "tile-"
                  +       std::to_string(tinfo.i) 
                  + "_" + std::to_string(tinfo.j) 
                  + "_" + std::to_string(tinfo.k);
      auto gr1 = file[name1];

      std::string name2 = "loc-"
                  +       std::to_string(tinfo.q) 
                  + "_" + std::to_string(tinfo.r) 
                  + "_" + std::to_string(tinfo.s);
      auto gr2 = gr1[name2];

      std::string name3 = "sp-"
                  +       std::to_string(tinfo.sp);
      auto gr = gr2[name3];

      gr["i"] = tinfo.i;
      gr["j"] = tinfo.j;
      gr["k"] = tinfo.k;

      gr["q"] = tinfo.q;
      gr["r"] = tinfo.r;
      gr["s"] = tinfo.s;

      gr["sp"] = tinfo.sp;
      

      // mesh metadata
      gr["maximum_refinement_level"] = mesh.maximum_refinement_level;
      gr["error_cid"]   = mesh.error_cid;
      gr["error_index"] = mesh.error_index;
      gr["top_refinement_level"] = mesh.top_refinement_level;

      gr["length"] = mesh.length;
      gr["mins"] = mesh.mins;
      gr["maxs"] = mesh.maxs;

      gr["number_of_tiles"] = len;

      //--------------------------------------------------
      // save actual mesh points
      std::vector<uint64_t> cids;
      std::vector<T       > vals;
      cids.reserve( len );
      vals.reserve( len );

      // serialize into 1D containers
      for(uint64_t cid : mesh.get_cells(true) ) {
        cids.push_back(cid);
        vals.push_back(mesh.data.at(cid));
      }

      gr["cids"] = cids;
      gr["vals"] = vals;


      return true;
    }




};




} // end of namespace io
