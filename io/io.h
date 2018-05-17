#pragma once

#include <vector>
#include <string>


#include "../em-fields/fields.h"
#include "../tools/mesh.h"

//#include "../tools/HighFive/include/highfive/H5File.hpp"
//#include "../tools/HighFive/include/highfive/H5DataSet.hpp"
//#include "../tools/HighFive/include/highfive/H5DataSpace.hpp"
#include "../tools/ezh5/src/ezh5.hpp"



namespace h5io {

using ezh5::File;

using std::string;
using std::to_string;


class Namer {

  private:
    const string extension = ".h5";


  public:
    string name;

    Namer(string& prefix) 
    {
      name = prefix + extension;
    }
      
    //standard numbering scheme
    // TODO generalize to variable argument
    string numbering(size_t i, size_t j, size_t k)
    {
      return to_string(i) + "_" + to_string(j) + "_" + to_string(k);
    }


};


class Writer {

  private:
    
    // Object to handle file names and extensions
    Namer fname;

    // actual iostream of hdf5 file
    File file;

  public:

    // constructor that creates a name and opens the file handle
    Writer(string& prefix) : 
      fname(prefix),
      file(fname.name, H5F_ACC_TRUNC) 
    {};


    // Destructor that explicitly closes the file handle
    ~Writer() {
      file.~File(); // call destructor explicitly
    }



    bool write( fields::PlasmaCell& cell )
    {
      auto& yee = cell.getYee();

      // internal cell numbering 
      string numbering = fname.numbering(cell.my_i, cell.my_j, 0);

      // open individual group for the data
      auto gr = file["yee_"+numbering];

      //Nx, Ny, Nz
      gr["Nx"] = yee.Nx;
      gr["Ny"] = yee.Ny;
      gr["Nz"] = yee.Nz;

      gr["jx"] = yee.jx.serialize();
      //gr["jy"] = yee.jy.serialize();
      //gr["jz"] = yee.jz.serialize();

      gr["ex"] = yee.ex.serialize();
      //gr["ey"] = yee.ey.serialize();
      //gr["ez"] = yee.ez.serialize();

      gr["bx"] = yee.bx.serialize();
      //gr["by"] = yee.by.serialize();
      //gr["bz"] = yee.bz.serialize();

      gr["rho"] = yee.rho.serialize();


      return true;
    }







};




} // end of namespace io
