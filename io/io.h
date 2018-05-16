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


class Writer {

  private:
    const string extension = ".h5";


  public:

    bool write( fields::PlasmaCell& cell )
    {

      // file naming
      string prefix ("fields_");
      string numbering = 
                to_string(cell.my_i) 
        + "_" + to_string(cell.my_j) 
        + "-" + to_string(cell.owner);
      string name = prefix + extension;

      // open
      File file (name, H5F_ACC_TRUNC);


      //write
      auto& yee = cell.getYee();
      write(file, yee);


      return true;
    }




    bool write(
        File& file, 
        fields::YeeLattice& yee
        )
    {
      auto gr = file["Yee"];

      //Nx, Ny, Nz
      gr["Nx"] = yee.Nx;
      gr["Ny"] = yee.Ny;
      gr["Nz"] = yee.Nz;

      //ex,ey,ez
      //bx,by,bz
      //jx,jy,jz
      //rho



      return true;
    }







};




} // end of namespace io
