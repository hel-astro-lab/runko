#pragma once

#include <string>

#include "corgi/cell.h"
#include "velomesh.h"
#include "dataContainer.h"
#include "dataContainer.c++"


namespace vlasov {


/*! \brief Vlasov cell 
 *
 * Cell infrastructure methods are inherited from corgi::Cell
 */
class VlasovCell : public corgi::Cell {
  public:

    /// constructor
    VlasovCell(size_t i, size_t j, 
               int o, 
               size_t nx, size_t ny
               ) : corgi::Cell(i, j, o, nx, ny) { }


    /// destructor
    ~VlasovCell() { };


    // Purely for testing class expansion
    // void bark();
    std::string bark() { return std::string("Wuff!"); };

    /// Simulation data container
    datarotators::DataContainer<vmesh::VeloMesh> data;

    /// Add data to the container
    void addData(vmesh::VeloMesh m) {
      data.push_back(m);
    }

    // XXX defined only for python API
    vmesh::VeloMesh getData() {
      return *data.get();
    };

};







} // end of namespace vlasov
