#pragma once

#include "corgi/corgi.h"
#include "velomesh.h"
#include "dataContainer.h"
#include "dataContainer.c++"


namespace vlasov {


/*! \brief Vlasov cell 
 *
 * Cell infrastructure methods are inherited from corgi::Cell
 */
class VCell : public corgi::Cell {
  public:
    VCell(size_t i, size_t j, int o) : corgi::Cell(i, j, o) { }

    // Purely for testing class expansion
    // void bark();
    void bark() { fmt::print("Woof!\n"); };

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
