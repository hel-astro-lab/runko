#pragma once

#include <string>

#include "corgi/cell.h"
#include "velomesh.h"
#include "dataContainer.h"
#include "maxwell.h"


namespace vlasov {


/*! \brief Vlasov cell 
 *
 * Cell infrastructure methods are inherited from corgi::Cell
 */
class VlasovCell : virtual public maxwell::PlasmaCell, virtual public corgi::Cell {
  public:

    /// constructor
    VlasovCell(size_t i, size_t j, 
               int o, 
               size_t NxG, size_t NyG,
               size_t NxMesh, size_t NyMesh
               ) : 
      corgi::Cell(i, j, o, NxG, NyG),
      maxwell::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh)
      { }


    /// destructor
    ~VlasovCell() { };


    // Purely for testing class expansion
    // void bark();
    std::string bark() { return std::string("Wuff!"); };

    /// Simulation data container
    datarotators::DataContainer<vmesh::VeloMesh> vmeshes;


    /// Add data to the container
    void addData(vmesh::VeloMesh m) {
      vmeshes.push_back(m);
    }

    // XXX defined only for python API
    // vmesh::VeloMesh getData() {
    //   return *data.get();
    // };

    vmesh::VeloMesh& getData() {
      return *vmeshes.get();
    };

    // get pointer to the data
    // TODO: should be shared_ptr explicitly to make it memory save
    vmesh::VeloMesh* getDataPtr() {
      return vmeshes.get();
    };


    /// Clip all the meshes inside cell
    void clip() {
      vmeshes.get()->clip();
    };

};



} // end of namespace vlasov
