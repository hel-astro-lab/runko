#pragma once

#include <string>

#include "../definitions.h"
#include "../corgi/cell.h"
#include "amr_mesh.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../em-fields/fields.h"


namespace vlasov {


template<typename M, typename T>
class PairPlasmaIterator {

  private:

    /// First species 
    size_t _beginning = 0;

    /// beyond last species 
    size_t _ending    = 2;

    /// internal pointer to the parent object/container (Mother object)
    M* ptr;

    /// Iterators own internal counting system
    size_t spcs = 0; 


  public:

    PairPlasmaIterator(M& rhs) : ptr(&rhs) {}
    PairPlasmaIterator(const M& rhs) : ptr(&rhs) {}

    PairPlasmaIterator(const PairPlasmaIterator& rhs) : ptr(rhs.ptr) {}


    /// Assignment
    PairPlasmaIterator& operator= (const PairPlasmaIterator& rhs) = default;

    /// iterate
    PairPlasmaIterator& operator++ () {
      ++this->spcs;
      return *this;
    }

    /// Referencing cell interiors
    T& operator *() {
      if(spcs == 0) return (T&) (ptr->electrons);
      else if(spcs == 1) return (T&) (ptr->positrons);
      else throw std::range_error("iterator goes beyond electrons (0) or positrons (1)");
    }


    /// iterate with steps
    // PairPlasmaIterator operator++ (int) {
    //   PairPlasmaIterator temp(*ptr);
    //   ++*this;
    //   return (temp);
    // }

    /// equal comparison done by comparing internal spcs value
    bool operator== (PairPlasmaIterator& rhs) const {
      return (ptr == rhs.ptr) && (spcs == rhs.spcs);
    }

    /// unequal comparison done by comparing internal spcs value
    bool operator!= (PairPlasmaIterator& rhs) const {
      return (ptr != rhs.ptr) || (spcs != rhs.spcs);
    }

    /// Returns an iterator pointing to the first element in the sequence
    PairPlasmaIterator begin() {
      PairPlasmaIterator temp(*ptr);
      temp.spcs = _beginning;
      return temp;
    }

    /// Returns an iterator pointing to the past-the-end element in the sequence
    PairPlasmaIterator end() {
      PairPlasmaIterator temp(*ptr);
      temp.spcs = _ending;
      return temp ;
    }


};




/*! \brief Vlasov fluid inside the cell
 *
 * Container to hold different plasma species.
 */
class VlasovFluid {

  private:
    typedef toolbox::Mesh<toolbox::AdaptiveMesh<Realf, 3>, 3> T;

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  T electrons;
  T positrons;

  VlasovFluid(size_t Nx, size_t Ny, size_t Nz) : Nx(Nx), Ny(Ny), Nz(Nz),
    electrons(Nx, Ny, Nz),
    positrons(Nx, Ny, Nz) { }


  PairPlasmaIterator<VlasovFluid, T> species() {
    PairPlasmaIterator<VlasovFluid, T> ret(*this);
    return ret;
  };

  std::vector<Realf> qms;

  Realf getQ(size_t i) {
    return qms[i];
  }


};



/*! \brief Vlasov cell 
 *
 * Cell infrastructure methods are inherited from corgi::Cell
 * Maxwell field equation solver is inherited from maxwell::PlasmaCell
 */
class VlasovCell : 
                   virtual public fields::PlasmaCell, 
                   virtual public corgi::Cell {

  public:
    
    /// Size of the internal grid
    size_t NxGrid;
    size_t NyGrid;
    size_t NzGrid;


    /// constructor
    VlasovCell(size_t i, size_t j, 
               int o, 
               size_t NxG, size_t NyG,
               size_t NxMesh, size_t NyMesh
               ) : 
      corgi::Cell(i, j, o, NxG, NyG),
      fields::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh),
      NxGrid(NxMesh),
      NyGrid(NyMesh),
      NzGrid(1)
      { 
        // add also initial velomeshes
        plasma.push_back( VlasovFluid(NxGrid, NyGrid, 1) );
        plasma.push_back( VlasovFluid(NxGrid, NyGrid, 1) );
      
      }

    /// destructor
    ~VlasovCell() { };

    /// cell temporal and spatial scales
    Realf dt = 0.0;
    Realf dx = 0.0;
    Realf dy = 0.0;
    Realf dz = 0.0;


    /// General clipping threshold
    Realf threshold = 1.0e-5;

    // Purely for testing class expansion
    // void bark();
    std::string bark() { return std::string("Wuff!"); };


    /// Simulation data container
    toolbox::Rotator<VlasovFluid> plasma;


    /// Add data to the container
    /*
    void addData(vmesh::VeloMesh m) {
      vmeshes.push_back(m);
    }
    */

    // XXX defined only for python API
    // vmesh::VeloMesh getData() {
    //   return *data.get();
    // };

    VlasovFluid& getPlasmaGrid() {
      return plasma.getRef();
    };

    VlasovFluid& getNewPlasmaGrid() {
      return plasma.getNewRef();
    };


    // get pointer to the data
    // TODO: should be shared_ptr explicitly to make it memory save
    /*
    vmesh::VeloMesh* getDataPtr() {
      return vmeshes.get();
    };
    */


    /// Clip all the meshes inside cell
    void clip() {
      VlasovFluid& gr = plasma.getRef();

      for (size_t k=0; k<NzGrid; k++) {
        for (size_t j=0; j<NyGrid; j++) {
          for (size_t i=0; i<NxGrid; i++) {
            gr.electrons(i,j,k).clip_cells(threshold);
            gr.positrons(i,j,k).clip_cells(threshold);
          }
        }
      }

    }


    /// Cycle internal plasma container to another solution step
    void cycle() {
      plasma.cycle();
    }



};


} // end of namespace vlasov
