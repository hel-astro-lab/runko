#pragma once

#include <string>
#include <array>
#include <vector>


#include "../definitions.h"
#include "../corgi/cell.h"
#include "amr/mesh.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../em-fields/fields.h"

#include "grid.h"


namespace vlasov {


/*! \brief Block of Vlasov fluid's inside the cell
*
* Container to hold a plasma species block
*/
class PlasmaBlock {
  typedef toolbox::Mesh< toolbox::AdaptiveMesh<Realf,3>, 3> T;

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  T block;

  PlasmaBlock(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    block(Nx, Ny, Nz) 
  { }


  Realf qm; 
};



/*! \brief Vlasov cell 
*
* Cell infrastructure methods are inherited from corgi::Cell
* Maxwell field equation solver is inherited from fields::PlasmaCell
*/
class VlasovCell : 
                 virtual public fields::PlasmaCell, 
                 virtual public corgi::Cell {

public:
  
  /// Size of the internal grid
  size_t NxGrid;
  size_t NyGrid;
  size_t NzGrid;

  size_t Nspecies = 2;
  size_t Nsteps   = 2;


  /// Simulation data container (2 steps long)
  // 
  // This container has multiple snapshots of the simulation such that:
  //  - inside is a container holding various different particle species
  //  - inside these is a PlasmaBlock that has the actual velocity grids
  //    stored in a local block
  toolbox::Rotator< std::vector<PlasmaBlock>, 2 > steps;


  /// constructor
  VlasovCell(size_t i, size_t j, 
             int o, 
             size_t NxG, size_t NyG,
             size_t NxMesh, size_t NyMesh
             ) : 
    corgi::Cell(i, j, o, NxG, NyG),
    fields::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh, 1),
    NxGrid(NxMesh),
    NyGrid(NyMesh),
    NzGrid(1)
    { 

      // XXX pre-initialize; now done via python
      /*
      for(size_t t=0; t<Nsteps; t++) {
        std::vector< PlasmaBlock > particles;

        fmt::print("pushing step {}\n",t);
        for(size_t p=0; p<Nspecies; p++) {
          fmt::print("pushing species {}\n",p);
          PlasmaBlock bl(NxGrid, NyGrid, 1);
          particles.push_back( bl );
        }
        steps.push_back( particles );
      }
      */

    }


  /// destructor
  ~VlasovCell() { };

  /// cell temporal and spatial scales
  //Realf dt = 0.0;
  //Realf dx = 0.0;
  //Realf dy = 0.0;
  //Realf dz = 0.0;
  Realf cfl = 0.45;


  /// General clipping threshold
  Realf threshold = 1.0e-5;


  /// Clip all the meshes inside cell
  void clip() {
    auto& species = steps.get();

    for(auto&& internal_mesh : species) {
      for (size_t k=0; k<NzGrid; k++) {
        for (size_t j=0; j<NyGrid; j++) {
          for (size_t i=0; i<NxGrid; i++) {
            internal_mesh.block(i,j,k).clip_cells(threshold);
          }
        }
      }
    }
  }

  /// Cycle internal plasma container to another solution step
  void cycle() { steps.cycle(); }


  /// advance location 
  virtual void stepLocation(vlasov::Grid& grid);
  /*
  {
    vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
    ssol.solve(*this, grid);
  }
  */

  /// get neighboring cell from grid
  // TODO: separate into own communication module/header
  vlasov::PlasmaBlock& get_external_data(
      int i, int j, int ispc,
      vlasov::Grid& grid);


};







} // end of namespace vlasov
