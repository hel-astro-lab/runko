#pragma once

#include <string>
#include <array>
#include <vector>


#include "../definitions.h"
#include "../corgi/tile.h"
#include "amr/mesh.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../em-fields/tile.h"

#include "grid.h"


namespace vlv {


/*! \brief Block of Vlasov fluid's inside the tile
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



/*! \brief Vlasov tile 
*
* Tile infrastructure methods are inherited from corgi::Tile
* Maxwell field equation solver is inherited from fields::Tile
*/

template<std::size_t D>
class Tile : 
  virtual public fields::Tile<D>, 
  virtual public  corgi::Tile<D> 
{

  public:
  
  /// Size of the internal grid
  //size_t NxGrid;
  //size_t NyGrid;
  //size_t NzGrid;
  using fields::Tile<D>::mesh_lengths;

  const size_t Nspecies = 2;

  const size_t Nsteps   = 2;

  /// Simulation data container (2 steps long)
  // 
  // This container has multiple snapshots of the simulation such that:
  //  - inside is a container holding various different particle species
  //  - inside these is a PlasmaBlock that has the actual velocity grids
  //    stored in a local block
  toolbox::Rotator< std::vector<PlasmaBlock>, 2 > steps;


  /// constructor
  template< typename... Dims,
    typename = corgi::internals::enable_if_t< (sizeof...(Dims) == D) && 
               corgi::internals::are_integral<Dims...>::value, void >
  > 
  Tile(Dims... mesh_lens) :
     corgi::Tile<D>(),
    fields::Tile<D>(mesh_lens...)
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
  ~Tile() override = default;


  /// tile temporal and spatial scales
  using fields::Tile<D>::cfl;
  using fields::Tile<D>::dt;
  using fields::Tile<D>::dx;


  /// General clipping threshold
  Realf threshold = 1.0e-5;


  /// Clip all the meshes inside tile
  void clip() {
    auto& species = steps.get();

    for(auto&& internal_mesh : species) {

      for (size_t k=0; k<mesh_lengths[2]; k++)
      for (size_t j=0; j<mesh_lengths[1]; j++)
      for (size_t i=0; i<mesh_lengths[0]; i++) {
        internal_mesh.block(i,j,k).clip_cells(threshold);
      }
    }
  }

  /// Cycle internal plasma container to another solution step
  void cycle() { steps.cycle(); }


  /// advance location 
  virtual void stepLocation(Grid<D>& grid);
  /*
  {
    vlv::AmrSpatialLagrangianSolver<Realf> ssol;
    ssol.solve(*this, grid);
  }
  */

  /// get neighboring tile from grid
  // TODO: separate into own communication module/header
  PlasmaBlock& get_external_data(
      int i, int j, int ispc,
      Grid<D>& grid);


};



} // end of namespace vlv
