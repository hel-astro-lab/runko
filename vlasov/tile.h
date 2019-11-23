#pragma once

#include <string>
#include <array>
#include <vector>


#include "../definitions.h"
#include "../corgi/tile.h"
#include "amr/mesh.h"
#include "amr/integrate.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../em-fields/tile.h"

//#include "grid.h"


namespace vlv {


/*! \brief Block of Vlasov fluid's inside the tile
*
* Container to hold a plasma species block
*/
class PlasmaBlock {
  typedef toolbox::Mesh< toolbox::AdaptiveMesh<Realf,3>, 3> T;

  public:

  int Nx;
  int Ny;
  int Nz;

  T block;

  PlasmaBlock(int Nx, int Ny, int Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    block(Nx, Ny, Nz) 
  { }

  //virtual ~PlasmaBlock() = default;

  Realf qm; 
};



/*! \brief Vlasov tile 
*
* Tile infrastructure methods are inherited from corgi::Tile
* Maxwell field equation solver is inherited from fields::Tile
*/

template<std::size_t D>
class Tile : 
  virtual public fields::Tile<D>
{

  public:
  
  /// Size of the internal grid
  using fields::Tile<D>::mesh_lengths;

  int Nspecies = 2;

  int Nsteps   = 2;

  /// Simulation data container (2 steps long)
  // 
  // This container has multiple snapshots of the simulation such that:
  //  - inside is a container holding various different particle species
  //  - inside these is a PlasmaBlock that has the actual velocity grids
  //    stored in a local block
  toolbox::Rotator< std::vector<PlasmaBlock>, 2 > steps;

  /// temporary current
  toolbox::Mesh<real_short, 3> jx1;
  toolbox::Mesh<real_short, 3> jy1;
  toolbox::Mesh<real_short, 3> jz1;


  /// constructor
  Tile(int nx, int ny, int nz) :
    fields::Tile<D>(nx,ny,nz)
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


  /// explicitly avoid copies; TODO: is this needed?
  Tile(Tile& ) = delete;

  /// tile temporal and spatial scales
  using fields::Tile<D>::cfl;
  //using fields::Tile<D>::dx;


  /// General clipping threshold
  Realf threshold = 1.0e-5;


  /// Clip all the meshes inside tile
  void clip() {
    auto& species = steps.get();
    Realf norm0, norm1;
    size_t nclips;

    for(auto&& internal_mesh : species) {

      for (int k=0; k<mesh_lengths[2]; k++)
      for (int j=0; j<mesh_lengths[1]; j++)
      for (int i=0; i<mesh_lengths[0]; i++) {
        auto& mesh = internal_mesh.block(i,j,k);
        if( mesh.size() < 10 ) continue; // do not clip small meshes

        // normalize
        norm0 = integrate_moment(mesh);
        norm0 = norm0 <= 0 ? 1 : norm0;

        // actual clipping
        nclips = mesh.clip_cells(threshold);
        if (nclips == 0) continue;

        // normalize back to original weight
        // simulates collisions/diffusion
        norm1 = integrate_moment(mesh);
        norm1 = norm1 <= 0 ? 1 : norm1;

        mesh *= norm0/norm1;
      }
    }
  }
    
  /// Clip all the meshes inside tile with less aggressive clip_neighbors
  void clip_neighbors() {
    Realf norm0, norm1;
    size_t nclips;

    auto& species = steps.get();
    for(auto&& internal_mesh : species) {
      for (int k=0; k<mesh_lengths[2]; k++)
      for (int j=0; j<mesh_lengths[1]; j++)
      for (int i=0; i<mesh_lengths[0]; i++) {
        auto& mesh = internal_mesh.block(i,j,k);
        if( mesh.size() < 10 ) continue; // do not clip small meshes

        // normalize
        norm0 = integrate_moment( mesh );
        norm0 = norm0 <= 0 ? 1 : norm0;

        nclips = mesh.clip_neighbors(threshold);
        if (nclips == 0) continue;

        // normalize back to original weight
        // simulates collisions/diffusion
        norm1 = integrate_moment( mesh );
        norm1 = norm1 <= 0 ? 1 : norm1;

        mesh *= norm0/norm1;
      }
    }
  }


  /// Cycle internal plasma container to another solution step
  void cycle() { steps.cycle(); }


  /// advance location 
  virtual void step_location(corgi::Grid<D>& grid);
  /*
  {
    vlv::AmrSpatialLagrangianSolver<Realf> ssol;
    ssol.solve(*this, grid);
  }
  */

  /// get neighboring tile from grid
  // TODO: separate into own communication module/header
  PlasmaBlock& get_external_data(corgi::Grid<D>&, int, int);

  void cycle_current();

};



} // end of namespace vlv
