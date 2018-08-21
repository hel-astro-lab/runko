#pragma once

#include "../definitions.h"
#include "../corgi/tile.h"

#include "particle.h"

#include "../em-fields/fields.h"

namespace pic {


/*! \brief PiC tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from fields::PlasmaTile
*/
class PicTile :
  virtual public fields::PlasmaTile, 
  virtual public corgi::Tile {

public:

  /// Size of the internal grid
  size_t NxGrid;
  size_t NyGrid;
  size_t NzGrid;



  //ParticleBlock container;
  std::vector<ParticleBlock> containers;

  /// get i:th container
  ParticleBlock& get_container(size_t i) { return containers[i]; };

  /// set i:th container
  void set_container(const ParticleBlock& block) {containers.push_back(block);};

  size_t Nspecies() { return containers.size(); };


  /// constructor
  PicTile(size_t i, size_t j, 
             int o, 
             size_t NxG, size_t NyG,
             size_t NxMesh, size_t NyMesh
             ) : 
    corgi::Tile(i, j, o, NxG, NyG),
    fields::PlasmaTile(i,j,o,NxG, NyG, NxMesh,NyMesh, 1),
    NxGrid(NxMesh), NyGrid(NyMesh), NzGrid(1)
  { }


  /// destructor
  ~PicTile() override = default;

  /// tile temporal and spatial scales
  using fields::PlasmaTile::cfl;
  using fields::PlasmaTile::dt;
  using fields::PlasmaTile::dx;


};



} // end of namespace pic
