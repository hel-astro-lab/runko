#pragma once

#include "array"

#include "../definitions.h"
#include "../corgi/tile.h"

#include "particle.h"

#include "../em-fields/tile.h"

namespace pic {


/*! \brief PiC tile
 *
 * Tile infrastructures are inherited from corgi::Tile
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


  //ParticleBlock container;
  std::vector<ParticleBlock> containers;

  /// get i:th container
  ParticleBlock& get_container(size_t i) { return containers[i]; };

  /// set i:th container
  void set_container(const ParticleBlock& block) {containers.push_back(block);};

  size_t Nspecies() { return containers.size(); };


  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz)
  { }


  /// destructor
  ~Tile() override = default;

  /// tile temporal and spatial scales

  using fields::Tile<D>::cfl;
  using fields::Tile<D>::dx;

};



} // end of namespace pic
