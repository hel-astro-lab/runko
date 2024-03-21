#pragma once

#include "../definitions.h"
#include "../corgi/tile.h"
#include "../emf/tile.h"
#include "../pic/tile.h"
#include "photon.h"



namespace qed {

/*! \brief Radiation tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from emf::Tile
 * Particle infrastructure comes from pic::Tile.
 *
 * In addition to these, we have photon containers that we 
 * call buckets (to make a distinction between particle containers).
*/

template<std::size_t D>
class Tile :
  virtual public pic::Tile<D>,
  virtual public emf::Tile<D>, 
  virtual public corgi::Tile<D> 
{

public:

  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;
  using emf::Tile<D>::mesh_lengths;

  using emf::Tile<D>::grids;
  using emf::Tile<D>::cfl;

  using pic::Tile<D>::containers;

  /// constructor
  Tile(int nx, int ny, int nz) :
     corgi::Tile<D>(),
    emf::Tile<D>{nx,ny,nz},
       pic::Tile<D>{nx,ny,nz}
  { }


  //PhotonContainer container/storage
  std::vector<PhotonContainer> buckets;

  /// get i:th bucket
  PhotonContainer& get_bucket(int i) { return buckets[i]; };

  /// push_back new bucket
  void push_back(const PhotonContainer& bucket) {buckets.push_back(bucket);};


}; // end of tile class


} // ns qed
