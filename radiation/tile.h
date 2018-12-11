#pragma once

#include "../definitions.h"
#include "../corgi/tile.h"
#include "../em-fields/tile.h"
#include "../pic/tile.h"
#include "photon.h"



namespace rad {

/*! \brief Radiation tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from fields::Tile
 * Particle infrastructure comes from pic::Tile.
 *
 * In addition to these, we have photon containers that we 
 * call buckets (to make a distinction between particle containers).
*/

template<std::size_t D>
class Tile :
  virtual public    pic::Tile<D>,
  virtual public fields::Tile<D>, 
  virtual public  corgi::Tile<D> 
{

public:

  using fields::Tile<D>::mesh_lengths;
  using fields::Tile<D>::cfl;
  using fields::Tile<D>::dx;


  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz),
       pic::Tile<D>(nx,ny,nz)
  { }

  /// destructor
  ~Tile() override = default;

  //PhotonContainer container/storage
  std::vector<PhotonContainer> buckets;

  /// get i:th bucket
  PhotonContainer& get_bucket(size_t i) { return buckets[i]; };

  /// push_back new bucket
  void push_back(const PhotonContainer& bucket) {buckets.push_back(bucket);};



}; // end of tile class


} // ns rad
