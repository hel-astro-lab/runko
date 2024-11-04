#pragma once

#include <array>
#include <mpi4cpp/mpi.h>

#include "definitions.h"
#include "external/corgi/tile.h"
#include "external/corgi/corgi.h"
#include "core/emf/tile.h"
#include "core/pic/particle.h"
#include "external/iter/allocator.h"


namespace pic {

using namespace mpi4cpp;

/*! \brief PiC tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from emf::Tile
*/

template<std::size_t D>
class Tile :
  virtual public emf::Tile<D>, 
  virtual public  corgi::Tile<D>,
  virtual public ManagedParent
{


public:

  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;
  using corgi::Tile<D>::cid;

  using emf::Tile<D>::mesh_lengths;

  using emf::Tile<D>::grids;

  using emf::Tile<D>::cfl;

  std::vector<ParticleContainer<D> , ManagedAlloc<ParticleContainer<D> >> containers;

  //--------------------------------------------------
  // normal container methods
     
  /// get i:th container
  ParticleContainer<D>& get_container(int i) { return containers[i]; }

  const ParticleContainer<D>& get_const_container(int i) const { return containers[i]; };

  /// set i:th container
  void set_container(ParticleContainer<D>& block) 
  { 
    // label container with tile cid; used for identifying containers when debugging
    block.cid = cid;

    // label container with exterior tile limits
    block.mins = mins;
    block.maxs = maxs;

    // add 
    containers.push_back(block); 
  };

  int Nspecies() const { return containers.size(); };


  /// constructor
  Tile(int nx, int ny, int nz) :
     corgi::Tile<D>(),
    emf::Tile<D>{nx,ny,nz}
  { }


  //--------------------------------------------------
  // MPI send
  std::vector<mpi::request> 
  send_data( mpi::communicator& /*comm*/, int dest, int mode, int tag) override;

  /// actual tag=0 send
  std::vector<mpi::request> 
  send_particle_data( mpi::communicator& /*comm*/, int dest, int tag);

  /// actual tag=1 send
  std::vector<mpi::request> 
  send_particle_extra_data( mpi::communicator& /*comm*/, int dest, int tag);


  //--------------------------------------------------
  // MPI recv
  std::vector<mpi::request> 
  recv_data(mpi::communicator& /*comm*/, int orig, int mode, int tag) override;

  /// actual tag=0 recv
  std::vector<mpi::request> 
  recv_particle_data(mpi::communicator& /*comm*/, int orig, int tag);

  /// actual tag=1 recv
  std::vector<mpi::request> 
  recv_particle_extra_data(mpi::communicator& /*comm*/, int orig, int tag);
  //--------------------------------------------------


  /// check all particle containers for particles
  // exceeding limits
  void check_outgoing_particles();

  /// delete particles from each container that are exceeding
  // the boundaries
  void delete_transferred_particles();

  /// get particles flowing into this tile
  void get_incoming_particles(corgi::Grid<D>& grid);

  /// pack all particles for MPI message
  void pack_all_particles();

  /// pack particles for MPI message
  void pack_outgoing_particles();

  /// unpack received MPI message particles
  void unpack_incoming_particles();

  /// delete all particles from each container
  void delete_all_particles();

  /// shrink to fit all internal containers
  void shrink_to_fit_all_particles();


private:
  std::size_t dim = D;
};



} // end of namespace pic
