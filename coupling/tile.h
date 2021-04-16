#pragma once

#include <vector>
#include <array>
#include <mpi4cpp/mpi.h>

#include "../definitions.h"

//#include "../corgi/tile.h"
//#include "../corgi/corgi.h"
//#include "../em-fields/tile.h"
//#include "../ffe/skinny_yee.h"

#include "../ffe/tile.h"
#include "../pic/tile.h"

#include"../tools/iter/allocator.h"


namespace cpl {

/*! \brief Force-Free electrodynamics methods
 *
 */
template<std::size_t D>
class Tile : 
  virtual public pic::Tile<D>,
  virtual public ffe::Tile<D>,
  virtual public fields::Tile<D>, 
  virtual public corgi::Tile<D>, 
  virtual public ManagedParent
{

  public:

  // corgi
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;

  // em fields
  using fields::Tile<D>::mesh_lengths;
  using fields::Tile<D>::yee;
  using fields::Tile<D>::cfl;

  // ffe
  using ffe::Tile<D>::dF;
  using ffe::Tile<D>::Fn;

  //using ffe::Tile<D>::rk3_update;
  //using ffe::Tile<D>::copy_eb;

  // pic
  using pic::Tile<D>::containers;

  // enumerate and diferentiate functoinality modes with this switch
  int active_mode = 0; 


  /// constructor
  Tile(int nx, int ny, int nz) :
     corgi::Tile<D>(),
     fields::Tile<D>{nx,ny,nz},
     ffe::Tile<D>{nx,ny,nz},
     pic::Tile<D>{nx,ny,nz}
  { }

  /// destructor
  ~Tile() override = default;

  //--------------------------------------------------
  // MPI send/recv
  //std::vector<mpi::request> send_data(mpi::communicator& /*comm*/, int dest, int mode, int tag) override;
  //std::vector<mpi::request> recv_data(mpi::communicator& /*comm*/, int orig, int mode, int tag) override;

};


} // end of namespace cpl
