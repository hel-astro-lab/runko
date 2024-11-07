#include <cmath>

#include "core/pic/tile.h"
#include "core/pic/communicate.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


namespace pic {

using namespace mpi4cpp;


// create MPI tag given tile id and extra layer of differentiation
int get_tag(int tag, int extra_param)
{
  assert(extra_param <= 8); // max 8 species
  assert(tag < (pow(2,16) - 1)); // cray-mpich supports maximum of 2^22-1 tag value

  return tag + (9+extra_param)*pow(2,16);
}

// create MPI tag for extra communication given tile 
// id and extra layer of differentiation
int get_extra_tag(int tag, int extra_param)
{
  assert(extra_param <= 8); // max 8 species
  assert(tag < (pow(2,16) - 1)); // cray-mpich supports maximum of 2^22-1 tag value

  return tag + (9+8+extra_param)*pow(2,16);
}




//--------------------------------------------------

template<std::size_t D>
void Tile<D>::check_outgoing_particles()
{
  std::array<double,3> 
    tile_mins = {{0,0,0}},
    tile_maxs = {{1,1,1}};

  for(size_t i=0; i<D; i++) tile_mins[i] = corgi::Tile<D>::mins[i];
  for(size_t i=0; i<D; i++) tile_maxs[i] = corgi::Tile<D>::maxs[i];

  for(auto&& container : containers)
    container.check_outgoing_particles(tile_mins, tile_maxs);
}

template<std::size_t D>
void Tile<D>::delete_transferred_particles()
{
  for(auto&& container : containers) 
    container.delete_transferred_particles();
}

//--------------------------------------------------

template<>
void Tile<1>::get_incoming_particles(
    corgi::Grid<1>& grid)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::array<double,3> global_mins = {
    static_cast<double>( grid.get_xmin() ),
    static_cast<double>( 0.0 ),
    static_cast<double>( 0.0 )
  };

  std::array<double,3> global_maxs = {
    static_cast<double>( grid.get_xmax() ),
    static_cast<double>( 1.0 ),
    static_cast<double>( 1.0 )
  };

  // fetch incoming particles from neighbors around me
  int j = 0;
  int k = 0;
  for(int i=-1; i<=1; i++) {
        // get neighboring tile
        auto ind = this->neighs(i); 
        uint64_t cid = grid.id( std::get<0>(ind));
        Tile& external_tile = dynamic_cast<Tile&>( grid.get_tile(cid) );

        // loop over all containers
        for(int ispc=0; ispc<Nspecies(); ispc++) {
          auto& container = get_container(ispc);
          auto& neigh = external_tile.get_container(ispc);

          container.transfer_and_wrap_particles(
              neigh, {i,j,k}, global_mins, global_maxs);
        }

  }

#ifdef GPU
  nvtxRangePop();
#endif

}

template<>
void Tile<2>::get_incoming_particles(
    corgi::Grid<2>& grid)
{

  std::array<double,3> global_mins = {
    static_cast<double>( grid.get_xmin() ),
    static_cast<double>( grid.get_ymin() ),
    static_cast<double>( 0.0             )
  };

  std::array<double,3> global_maxs = {
    static_cast<double>( grid.get_xmax() ),
    static_cast<double>( grid.get_ymax() ),
    static_cast<double>( 1.0             )
  };

  // fetch incoming particles from neighbors around me
  int k = 0;
  for(int i=-1; i<=1; i++) {
    for(int j=-1; j<=1; j++) {

      // get neighboring tile
      auto ind = this->neighs(i, j); 
      uint64_t cid = grid.id( std::get<0>(ind), std::get<1>(ind) );
      Tile& external_tile = 
        dynamic_cast<Tile&>( grid.get_tile(cid) );

      // loop over all containers
        
      for(int ispc=0; ispc<Nspecies(); ispc++) {
        auto& container = get_container(ispc);
        auto& neigh = external_tile.get_container(ispc);

        container.transfer_and_wrap_particles(
            neigh, {i,j,k}, global_mins, global_maxs);
      }
    }
  }
}


template<>
void Tile<3>::get_incoming_particles(
    corgi::Grid<3>& grid)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::array<double,3> global_mins = {
    static_cast<double>( grid.get_xmin() ),
    static_cast<double>( grid.get_ymin() ),
    static_cast<double>( grid.get_zmin() )
  };

  std::array<double,3> global_maxs = {
    static_cast<double>( grid.get_xmax() ),
    static_cast<double>( grid.get_ymax() ),
    static_cast<double>( grid.get_zmax() )
  };

  // fetch incoming particles from neighbors around me
  for(int i=-1; i<=1; i++) {
    for(int j=-1; j<=1; j++) {
      for(int k=-1; k<=1; k++) {

        if( i==0 && j==0 && k==0 ) continue;
          
        // get neighboring tile
        auto ind = this->neighs(i, j, k); 
        uint64_t cid = grid.id( std::get<0>(ind), std::get<1>(ind), std::get<2>(ind) );
        Tile& external_tile = dynamic_cast<Tile&>( grid.get_tile(cid) );

        // loop over all containers
        for(int ispc=0; ispc<Nspecies(); ispc++) {
          auto& container = get_container(ispc);
          auto& neigh = external_tile.get_container(ispc);

          container.transfer_and_wrap_particles(
              neigh, {i,j,k}, global_mins, global_maxs);
        }

        }
    }
  }

#ifdef GPU
  nvtxRangePop();
#endif

}

template<std::size_t D>
std::vector<mpi::request> Tile<D>::send_data( 
    mpi::communicator& comm, 
    int dest, 
    int mode,
    int tag)
{
  if(mode == 0) return emf::Tile<D>::send_data(comm, dest, mode, tag);
  if(mode == 1) return emf::Tile<D>::send_data(comm, dest, mode, tag);
  if(mode == 2) return emf::Tile<D>::send_data(comm, dest, mode, tag);

  if(mode == 3) return send_particle_data(comm,dest,tag);
  if(mode == 4) return send_particle_extra_data(comm,dest,tag);

  assert(false);
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::send_particle_data( 
    mpi::communicator& comm, 
    int dest,
    int tag)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::vector<mpi::request> reqs;
  for(int ispc=0; ispc<Nspecies(); ispc++) {
    auto& container = get_container(ispc);

    reqs.emplace_back(
        comm.isend(dest, get_tag(tag, ispc), 
          container.outgoing_particles.data(), 
          container.first_message_size)
        );
  }

#ifdef GPU
  nvtxRangePop();
#endif

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::send_particle_extra_data( 
    mpi::communicator& comm, 
    int dest,
    int tag)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::vector<mpi::request> reqs;
  for(int ispc=0; ispc<Nspecies(); ispc++) {
    auto& container = get_container(ispc);

    if(!container.outgoing_extra_particles.empty()) {
      reqs.emplace_back(
          comm.isend(dest, get_extra_tag(tag, ispc), 
            container.outgoing_extra_particles.data(), 
            container.outgoing_extra_particles.size())
          );
    }

    //std::cout << this->communication.cid 
    //  << " send " << container.outgoing_particles.size() 
    //  << " / " << container.outgoing_particles[0].id
    //  << " + " << container.outgoing_extra_particles.size() 
    //  << " / " << container.outgoing_extra_particles[0].id
    //  << " particles\n";
  }

#ifdef GPU
  nvtxRangePop();
#endif

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_data( 
    mpi::communicator& comm, 
    int orig, 
    int mode,
    int tag)
{
  if(mode == 0) return emf::Tile<D>::recv_data(comm, orig, mode, tag);
  if(mode == 1) return emf::Tile<D>::recv_data(comm, orig, mode, tag);
  if(mode == 2) return emf::Tile<D>::recv_data(comm, orig, mode, tag);

  if(mode == 3) return recv_particle_data(comm,orig,tag);
  if(mode == 4) return recv_particle_extra_data(comm,orig,tag);

  assert(false);
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_particle_data( 
    mpi::communicator& comm, 
    int orig,
    int tag)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::vector<mpi::request> reqs;
  for (int ispc=0; ispc<Nspecies(); ispc++) {
    auto& container = get_container(ispc);

    // static array
    //container.incoming_particles.resize( container.first_message_size );

    reqs.emplace_back(
        comm.irecv(orig, get_tag(tag, ispc),
          container.incoming_particles.data(),
          container.first_message_size)
        );
  }

#ifdef GPU
  nvtxRangePop();
#endif

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_particle_extra_data( 
    mpi::communicator& comm, 
    int orig,
    int tag)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::vector<mpi::request> reqs;

  // this assumes that wait for the first message is already called and passed.

  // normal particles
  for (int ispc=0; ispc<Nspecies(); ispc++) {
    auto& container = get_container(ispc);

    //std::cout << "recv_prtcl1: got " << container.incoming_particles[0].id << " by mpi\n";

    int first_message_size = container.first_message_size;
    int np_tot = container.incoming_particles[0].id;
    int np_extra = np_tot > first_message_size ? np_tot - first_message_size : 0;

    if(np_extra > 0) {
      container.incoming_extra_particles.resize(np_extra);

      reqs.emplace_back(
          comm.irecv(orig, get_extra_tag(tag, ispc),
            container.incoming_extra_particles.data(),
            np_extra)
          );
    } else {
      container.incoming_extra_particles.clear();
    }

    //std::cout << this->communication.cid << " recv " << container.incoming_particles.size() << " + " << container.incoming_extra_particles.size() << " particles\n";
  }

#ifdef GPU
  nvtxRangePop();
#endif

  return reqs;
}

template<std::size_t D>
void Tile<D>::pack_all_particles()
{
  for(auto&& container : containers) 
    container.pack_all_particles();
}



template<std::size_t D>
void Tile<D>::pack_outgoing_particles()
{
  for(auto&& container : containers) 
    container.pack_outgoing_particles();

}


template<std::size_t D>
void Tile<D>::unpack_incoming_particles()
{
  for(auto&& container : containers) 
    container.unpack_incoming_particles();

}


template<std::size_t D>
void Tile<D>::delete_all_particles()
{
  for(auto&& container : containers) 
  {
    container.resize(0);
    container.Nprtcls = 0;
  }

}


template<std::size_t D>
void Tile<D>::shrink_to_fit_all_particles()
{
  for(auto&& container : containers) {

    // mpi extra message containers
    container.incoming_extra_particles.shrink_to_fit();
    container.outgoing_extra_particles.shrink_to_fit();

    // internal main particle containers
    container.shrink_to_fit();
  }
}



} // end of ns pic



template class pic::Tile<1>;
template class pic::Tile<2>;
template class pic::Tile<3>;
