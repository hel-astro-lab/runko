
#include "particle.h"
#include "../tools/wrap.h"
#include "../tools/iter/devcall.h"
#include "../tools/iter/iter.h"

#include <algorithm>
#include <map>
#include <utility>
#include <mpi.h>
#include <functional>

#ifdef GPU
#include <cuda_runtime_api.h>
#include <nvtx3/nvToolsExt.h> 
#include "../tools/cub/cub.cuh"
#endif


namespace pic {

inline Particle::Particle(
    float_p x,  float_p y,  float_p z,
    float_p ux, float_p uy, float_p uz, 
    float_p wgt,
    int __ind, int __proc
    ) : 
  _id(__ind),
  _proc(__proc)
{
  data[0] = x;
  data[1] = y;
  data[2] = z;
  data[3] = ux;
  data[4] = uy;
  data[5] = uz;
  data[6] = wgt;
}


inline Particle::Particle( size_t number_of_particles) 
{
  data[0] = static_cast<float_p>(number_of_particles);
}

/// special method for info particle that re-uses x mem location
size_t Particle::number_of_particles() {
  return static_cast<size_t>( data[0] );
}



//--------------------------------------------------
//--------------------------------------------------
// ParticleContainer methods

template<std::size_t D>
ParticleContainer<D>::ParticleContainer()
{ 

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // Get the number of processes
  //MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);

  // Get the rank of the process
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  incoming_particles.resize(optimal_message_size);
  incoming_extra_particles.resize(optimal_message_size);

  outgoing_particles.resize(optimal_message_size);
  outgoing_extra_particles.resize(optimal_message_size);

#ifdef GPU
  //DEV_REGISTER
  temp_storage_bytes = 10000;
  //std::cout << temp_storage_bytes << std::endl;
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
#endif


#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::reserve(size_t N) {

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // always reserve at least 1 element to ensure proper array initialization
  if (N <= 0) N = 1;

  for(size_t i=0; i<3; i++) locArr[i].reserve(N);
  for(size_t i=0; i<3; i++) velArr[i].reserve(N);
  for(size_t i=0; i<2; i++) indArr[i].reserve(N);
  wgtArr.reserve(N);
    
  // reserve 1d N x D array for particle-specific fields
  Epart.reserve(N*3);
  Bpart.reserve(N*3);

#ifdef GPU
  nvtxRangePop();
#endif
}

template<std::size_t D>
void ParticleContainer<D>::resize(size_t N)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  for(size_t i=0; i<3; i++) locArr[i].resize(N);
  for(size_t i=0; i<3; i++) velArr[i].resize(N);
  for(size_t i=0; i<2; i++) indArr[i].resize(N);
  wgtArr.resize(N);

  Epart.resize(N*3);
  Bpart.resize(N*3);

  Nprtcls = N;

#ifdef GPU
  nvtxRangePop();
#endif
}

template<std::size_t D>
void ParticleContainer<D>::shrink_to_fit()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  for(size_t i=0; i<3; i++) locArr[i].shrink_to_fit();
  for(size_t i=0; i<3; i++) velArr[i].shrink_to_fit();
  for(size_t i=0; i<2; i++) indArr[i].shrink_to_fit();
  wgtArr.shrink_to_fit();

#ifdef GPU
  nvtxRangePop();
#endif
}


//template<std::size_t D>
//size_t ParticleContainer<D>::size() 
//{ 
//  assert(locArr[0].size() == Nprtcls);
//  assert(locArr[1].size() == Nprtcls);
//  assert(locArr[2].size() == Nprtcls);
//
//  assert(velArr[0].size() == Nprtcls);
//  assert(velArr[1].size() == Nprtcls);
//  assert(velArr[2].size() == Nprtcls);
//
//  assert(indArr[0].size() == Nprtcls);
//  assert(indArr[1].size() == Nprtcls);
//
//  assert(wgtArr.size() == Nprtcls);
//
//  return Nprtcls; 
//}


template<std::size_t D>
std::pair<int,int> pic::ParticleContainer<D>::keygen() 
{
    // get running key and increment internal counter
  int unique_key = _key;
  _key++; // TODO: add atomic around this to assure non-overlapping keys

  return std::make_pair(unique_key, _rank);
}


template<std::size_t D>
void ParticleContainer<D>::add_particle (
    std::vector<float_p> prtcl_loc,
    std::vector<float_p> prtcl_vel,
    float_p prtcl_wgt)
{

  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  // get unique running key
  auto unique_key = keygen();
  indArr[0].push_back(std::get<0>(unique_key));
  indArr[1].push_back(std::get<1>(unique_key));

  Nprtcls++;
}

template<std::size_t D>
void ParticleContainer<D>::add_identified_particle (
    std::vector<float_p> prtcl_loc,
    std::vector<float_p> prtcl_vel,
    float_p prtcl_wgt,
    int _id, int _proc)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  indArr[0].push_back(_id);
  indArr[1].push_back(_proc);

  Nprtcls++;

#ifdef GPU
  nvtxRangePop();
#endif
}

template<>
void ParticleContainer<1>::check_outgoing_particles(
    std::array<double,3>& /*mins*/,
    std::array<double,3>& /*maxs*/)
{
  // TODO implement
  assert(false);
}

template<>
void ParticleContainer<2>::check_outgoing_particles(
    std::array<double,3>& mins,
    std::array<double,3>& maxs)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  to_other_tiles.clear();

  // unpack limits
  double 
    xmin = mins[0],
    ymin = mins[1],

    xmax = maxs[0],
    ymax = maxs[1];

  int lenx = static_cast<int>( xmax - xmin );
  int leny = static_cast<int>( ymax - ymin );

  int i0, j0;

  // shortcut for particle locations
  float_p* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  int i,j,k; // relative indices
  for(size_t n=0; n<size(); n++) {
    i = 0;
    j = 0;

    i0 = static_cast<int>( floor(locn[0][n] - mins[0]) );
    j0 = static_cast<int>( floor(locn[1][n] - mins[1]) );

    if(i0 <  0)    i--; // left wrap
    if(i0 >= lenx) i++; // right wrap

    if(j0 <  0)    j--; // bottom wrap
    if(j0 >= leny) j++; // top wrap

    // collapse z dimension
    if ((i == 0) && (j == 0)) continue; 

    if ( (i != 0) || (j != 0) || (k != 0) ) to_other_tiles.push_back( {i,j,k,n} );
  }

#ifdef GPU
  nvtxRangePop();
#endif
}



template<std::size_t D>
void ParticleContainer<D>::delete_particles(std::vector<int> to_be_deleted) 
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::sort(to_be_deleted.begin(), to_be_deleted.end(), std::greater<int>() );

  float_p* locn[3];
  for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  float_p* veln[3];
  for(int i=0; i<3; i++) veln[i] = &( vel(i,0) );

  int* idn[2];
  for(int i=0; i<2; i++) idn[i] = &( id(i,0) );


  // overwrite particles with the last one on the array and 
  // then resize the array
  int last = size();
  for(int indx : to_be_deleted) {
    last--;
    if(indx == last) continue;
    //std::cout << "deleting " << indx 
    //          << " by putting it to " << last << '\n';
    for(int i=0; i<3; i++) locn[i][indx] = locn[i][last];
    for(int i=0; i<3; i++) veln[i][indx] = veln[i][last];
    wgtArr[indx] = wgtArr[last];
    for(int i=0; i<2; i++) idn[i][indx] = idn[i][last];
  }

  // resize if needed and take care of the size
  last = last < 0 ? 0 : last;
  if ((last != (int)size()) && (size() > 0)) resize(last);

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::delete_transferred_particles()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::vector<int> to_be_deleted;

  // get transferred 

  //std::sort(to_other_tiles.begin(), to_other_tiles.end(), [](auto a, auto b){return a.n > b.n;} );

  float_p* locn[3];
  for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  float_p* veln[3];
  for(int i=0; i<3; i++) veln[i] = &( vel(i,0) );

  int* idn[2];
  for(int i=0; i<2; i++) idn[i] = &( id(i,0) );


  // overwrite particles with the last one on the array and 
  // then resize the array
  int last = size()-to_other_tiles.size();

  UniIter::iterate([=] DEVCALLABLE (int i, ManVec<to_other_tiles_struct> &to_other_tiles){
    int other = last+i;////size() - 1 - i;
    int indx = to_other_tiles[i].n;

  //int last = size();
  //for(auto elem : to_other_tiles) {
  //  int indx = elem.n;
  //  last--;
    if(indx >= last) return;
    //std::cout << "deleting " << indx << " by putting it to " << last << '\n';
    for(int i=0; i<3; i++) locn[i][indx] = locn[i][other];
    for(int i=0; i<3; i++) veln[i][indx] = veln[i][other];
    wgtArr[indx] = wgtArr[other];
    for(int i=0; i<2; i++) idn[i][indx] = idn[i][other];
  }, to_other_tiles.size(), to_other_tiles);

  UniIter::sync();

  // resize if needed and take care of the size
  last = last < 0 ? 0 : last;
  if ((last != (int)size()) && (size() > 0)) resize(last);

  /*
  for (size_t ii = 0; ii < to_other_tiles.size(); ii++)
  {
    const auto &elem = to_other_tiles[ii];
    to_be_deleted.push_back( elem.n );
  }
  delete_particles(to_be_deleted);
  */


#ifdef GPU
  nvtxRangePop();
#endif
}

template<>
void ParticleContainer<1>::transfer_and_wrap_particles( 
    ParticleContainer&    /*neigh*/,
    std::array<int,3>     /*dirs*/, 
    std::array<double,3>& /*global_mins*/, 
    std::array<double,3>& /*global_maxs*/
    )
{
  // TODO not implemented
  assert(false);
}


template<>
void ParticleContainer<2>::transfer_and_wrap_particles( 
    ParticleContainer& neigh,
    std::array<int,3>    dirs, 
    std::array<double,3>& global_mins, 
    std::array<double,3>& global_maxs
    )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // particle overflow from tiles is done in shortest precision
  // to avoid rounding off errors and particles left in a limbo
  // between tiles.
  float_p locx, locy, locz, velx, vely, velz, wgt;
  int id, proc;

  int i;
  //for (auto&& elem : neigh.to_other_tiles) {
  for (size_t ii = 0; ii < neigh.to_other_tiles.size(); ii++)
  {
    const auto &elem = neigh.to_other_tiles[ii];

    if(elem.i == 0 && elem.j == 0) continue; 

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile

    if (elem.i == -dirs[0] &&
        elem.j == -dirs[1] ) {

      i = elem.n;

      // NOTE: wrap bounds to [min, max)
      locx = wrap( neigh.loc(0, i), static_cast<float_p>(global_mins[0]), static_cast<float_p>(global_maxs[0]) );
      locy = wrap( neigh.loc(1, i), static_cast<float_p>(global_mins[1]), static_cast<float_p>(global_maxs[1]) );
      locz = wrap( neigh.loc(2, i), static_cast<float_p>(global_mins[2]), static_cast<float_p>(global_maxs[2]) );

      velx = neigh.vel(0, i);
      vely = neigh.vel(1, i);
      velz = neigh.vel(2, i);

      wgt  = neigh.wgt(i);

      id   = neigh.id(0,i);
      proc = neigh.id(1,i);

      add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgt, id, proc);
    }
  }

#ifdef GPU
  nvtxRangePop();
#endif
}




template<std::size_t D>
void ParticleContainer<D>::pack_all_particles()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif


  outgoing_particles.clear();
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  size_t np = size() + 1;

  // FIXME
  //outgoing_particles.reserve(optimal_message_size);
  if(np-optimal_message_size > 0) {
    outgoing_extra_particles.reserve( np-optimal_message_size );
  }

  // first particle is always the message info
  outgoing_particles.push_back({np});

  // next, pack all other particles
  int i=1;
  for(size_t ind=0; ind < size(); ind++) {
    if(i < optimal_message_size) {
      outgoing_particles.push_back({ 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) });
    } else {
      outgoing_extra_particles.push_back({ 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) });
    }
    i++;
  }
  //outgoing_extra_particles.shrink_to_fit();

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::pack_outgoing_particles()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  outgoing_particles.clear();
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  size_t np = to_other_tiles.size() + 1;

  //std::cout << "reserving1: " << optimal_message_size << "\n";
  //std::cout << "reserving2: " << np <<" minus " << np-optimal_message_size << "\n";


  // FIXME altered here from v0 w/ reserve to new ver w/ resize
  //outgoing_particles.reserve(optimal_message_size);
  if (np > optimal_message_size) {
    // FIXME altered here from v0 w/ reserve to new ver w/ resize
    outgoing_extra_particles.resize( np-optimal_message_size );
  }

  // first particle is always the message info
  outgoing_particles.push_back({np});

  // next, pack all other particles
  int i=1, ind;
  
//  for (auto&& elem : to_other_tiles) {
  for (size_t ii = 0; ii < to_other_tiles.size(); ii++)
  {
    const auto &elem = to_other_tiles[ii];
    ind = elem.n;

    if(i < optimal_message_size) {
      outgoing_particles.push_back({ 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) });
    } else {
      outgoing_extra_particles.push_back({ 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) });
    }

    i++;
  }

  //outgoing_extra_particles.shrink_to_fit();

  // TODO: set next message size dynamically according to history
  //optimal_message_size = np;
  //
#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::unpack_incoming_particles()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  float_p locx, locy, locz, velx, vely, velz, wgts;
  int ids, proc;

  // get real number of incoming particles
  int number_of_incoming_particles = incoming_particles[0].number_of_particles();

  int number_of_primary_particles = 
    number_of_incoming_particles > optimal_message_size 
    ? optimal_message_size : number_of_incoming_particles;

  int number_of_secondary_particles = incoming_extra_particles.size();

  // skipping 1st info particle
  for(int i=1; i<number_of_primary_particles; i++){
    locx = incoming_particles[i].x();
    locy = incoming_particles[i].y();
    locz = incoming_particles[i].z();

    velx = incoming_particles[i].ux();
    vely = incoming_particles[i].uy();
    velz = incoming_particles[i].uz();
    wgts = incoming_particles[i].wgt();

    ids  = incoming_particles[i].id();
    proc = incoming_particles[i].proc();

    add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc);
  }

  for(int i=0; i<number_of_secondary_particles; i++){
    locx = incoming_extra_particles[i].x();
    locy = incoming_extra_particles[i].y();
    locz = incoming_extra_particles[i].z();

    velx = incoming_extra_particles[i].ux();
    vely = incoming_extra_particles[i].uy();
    velz = incoming_extra_particles[i].uz();
    wgts = incoming_extra_particles[i].wgt();

    ids  = incoming_extra_particles[i].id();
    proc = incoming_extra_particles[i].proc();

    add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc);
  }

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::set_keygen_state(int __key, int __rank)
{
  _key  = __key;
  _rank = __rank;
}


template<>
void ParticleContainer<3>::check_outgoing_particles(
    std::array<double,3>& mins,
    std::array<double,3>& maxs)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  to_other_tiles.clear();
  outgoing_count = 0;

  // shortcut for particle locations
  float_p* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );


/*
  int maxCap = to_other_tiles.capacity();
  to_other_tiles.resize(to_other_tiles.capacity());
  auto listPtr = to_other_tiles.data();
  
  // redesigned for parallelism
  // first try and add particle to to_other_tiles as long as there is space
  // if we run out of space keep counting
  // if count is larger than the space, reallocate for count, clear the whole thign and rerun.
  //
  UniIter::iterate(
    [=] DEVFUNIFGPU (int n, int &count)
    {
      int i = 0;
      int j = 0;
      int k = 0;

      int i0, j0, k0;
      i0 = static_cast<int>( floor(locn[0][n] - mins[0]) );
      j0 = static_cast<int>( floor(locn[1][n] - mins[1]) );
      k0 = static_cast<int>( floor(locn[2][n] - mins[2]) );

      if(i0 <  0)    i--; // left wrap
      if(i0 >= lenx) i++; // right wrap

      if(j0 <  0)    j--; // bottom wrap
      if(j0 >= leny) j++; // top wrap

      if(k0 <  0)    k--; // back
      if(k0 >= lenz) k++; // front

      if ( (i != 0) || (j != 0) || (k != 0) ) 
      {
        //to_other_tiles.push_back( {i,j,k,n} );
        #ifdef GPU
          int pos = atomicAdd(&count, 1);
        #else
        int pos;
          #pragma omp atomic capture 
          {
            pos = count;
            count++;
          }
        #endif
        
        if(pos < maxCap)
          listPtr[pos] = {i,j,k,n};
      }
    },size(), outgoing_count);
    UniIter::sync();
    // check outgoing_count and react to it...
    if(outgoing_count > maxCap)
    {
      UniIter::sync();
      //std::cout << outgoing_count << " outgoing count over cap " << maxCap << std::endl;
      to_other_tiles.clear();
      to_other_tiles.reserve((outgoing_count));

      // clear, realloc and recall
      check_outgoing_particles(mins, maxs);
    }
    else{
      to_other_tiles.resize(outgoing_count);
    }
    */


#ifdef GPU
  particleIndexesA.resize(size());
  particleIndexesB.resize(size());

  UniIter::iterate([=] DEVCALLABLE (int ii, ParticleContainer<3> &self){
    self.particleIndexesA[ii] = ii;
  }, size(), *this);

  cub::DeviceSelect::If(
          d_temp_storage, temp_storage_bytes, 
          particleIndexesA.data(), 
          particleIndexesB.data(), 
  &pCount, size(), 
   [=]__device__ (int n){
    int i,j,k; // relative indices
    i = 0;
    j = 0;
    k = 0;

    if( locn[0][n]-mins[0] <  0.0 ) i--; // left wrap
    if( locn[0][n]-maxs[0] >= 0.0 ) i++; // right wrap

    if( locn[1][n]-mins[1] <  0.0 ) j--; // bottom wrap
    if( locn[1][n]-maxs[1] >= 0.0 ) j++; // top wrap

    if( locn[2][n]-mins[2] <  0.0 ) k--; // back wrap
    if( locn[2][n]-maxs[2] >= 0.0 ) k++; // front wrap

    return ( (i != 0) || (j != 0) || (k != 0) ) ;
  });

  cudaDeviceSynchronize();
  
  to_other_tiles.resize(pCount);


  UniIter::iterate([=] DEVCALLABLE (int ii, ParticleContainer<3> &self){
    int n = self.particleIndexesB[ii];
    int i=0,j=0,k=0; // relative indices

    if( locn[0][n]-mins[0] <  0.0 ) i--; // left wrap
    if( locn[0][n]-maxs[0] >= 0.0 ) i++; // right wrap

    if( locn[1][n]-mins[1] <  0.0 ) j--; // bottom wrap
    if( locn[1][n]-maxs[1] >= 0.0 ) j++; // top wrap

    if( locn[2][n]-mins[2] <  0.0 ) k--; // back wrap
    if( locn[2][n]-maxs[2] >= 0.0 ) k++; // front wrap

    self.to_other_tiles[ii] =  {i,j,k,n};
  }, pCount, *this);


#else
  for(size_t n=0; n<size(); n++) {
    int i=0,j=0,k=0; // relative indices

    if( loc(0,n) - float_p( mins[0] ) <  0.0 ) i--; // left wrap
    if( loc(0,n) - float_p( maxs[0] ) >= 0.0 ) i++; // right wrap
    if( loc(1,n) - float_p( mins[1] ) <  0.0 ) j--; // bottom wrap
    if( loc(1,n) - float_p( maxs[1] ) >= 0.0 ) j++; // top wrap
    if( loc(2,n) - float_p( mins[2] ) <  0.0 ) k--; // back wrap
    if( loc(2,n) - float_p( maxs[2] ) >= 0.0 ) k++; // front wrap

    if ( (i != 0) || (j != 0) || (k != 0) ) {
        to_other_tiles.push_back( {i,j,k,n} );
        outgoing_count++;
    }
  }

  //std::cout << "outgoing count:" << outgoing_count << "\n";
#endif


#ifdef GPU
  nvtxRangePop();
#endif
}


template<>
void ParticleContainer<3>::transfer_and_wrap_particles( 
    ParticleContainer& neigh,
    std::array<int,3>    dirs, 
    std::array<double,3>& global_mins, 
    std::array<double,3>& global_maxs
    )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // particle overflow from tiles is done in shortest precision
  // to avoid rounding off errors and particles left in a limbo
  // between tiles.
  float_p locx, locy, locz;

  int ind;
//  for (size_t ii = 0; ii < neigh.to_other_tiles.size(); ii++)
//  {
//    const auto &elem = neigh.to_other_tiles[ii];
  for (auto&& elem : neigh.to_other_tiles) {
      
    if(elem.i == 0 && 
       elem.j == 0 &&
       elem.k == 0) continue; 

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile

    if (elem.i == -dirs[0] &&
        elem.j == -dirs[1] &&
        elem.k == -dirs[2] ) {

      ind = elem.n;

      locx = wrap( neigh.loc(0, ind), static_cast<float_p>(global_mins[0]), static_cast<float_p>(global_maxs[0]) );
      locy = wrap( neigh.loc(1, ind), static_cast<float_p>(global_mins[1]), static_cast<float_p>(global_maxs[1]) );
      locz = wrap( neigh.loc(2, ind), static_cast<float_p>(global_mins[2]), static_cast<float_p>(global_maxs[2]) );

      locArr[0].push_back(locx);
      locArr[1].push_back(locy);
      locArr[2].push_back(locz);

      //std::cout << locx << " " << locy << " " << locz << " " <<  velx << " " << vely << " " << velz << " " <<  wgt << " " <<  id << " " <<  proc << std::endl;
      //add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgt, id, proc);

      velArr[0].push_back( neigh.vel(0, ind) );
      velArr[1].push_back( neigh.vel(1, ind) );
      velArr[2].push_back( neigh.vel(2, ind) );

      wgtArr.push_back( neigh.wgt(ind));

      indArr[0].push_back(neigh.id(0,ind));
      indArr[1].push_back(neigh.id(1,ind));

      Nprtcls++;
    }
  }


#ifdef GPU
  nvtxRangePop();
#endif
}


} // end ns pic


template class pic::ParticleContainer<1>;
template class pic::ParticleContainer<2>;
template class pic::ParticleContainer<3>;
