#include <algorithm>
#include <map>
#include <utility>
#include <mpi.h>
#include <functional>
#include <type_traits>

#include "core/pic/particle.h"
#include "tools/wrap.h"
#include "external/iter/devcall.h"
#include "external/iter/iter.h"


#ifdef GPU
#include <cuda_runtime_api.h>
#include <nvtx3/nvToolsExt.h> 
#include "../tools/cub/cub.cuh"
#endif


namespace pic {


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

#ifdef DEBUG
  // NOTE: MPI messaging and storing of particles to our own manual vector 
  // requires that Particle is a POD and trivially copyable
  static_assert( std::is_pod_v<Particle>                == true );
  static_assert( std::is_trivially_copyable_v<Particle> == true );
  static_assert( std::is_trivial_v<Particle>            == true );
  static_assert( std::is_standard_layout_v<Particle>    == true );
#endif

  //incoming_particles.resize(first_message_size);
  incoming_extra_particles.reserve(first_message_size); // pre-allocating 

  //outgoing_particles.resize(first_message_size);
  outgoing_extra_particles.reserve(first_message_size); // pre-allocating

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
  infoArr.reserve(N);
    
  // reserve 1d N x D array for particle-specific emf
  Epart.reserve(N*3);
  Bpart.reserve(N*3);

  Nprtcls_cap = N; // mark the capacity of the array

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
  infoArr.resize(N);

  Epart.resize(N*3);
  Bpart.resize(N*3);

  //std::cout << " INFO: " << cid << " resizing container from " << Nprtcls << " to  " << N << std::endl;

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
  infoArr.shrink_to_fit();

#ifdef GPU
  nvtxRangePop();
#endif
}


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
    std::vector<float> prtcl_loc,
    std::vector<float> prtcl_vel,
    float prtcl_wgt)
{

#ifdef DEBUG
  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  assert(!std::isnan(prtcl_loc[0]));
  assert(!std::isnan(prtcl_loc[1]));
  assert(!std::isnan(prtcl_loc[2]));

  assert(!std::isnan(prtcl_vel[0]));
  assert(!std::isnan(prtcl_vel[1]));
  assert(!std::isnan(prtcl_vel[2]));

  assert(!std::isnan(prtcl_wgt));
#endif


  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  // get unique running key
  auto unique_key = keygen();
  indArr[0].push_back(std::get<0>(unique_key));
  indArr[1].push_back(std::get<1>(unique_key));

  infoArr.push_back(0); // extra

  Nprtcls++;
}


template<std::size_t D>
void ParticleContainer<D>::add_identified_particle (
    std::vector<float> prtcl_loc,
    std::vector<float> prtcl_vel,
    float prtcl_wgt,
    int _id, int _proc)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

#ifdef DEBUG
  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  assert(!std::isnan(prtcl_loc[0]));
  assert(!std::isnan(prtcl_loc[1]));
  assert(!std::isnan(prtcl_loc[2]));

  assert(!std::isnan(prtcl_vel[0]));
  assert(!std::isnan(prtcl_vel[1]));
  assert(!std::isnan(prtcl_vel[2]));

  assert(!std::isnan(prtcl_wgt));

  assert(!std::isnan(_id));
  assert(!std::isnan(_proc));
#endif

  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  indArr[0].push_back(_id);
  indArr[1].push_back(_proc);

  infoArr.push_back(0); // extra

  Nprtcls++;

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::insert_identified_particle (
    std::vector<float> prtcl_loc,
    std::vector<float> prtcl_vel,
    float prtcl_wgt,
    int _id, int _proc, int ind)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

#ifdef DEBUG
  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  assert(!std::isnan(prtcl_loc[0]));
  assert(!std::isnan(prtcl_loc[1]));
  assert(!std::isnan(prtcl_loc[2]));

  assert(!std::isnan(prtcl_vel[0]));
  assert(!std::isnan(prtcl_vel[1]));
  assert(!std::isnan(prtcl_vel[2]));

  assert(!std::isnan(prtcl_wgt));

  assert(!std::isnan(_id));
  assert(!std::isnan(_proc));

  assert(Nprtcls_cap > ind); // require that there are free slots for 
#endif

  locArr[0][ind] = prtcl_loc[0];
  locArr[1][ind] = prtcl_loc[1];
  locArr[2][ind] = prtcl_loc[2];

  velArr[0][ind] = prtcl_vel[0];
  velArr[1][ind] = prtcl_vel[1];
  velArr[2][ind] = prtcl_vel[2];

  wgtArr[ind] = prtcl_wgt;

  indArr[0][ind] = _id;
  indArr[1][ind] = _proc;

  infoArr[ind] = 0; // extra
                          
  //Nprtcls++; // NOTE: insertion needs to be added manually 

#ifdef GPU
  nvtxRangePop();
#endif
}


#pragma omp declare simd
inline int dir2info(int i, int j, int k){
  return 1 + (i+1) + 3*( (j+1) + 3*(k+1) );
}

#pragma omp declare simd
inline auto info2dir(int n){
  if(n==0) return std::make_tuple(0,0,0);

  n--; // reserve 1 element (=0) for particles that are transferred and do not have incoming dir
         
  int k = n / (3*3);
  int j = (n - k*3*3) / 3;
  int i = (n - k*3*3) % 3;

  return std::make_tuple(i-1,j-1,k-1);
}

#pragma omp declare simd
inline bool is_prtcl_inside(int n){
  return (n == 0) || (n==14); // default val or ijk=0
}



// check outgoing particles and update internal markers if particles are overflowing
template<size_t D>
void ParticleContainer<D>::check_outgoing_particles(
std::array<double,3>& mins,
std::array<double,3>& maxs)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

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

  assert(false); // not implemented with template parameter D

  // shortcut for particle locations
  float* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  particleIndexesA.resize(size());
  particleIndexesB.resize(size());

  UniIter::iterate([=] DEVCALLABLE (int ii, ParticleContainer<D> &self){
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


  UniIter::iterate([=] DEVCALLABLE (int ii, ParticleContainer<D> &self){
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

    
  //--------------------------------------------------
  //cpu version with no usage of external storage

  outgoing_count = 0;

  //UniIter::iterate([=] DEVCALLABLE (int n, ParticleContainer<D> &self){
#pragma omp simd reduction(+:outgoing_count)
  for(int n=0; n<size(); n++){
    int i=0,j=0,k=0; // relative indices

    if( (D>=1) &&  ( loc(0,n) - float( mins[0] ) <  0.0 )) i--; // left wrap
    if( (D>=1) &&  ( loc(0,n) - float( maxs[0] ) >= 0.0 )) i++; // right wrap
    if( (D>=2) &&  ( loc(1,n) - float( mins[1] ) <  0.0 )) j--; // bottom wrap
    if( (D>=2) &&  ( loc(1,n) - float( maxs[1] ) >= 0.0 )) j++; // top wrap
    if( (D>=3) &&  ( loc(2,n) - float( mins[2] ) <  0.0 )) k--; // back wrap
    if( (D>=3) &&  ( loc(2,n) - float( maxs[2] ) >= 0.0 )) k++; // front wrap
                                                    
    int info = dir2info(i,j,k);
    infoArr[n] = info;
    outgoing_count += !is_prtcl_inside(info);
  } //, size(), *this);

  //std::cout << "INFO " << cid << " outgoing count:" << outgoing_count << "\n";
#endif


#ifdef GPU
  nvtxRangePop();
#endif
}

//--------------------------------------------------

template<size_t D>
inline DEVCALLABLE float ParticleContainer<D>::get_prtcl_ene(size_t n)
{

  //const float mass = (type == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon

  const float mass = m; // read mass from class
  return std::sqrt( 
    mass*mass + 
    vel(0,n)*vel(0,n) + 
    vel(1,n)*vel(1,n) + 
    vel(2,n)*vel(2,n) 
    );
}

template<size_t D>
void ParticleContainer<D>::sort_in_rev_energy()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif
  
  eneArr.resize( size() ); // NOTE: do not create here but assume that it is initialized in constructor

  UniIter::iterate([=] DEVCALLABLE (size_t n, ParticleContainer<D>& self){
      self.eneArr[n] = get_prtcl_ene( n );
  }, size(), *this);
  UniIter::sync();

  //--------------------------------------------------
  //sort and apply
  auto indices = argsort_rev(eneArr);
  apply_permutation(indices);

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::delete_transferred_particles()
{
  if(size() == 0) return; // nothing to do; early return
                            
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  //--------------------------------------------------
  float* locn[3];
  for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );
  
  float* veln[3];
  for(int i=0; i<3; i++) veln[i] = &( vel(i,0) );
  
  int* idn[2];
  for(int i=0; i<2; i++) idn[i] = &( id(i,0) );

  //--------------------------------------------------
  // deletion using the remove_if idiom from std::remove
  int first = 0;
  int last = size();
  int iter = first;

  while( first != last ) {
    
    // check if outside tile 
    int n = infoArr[first];

    // replace good value with a bad value
    if( is_prtcl_inside(n) ) {
      if( iter != first ){

        // should be move operation; same?
        locn[0][iter] = locn[0][first];
        locn[1][iter] = locn[1][first];
        locn[2][iter] = locn[2][first];

        veln[0][iter] = veln[0][first];
        veln[1][iter] = veln[1][first];
        veln[2][iter] = veln[2][first];

        idn[0][ iter] =  idn[0][first];
        idn[1][ iter] =  idn[1][first];

        wgtArr[ iter] =  wgtArr[first];
        infoArr[iter] = infoArr[first];
      }
      iter++;
    }
    first++;
  }
  int new_last = iter;

  // now remove anything between result and last
  resize(new_last);
  Nprtcls -= (last - new_last); // substract number of particles removed

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void ParticleContainer<D>::transfer_and_wrap_particles( 
    ParticleContainer&    neigh,
    std::array<int,3>     dirs, 
    std::array<double,3>& global_mins, 
    std::array<double,3>& global_maxs
    )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // global grid limits
  const float minx = global_mins[0];
  const float miny = global_mins[1];
  const float minz = global_mins[2];

  const float maxx = global_maxs[0];
  const float maxy = global_maxs[1];
  const float maxz = global_maxs[2];

  //--------------------------------------------------
  // v2: loop w/o external to_other_tiles storage

  for(int n=0; n<neigh.size(); n++){
    auto [i,j,k] = info2dir( neigh.infoArr[n] );

    // these need to be skipped 
    if( (D==1) && (i== 0) )                     continue;
    if( (D==2) && (i== 0) && (j==0) )           continue;
    if( (D==3) && (i== 0) && (j==0) && (k==0) ) continue;

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile
    bool add=false;
    if( (D==1) && (i== -dirs[0]) )                                   add = true;
    if( (D==2) && (i== -dirs[0]) && (j==-dirs[1]) )                  add = true;
    if( (D==3) && (i== -dirs[0]) && (j==-dirs[1]) && (k==-dirs[2]) ) add = true;
    
    if(add) {
      float locx = D >= 1 ? wrap( neigh.loc(0, n), minx, maxx ) : neigh.loc(0, n);
      float locy = D >= 2 ? wrap( neigh.loc(1, n), miny, maxy ) : neigh.loc(1, n);
      float locz = D >= 3 ? wrap( neigh.loc(2, n), minz, maxz ) : neigh.loc(2, n);

      add_identified_particle(
          {locx, locy, locz}, 
          {neigh.vel(0, n), neigh.vel(1, n), neigh.vel(2, n)},
          neigh.wgt(n),
          neigh.id(0,n), neigh.id(1,n));
    }
  }

#ifdef GPU
  nvtxRangePop();
#endif
}

//--------------------------------------------------

template<std::size_t D>
void ParticleContainer<D>::pack_all_particles()
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // clear extra array (for appending)
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  int np = size() + 1;

  if(np-first_message_size > 0) {
    // reserve is needed here; if size is less than capacity, we do nothing
    outgoing_extra_particles.reserve( np-first_message_size );
  }

  // first particle is always the message info
  outgoing_particles[0] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np, 0}; // store prtcl number in id slot

  // next, pack all other particles
  int i=1;
  for(size_t ind=0; ind < size(); ind++) {
    if(i < first_message_size) {
      //outgoing_particles.push_back({ 
      outgoing_particles[i] = {
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) };
    } else {
      //outgoing_extra_particles.push_back({ 
      outgoing_extra_particles.push_back({ 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) });
    }
    i++;
  }

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

  // initialize msg arrays 
  outgoing_extra_particles.clear(); 

  //add info prtcl that describes how many total particles are coming (and therefore are in the variable extra msg)
  outgoing_particles[0] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0}; // store prtcl number in id slot

                                                                              
  // -------------------------------------------------- 
  // on-the-fly calculation of escape condition

  // next, pack all other particles
  int count=1;
  for(int n=0; n<size(); n++){

    // check if moving out
    bool to_be_packed = !is_prtcl_inside(infoArr[n]);

    // pack if true; split between fixed primary and adaptive extra message containers
    if(to_be_packed) {
      if(count < first_message_size) {
        outgoing_particles[count] = {
          loc(0, n), loc(1, n), loc(2, n), 
          vel(0, n), vel(1, n), vel(2, n), 
          wgt(n), 
          id(0, n), id(1, n) };
      } else {
        outgoing_extra_particles.push_back({ 
          loc(0, n), loc(1, n), loc(2, n), 
          vel(0, n), vel(1, n), vel(2, n), 
          wgt(n), 
          id(0, n), id(1, n) });
      }
      count++;
    }
  }
  outgoing_particles[0].id = count; // update info
    
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

  // get real number of incoming particles
  int np_tot = incoming_particles[0].id; // number stored in id slot

  // first and extra msg sizes
  int np_first  = np_tot > first_message_size ? first_message_size : np_tot;
  int np_extra = np_tot > first_message_size ? np_tot - first_message_size : 0;
  int N = size(); 

  reserve( N + np_tot - 1);  //reserve for addition

  // skipping 1st info particle
  for(int i=1; i<np_first; i++){
    float locx = incoming_particles[i].x;
    float locy = incoming_particles[i].y;
    float locz = incoming_particles[i].z;

    float velx = incoming_particles[i].ux;
    float vely = incoming_particles[i].uy;
    float velz = incoming_particles[i].uz;
    float wgts = incoming_particles[i].w;

    int ids  = incoming_particles[i].id;
    int proc = incoming_particles[i].proc;

    add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc);
  }

  for(int i=0; i<np_extra; i++){
    float locx = incoming_extra_particles[i].x;
    float locy = incoming_extra_particles[i].y;
    float locz = incoming_extra_particles[i].z;

    float velx = incoming_extra_particles[i].ux;
    float vely = incoming_extra_particles[i].uy;
    float velz = incoming_extra_particles[i].uz;
    float wgts = incoming_extra_particles[i].w;

    int ids  = incoming_extra_particles[i].id;
    int proc = incoming_extra_particles[i].proc;
      
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



template<size_t D>
void ParticleContainer<D>::apply_permutation( ManVec<size_t>& indices )
{

  // check that sizes match
  assert( indices.size() == size() );

  // https://stackoverflow.com/questions/67751784/how-to-do-in-place-sorting-a-list-according-to-a-given-index-in-c
  // and
  // https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095
  for (size_t i=0; i<indices.size(); i++) {
    auto current = i;

    while (i != indices[current]) {
      auto next = indices[current];

      // operation we perform is:
      //swap(v[current], v[next]);

      std::swap( locArr[0][current], locArr[0][next] ); 
      std::swap( locArr[1][current], locArr[1][next] ); 
      std::swap( locArr[2][current], locArr[2][next] ); 

      std::swap( velArr[0][current], velArr[0][next] ); 
      std::swap( velArr[1][current], velArr[1][next] ); 
      std::swap( velArr[2][current], velArr[2][next] ); 

      std::swap( indArr[0][current], indArr[0][next] ); 
      std::swap( indArr[1][current], indArr[1][next] ); 

      std::swap( wgtArr[current], wgtArr[next] ); 

      std::swap( eneArr[current], eneArr[next] ); 

      std::swap( infoArr[current], infoArr[next] );  // TODO might be possible to skip this swap

      // NOTE: these can be omitted if interpolator is called *after* sort
      //std::swap( ex[current], ex[next] ); 
      //std::swap( ey[current], ey[next] ); 
      //std::swap( ez[current], ez[next] ); 

      //std::swap( bx[current], bx[next] ); 
      //std::swap( by[current], by[next] ); 
      //std::swap( bz[current], bz[next] ); 

      indices[current] = current;
      current = next;
    }

  indices[current] = current;
  }

  return;
}


template<size_t D>
void ParticleContainer<D>::update_cumulative_arrays()
{
  const size_t N = size(); // number of prtcls
  wgtCumArr.resize(N); 

  // sum over wgt
  float wsum = 0.0;

  #pragma omp simd reduction(+:wsum)
  for(size_t i=0; i<N; i++) {
    wsum += wgtArr[i];
  }

  // normalized cumulative sum
  wgtCumArr[0] = wgtArr[0];
  for(size_t i=1; i<N; i++) wgtCumArr[i] = wgtArr[i] + wgtCumArr[i-1];

  return;
}





} // end ns pic


template class pic::ParticleContainer<1>;
template class pic::ParticleContainer<2>;
template class pic::ParticleContainer<3>;
