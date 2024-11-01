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

  outgoing_count = 0;

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
  //v0 with to_other_tiles as a storage

  //to_other_tiles.reserve(0.3*size()); // estimate for the num of outflowing prtcls
  //for(size_t n=0; n<size(); n++) {
  //  int i=0,j=0,k=0; // relative indices
  //  
  //  if( loc(0,n) - float( mins[0] ) <  0.0 ) i--; // left wrap
  //  if( loc(0,n) - float( maxs[0] ) >= 0.0 ) i++; // right wrap
  //  if( loc(1,n) - float( mins[1] ) <  0.0 ) j--; // bottom wrap
  //  if( loc(1,n) - float( maxs[1] ) >= 0.0 ) j++; // top wrap
  //  if( loc(2,n) - float( mins[2] ) <  0.0 ) k--; // back wrap
  //  if( loc(2,n) - float( maxs[2] ) >= 0.0 ) k++; // front wrap

  //  if ( (i != 0) || (j != 0) || (k != 0) ) {
	//    to_other_tiles.push_back( {i,j,k,n} );
	//    outgoing_count++;
  //  }
  //}
    
  //--------------------------------------------------
  //v1 with no storage

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
    //outgoing_count += (i!=0) || (j!=0) || (k!=0 ) ? 1 : 0;
    outgoing_count += !is_prtcl_inside(info);
  } //, size(), *this);


  // second pass to reconstruct to_other_tiles
  // TODO remove
  //to_other_tiles.clear();
  //to_other_tiles.reserve(outgoing_count);
  //for(size_t n=0; n<size(); n++) {

  //  //int ii=0,jj=0,kk=0; // relative indices
  //  //if( loc(0,n) - float( mins[0] ) <  0.0 ) ii--; // left wrap
  //  //if( loc(0,n) - float( maxs[0] ) >= 0.0 ) ii++; // right wrap
  //  //if( loc(1,n) - float( mins[1] ) <  0.0 ) jj--; // bottom wrap
  //  //if( loc(1,n) - float( maxs[1] ) >= 0.0 ) jj++; // top wrap
  //  //if( loc(2,n) - float( mins[2] ) <  0.0 ) kk--; // back wrap
  //  //if( loc(2,n) - float( maxs[2] ) >= 0.0 ) kk++; // front wrap

  //  auto [i,j,k] = info2dir( infoArr[n] );

  //  //std::cout << "comp " << 
  //  //  " ijk v1: " <<i<<","<<j<<","<<k << " n:"<<n<<
  //  //  " ijk v0: " <<ii<<","<<jj<<","<<kk
  //  //  <<"\n";
  //  
  //  if ( (i != 0) || (j != 0) || (k != 0) ) {
	//    to_other_tiles.push_back( {i,j,k,n} );
  //  }
  //}

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

  //--------------------------------------------------
  // energy array for sorting
  //const float mass = (type == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
  
  // construct energy array
  //std::vector<float> eneArr( size() );  
  
  eneArr.resize( size() ); // NOTE: do not create here but assume that it is initialized in constructor
                             
  // regular non-simd version
  //for(size_t n=0; n<size(); n++) {
  //  eneArr[n] = get_prtcl_ene( n );
  //}

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
void ParticleContainer<D>::delete_particles(std::vector<int> to_be_deleted) 
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  std::sort(to_be_deleted.begin(), to_be_deleted.end(), std::greater<int>() );
  
  float* locn[3];
  for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );
  
  float* veln[3];
  for(int i=0; i<3; i++) veln[i] = &( vel(i,0) );
  
  int* idn[2];
  for(int i=0; i<2; i++) idn[i] = &( id(i,0) );


  // overwrite particles with the last one on the array and 
  // then resize the array
  int last = size()-to_be_deleted.size();

  // NOTE vectorizing the below loop leads to race conditions; hence it is turned off for now

  //UniIter::iterate([=] DEVCALLABLE (
  //      int ii, 
  //      std::vector<int>& to_be_deleted){
  for(int ii=0; ii<to_be_deleted.size(); ii++){

    int other = last+ii; //size() - 1 - i;
    int indx = to_be_deleted[ii];

    //if(indx >= last) return;
    //std::cout << " sw " << indx << " to " << other << " while last " << last << std::endl;
      
    //std::cout << "deleting " << indx << " by putting it to " << last << '\n';
    for(int i=0; i<3; i++) locn[i][indx] = locn[i][other];
    for(int i=0; i<3; i++) veln[i][indx] = veln[i][other];
    for(int i=0; i<2; i++) idn[ i][indx] = idn[ i][other];
    wgtArr[ indx] = wgtArr[ other];
    infoArr[indx] = infoArr[other];

  //}, to_be_deleted.size(), to_be_deleted);
  }

  UniIter::sync();
  
  // resize if needed and take care of the size
  last = last < 0 ? 0 : last;
  if ((last != (int)size()) && (size() > 0)) resize(last);

  //std::cout << " INFO: " << cid << " v1: removing prtcls :" << Nprtcls << " - " << to_be_deleted.size() << std::endl;
  Nprtcls -= to_be_deleted.size();

#ifdef GPU
  nvtxRangePop();
#endif
}



template<std::size_t D>
void ParticleContainer<D>::delete_transferred_particles()
{
  if(size() == 0) return; // nothing to do
                            
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

  //--------------------------------------------------
  // remove_if idiom
  int first = 0;
  int last = size();
  int iter = first;

  while( first != last ) {
    
    // check if outside tile 
    int n = infoArr[first];

    //auto [i,j,k] = info2dir(n);
    //bool inside_tile = (i == 0) && (j == 0) && (k == 0);
    //bool to_be_removed = !inside_tile;

    // replace good value with a bad value
    //if( !to_be_removed ) {
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

  //--------------------------------------------------
  //for(int n=0; n<size(); n++){
  //  int nn = infoArr[n];
  //  auto [i,j,k] = info2dir(nn);
  //  std::cout << "list after del: n: " << n << " ind: " << nn << " ijk " << i << "," << j << "," << k << " wgt: " << wgtArr[n] << "\n";
  //}

  // now remove anything between result and last
  resize(new_last);
  Nprtcls -= (last - new_last); // substract number of particles removed

  //std::cout << " INFO: " << cid << " v3: remove_if:" << Nprtcls << " iter " << iter << " last:" << last << " deleted:" << last-new_last << std::endl;



  //--------------------------------------------------
  // v0
    
  // do nothing if empty
  //if(to_other_tiles.size() == 0) return;
  //  
  //// reverse sort so that following algo works
  //std::sort(to_other_tiles.begin(), to_other_tiles.end(), [](const auto& a, const auto& b){return a.n > b.n;} );

  ////--------------------------------------------------
////#ifdef DEBUG
  //// ensure that the array to be removed is unique

  //auto uniq = std::unique( to_other_tiles.begin(), to_other_tiles.end(), [](const auto& a, const auto& b){return a.n == b.n;} );
  //bool contains_duplicate = uniq != to_other_tiles.end();

  //if( contains_duplicate ){
  //  std::cerr << " dupl:";
  //  for(auto& i : to_other_tiles) std::cerr << "," << i.n;
  //  assert(false);
  //}
////#endif
  ////--------------------------------------------------
  //

  //float* locn[3];
  //for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );
  //
  //float* veln[3];
  //for(int i=0; i<3; i++) veln[i] = &( vel(i,0) );
  //
  //int* idn[2];
  //for(int i=0; i<2; i++) idn[i] = &( id(i,0) );
  //
  //
  //// overwrite particles with the last one on the array and 
  //// then resize the array
  //int last = size()-to_other_tiles.size();
  ////std::cout << "del: " << size() << " to be deleted: " << to_other_tiles.size() << std::endl;
  //

  //// NOTE vectorizing the below loop leads to race conditions; hence it is turned off for now
  //  
  ////UniIter::iterate([=] DEVCALLABLE (int ii, ManVec<to_other_tiles_struct> &to_other_tiles){
  //for(int ii=0; ii<to_other_tiles.size(); ii++){
  //  int other = last+ii; //size() - 1 - i;
  //  int indx = to_other_tiles[ii].n;

  //  //if(indx >= last) return;
  //  //std::cout << " sw " << indx << " to " << other << " while last " << last << std::endl;
  //    
  //  //std::cout << "deleting " << indx << " by putting it to " << last << '\n';

  //  // v1 swap routine; more safe
  //  std::swap( locn[0][indx], locn[0][other] );
  //  std::swap( locn[1][indx], locn[1][other] );
  //  std::swap( locn[2][indx], locn[2][other] );

  //  std::swap( veln[0][indx], veln[0][other] );
  //  std::swap( veln[1][indx], veln[1][other] );
  //  std::swap( veln[2][indx], veln[2][other] );

  //  std::swap(  idn[0][indx],  idn[0][other] );
  //  std::swap(  idn[1][indx],  idn[1][other] );

  //  std::swap(  wgtArr[indx],  wgtArr[other] );

  //  // v0: replacement routine
  //  //for(int i=0; i<3; i++) locn[i][indx] = locn[i][other];
  //  //for(int i=0; i<3; i++) veln[i][indx] = veln[i][other];
  //  //for(int i=0; i<2; i++) idn[ i][indx] = idn[ i][other];
  //  //wgtArr[indx] = wgtArr[other]; 

  //}
  ////, to_other_tiles.size(), to_other_tiles);
  //
  //UniIter::sync();
  //
  //// resize if needed and take care of the size
  //last = last < 0 ? 0 : last;
  //if ((last != (int)size()) && (size() > 0)) resize(last);

  ////std::cout << " INFO: " << cid << " v2: removing prtcls :" << Nprtcls << " - " << to_other_tiles.size() << std::endl;
  //Nprtcls -= to_other_tiles.size();


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
      float locx = wrap( neigh.loc(0, n), static_cast<float>(global_mins[0]), static_cast<float>(global_maxs[0]) );
      float locy = wrap( neigh.loc(1, n), static_cast<float>(global_mins[1]), static_cast<float>(global_maxs[1]) );
      float locz = wrap( neigh.loc(2, n), static_cast<float>(global_mins[2]), static_cast<float>(global_maxs[2]) );

      add_identified_particle(
          {locx, locy, locz}, 
          {neigh.vel(0, n), neigh.vel(1, n), neigh.vel(2, n)},
          neigh.wgt(n),
          neigh.id(0,n), neigh.id(1,n));
    }
  }


  //--------------------------------------------------
  // v1

  // count how many prtcls are coming
  //std::vector<int> ind_incoming;  // array for incoming indices
  //ind_incoming.reserve(neigh.to_other_tiles.size() ); //reserve up to max size

  //// pre-loop to collect incoming particles from right directions
  //for(auto&& elem : neigh.to_other_tiles) {

  //  //if(elem.i == 0 && elem.j == 0 && elem.k == 0) continue; 
  //  if (elem.i == -dirs[0] && elem.j == -dirs[1] && elem.k == -dirs[2] ) {
  //    ind_incoming.push_back(elem.n); // add to list
  //  }
  //}

  //int N = this->size(); // current number of prtcls and tip of the storage arrays
  //this->resize( this->size() + ind_incoming.size() ); // pre-reserve the expected num of prtcls; 
  //                                                    // resize of ManVec changes capacity only upwards
  //                                                    // so this is safe to issue here

  ////std::cout << "trnsf_wrap N:" << N << " ind_incoming: " << ind_incoming.size() << "\n";

  //// second loop to add the particles
  //UniIter::iterate([=] DEVCALLABLE (int i, std::vector<int>& ind_incoming){
  //  int ind = ind_incoming[i];
  // 
  //  float locx = wrap( neigh.loc(0, ind), static_cast<float>(global_mins[0]), static_cast<float>(global_maxs[0]) );
  //  float locy = wrap( neigh.loc(1, ind), static_cast<float>(global_mins[1]), static_cast<float>(global_maxs[1]) );
  //  float locz = wrap( neigh.loc(2, ind), static_cast<float>(global_mins[2]), static_cast<float>(global_maxs[2]) );
  // 
  //  // manual insert routine (SIMD friendly)
  //  insert_identified_particle(
  //      {locx, locy, locz}, 
  //      {neigh.vel(0, ind), neigh.vel(1, ind), neigh.vel(2, ind)},
  //      neigh.wgt(ind),
  //      neigh.id(0,ind), neigh.id(1,ind),
  //      N + i );
  // 
  //}, ind_incoming.size(), ind_incoming);

  //Nprtcls = N + ind_incoming.size(); // NOTE: we need to manually remember to expand the array

  //--------------------------------------------------
  // v3: SIMD loops w/o external to_other_tiles storage

  // collect outgoing prtcls
  //ManVec<int> ind_incoming;
  //ind_incoming.reserve(neigh.outgoing_count);
  //for(int ind=0; ind<neigh.size(); ind++){
  //  auto [i,j,k] = info2dir( neigh.infoArr[ind] );
  //  if (i == -dirs[0] && j == -dirs[1] && k == -dirs[2] ) {
  //    ind_incoming.push_back( ind );
  //  }
  //};

  //// insert them to current container
  //int N = this->size(); // current number of prtcls and tip of the storage arrays
  //this->resize( this->size() + ind_incoming.size() ); // pre-reserve the expected num of prtcls; 
  //                                                    // resize of ManVec changes capacity only upwards
  //                                                    // so this is safe to issue here
  //
  //UniIter::iterate([=] DEVCALLABLE (int i, ManVec<int>& ind_incoming){
  //  int ind = ind_incoming[i];
  // 
  //  float locx = (D >= 1) ? wrap( neigh.loc(0, ind), static_cast<float>(global_mins[0]), static_cast<float>(global_maxs[0]) ) : 0.0f;
  //  float locy = (D >= 2) ? wrap( neigh.loc(1, ind), static_cast<float>(global_mins[1]), static_cast<float>(global_maxs[1]) ) : 0.0f;
  //  float locz = (D >= 3) ? wrap( neigh.loc(2, ind), static_cast<float>(global_mins[2]), static_cast<float>(global_maxs[2]) ) : 0.0f;
  // 
  //  // manual insert routine (SIMD friendly)
  //  insert_identified_particle(
  //      {locx, locy, locz}, 
  //      {neigh.vel(0, ind), neigh.vel(1, ind), neigh.vel(2, ind)},
  //      neigh.wgt(ind),
  //      neigh.id(0,ind), neigh.id(1,ind),
  //      N + i );
  // 
  //}, ind_incoming.size(), ind_incoming);
  //Nprtcls = N + ind_incoming.size(); // NOTE: we need to manually remember to expand the array


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

  assert(false); // not workign anymore since new extra_message_size needs info prtcl

  //outgoing_particles.clear();
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  int np = size() + 1;

  // FIXME
  //outgoing_particles.reserve(first_message_size);
  //if(np > first_message_size + extra_message_size) {
  //  std::cerr << "Number of particles in MPI message exceeds maximum message size. See documentation." << std::endl;
  //  exit(1);
  //} else 

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

  //outgoing_particles.clear();
  outgoing_extra_particles.clear();
    

  int num_of_outgoing = 0;
  #pragma omp simd reduction(+:num_of_outgoing)
  for(int n=0; n<size(); n++){
    int info = infoArr[n];

    //auto [i,j,k] = info2dir(info);
    //num_of_outgoing += (i!=0) || (j!=0) || (k!=0 ) ? 1 : 0;
    num_of_outgoing += !is_prtcl_inside(info);
  } 

  // numbers of particles in MPI payloads
  // primary msg is fixed/static and its size is always first_message_size
  //    it carries an info particle (as a first element) that announces total number of prtcls coming
  // secondary (extra) message carries the rest if they do not fit in primary message
  //    it carries an info particle (as a first element) that announces its own length. 
  //    This gives us 2 ways to calculate extra container length: 1) via primary msg info and 2) via secondary info prtcl
    
  int np_tot = num_of_outgoing + 2;  // 2 info particles; one with total num of prtcls, and second with num of extra particles
  int np_extra = np_tot-1 > first_message_size ? np_tot - first_message_size : 1;  // always 1 extra particle 

                                  
  // v0 versions
  //int np = outgoing_count + 1; // value calculated in transfer_and_wrap_particles
  //int np = to_other_tiles.size() + 1;
  //assert(to_other_tiles.size() == outgoing_count ); DEBUG

  //std::cout << "outgoing prtcls: " << outgoing_count << " vs num_of_outgoing: " << num_of_outgoing << " vs to_other_tiles: " << to_other_tiles.size()
  //          << " np_tot: " << np_tot
  //          << " np_extra: " << np_extra
  //          << " reserving1: " << first_message_size
  //          << " reserving2: " << np_extra << "\n";
    
  //if (np_tot - 1 > first_message_size) { // +1 for info particle in primary static msg; -1 for info prtcl in extra payload
  //  // reserve is needed here; if size is less than capacity, we do nothing
  //  outgoing_extra_particles.reserve( np_tot - first_message_size ); 
  //}
  outgoing_extra_particles.reserve( np_extra );


  // first particle is always the message info
  outgoing_particles[0] =            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np_tot, 0}; // store prtcl number in id slot
  outgoing_extra_particles.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np_extra, 0}); // store prtcl number in id slot
                                                                              
  // -------------------------------------------------- 
  // v1 with on-the-fly calculation of escape condition

  // next, pack all other particles
  int ind=1;
  for(int n=0; n<size(); n++){

    // check if moving out
    //auto [i,j,k] = info2dir( infoArr[n] );
    //bool inside_tile = (i == 0) && (j == 0) && (k == 0);
    //bool to_be_packed = !inside_tile;
    bool to_be_packed = !is_prtcl_inside(infoArr[n]);

    // pack if true; split between fixed primary and adaptive extra message containers
    if(to_be_packed) {
      if(ind < first_message_size) {
        //outgoing_particles.push_back({ 
        outgoing_particles[ind] = {
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
      ind++;
    }
  }
  ind++; // for extra prtcl payload

  //std::cout << " inserted " << ind << "\n";

  // -------------------------------------------------- 
  // v0 with to_other_tiles 
  //int i = 1;
  //for (size_t ii = 0; ii < to_other_tiles.size(); ii++)
  //{
  //  const auto &elem = to_other_tiles[ii];
  //  int ind = elem.n;

  //  if(i < first_message_size) {
  //    outgoing_particles.push_back({ 
  //      loc(0, ind), loc(1, ind), loc(2, ind), 
  //      vel(0, ind), vel(1, ind), vel(2, ind), 
  //      wgt(ind), 
  //      id(0, ind), id(1, ind) });
  //  } else {
  //    outgoing_extra_particles.push_back({ 
  //      loc(0, ind), loc(1, ind), loc(2, ind), 
  //      vel(0, ind), vel(1, ind), vel(2, ind), 
  //      wgt(ind), 
  //      id(0, ind), id(1, ind) });
  //  }

  //  i++;
  //}


  // -------------------------------------------------- 
  //outgoing_extra_particles.shrink_to_fit();

  // TODO: set next message size dynamically according to history
  //first_message_size = np;
    
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

  float locx, locy, locz, velx, vely, velz, wgts;
  int ids, proc;

  // get real number of incoming particles
  int number_of_incoming_particles = incoming_particles[0].id; // number stored in id slot
  int number_of_incoming_extra_particles = incoming_extra_particles[0].id;

  int number_of_primary_particles = 
    number_of_incoming_particles > first_message_size 
    ? first_message_size : number_of_incoming_particles;

  int number_of_secondary_particles = number_of_incoming_extra_particles;

#ifdef DEBUG
  assert(number_of_secondary_particles == incoming_extra_particles.size() );
#endif

  // reserve arrays
  int N = size();

  //std::cout 
  //  << " mpi inc1: " << number_of_incoming_particles 
  //  << " mpi inc2: " << number_of_incoming_extra_particles 
  //  << " mpi msg2 size " << incoming_extra_particles.size()
  //  << " prim:" << number_of_primary_particles 
  //  << " seco:" << number_of_secondary_particles 
  //  << " N:"    << N << "\n"; 

  //reserve( N + number_of_incoming_particles + number_of_secondary_particles );  //reserve for addition
  resize( N + number_of_incoming_particles - 2 );  // resize for insertion

  // skipping 1st info particle
  for(int i=1; i<number_of_primary_particles; i++){
    //std::cout << "inserting1 to slot" << N+i-1 << " out of " << N << "\n";
    locx = incoming_particles[i].x;
    locy = incoming_particles[i].y;
    locz = incoming_particles[i].z;

    velx = incoming_particles[i].ux;
    vely = incoming_particles[i].uy;
    velz = incoming_particles[i].uz;
    wgts = incoming_particles[i].w;

    ids  = incoming_particles[i].id;
    proc = incoming_particles[i].proc;

    //add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc);
    insert_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc, N+i-1 );
  }

  N += number_of_primary_particles-1;
  for(int i=1; i<number_of_secondary_particles; i++){
    //std::cout << "inserting2 to slot" << N+i-1 << " out of " << N << "\n";

    locx = incoming_extra_particles[i].x;
    locy = incoming_extra_particles[i].y;
    locz = incoming_extra_particles[i].z;

    velx = incoming_extra_particles[i].ux;
    vely = incoming_extra_particles[i].uy;
    velz = incoming_extra_particles[i].uz;
    wgts = incoming_extra_particles[i].w;

    ids  = incoming_extra_particles[i].id;
    proc = incoming_extra_particles[i].proc;
      
    //add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc);
    insert_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgts, ids, proc, N+i-1);
  }

  // update internal counter after insertion
  Nprtcls += number_of_incoming_particles-2;


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

  //std::cout << "apply_permutation:" << indices.size() << " " << size() << std::endl;

  // check that sizes match
  assert( indices.size() == size() );

  // https://stackoverflow.com/questions/67751784/how-to-do-in-place-sorting-a-list-according-to-a-given-index-in-c
  // and
  // https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095
  for (size_t i=0; i<indices.size(); i++) {
    auto current = i;

    while (i != indices[current]) {
      auto next = indices[current];

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
