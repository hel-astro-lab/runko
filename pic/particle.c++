
#include "particle.h"
#include "../tools/wrap.h"

#include <algorithm>
#include <map>
#include <utility>
#include <mpi.h>
#include <functional>


namespace pic {

inline Particle::Particle(
    real_prtcl x, real_prtcl y, real_prtcl z,
    real_prtcl ux, real_prtcl uy, real_prtcl uz, 
    real_prtcl wgt,
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
  data[0] = static_cast<real_prtcl>(number_of_particles);
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
  locArr.resize(3);
  velArr.resize(3);
  indArr.resize(2);

  // Get the number of processes
  //MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);

  // Get the rank of the process
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  incoming_particles.resize(optimal_message_size);
  incoming_extra_particles.resize(optimal_message_size);

  outgoing_particles.resize(optimal_message_size);
  outgoing_extra_particles.resize(optimal_message_size);
}


template<std::size_t D>
void ParticleContainer<D>::reserve(size_t N) {

  // always reserve at least 1 element to ensure proper array initialization
  if (N <= 0) N = 1;

  for(size_t i=0; i<3; i++) locArr[i].reserve(N);
  for(size_t i=0; i<3; i++) velArr[i].reserve(N);
  for(size_t i=0; i<2; i++) indArr[i].reserve(N);
  wgtArr.reserve(N);
    
  // reserve 1d N x D array for particle-specific fields
  Epart.reserve(N*3);
  Bpart.reserve(N*3);
}

template<std::size_t D>
void ParticleContainer<D>::resize(size_t N)
{
  for(size_t i=0; i<3; i++) locArr[i].resize(N);
  for(size_t i=0; i<3; i++) velArr[i].resize(N);
  for(size_t i=0; i<2; i++) indArr[i].resize(N);
  wgtArr.resize(N);
  Nprtcls = N;
}

template<std::size_t D>
void ParticleContainer<D>::shrink_to_fit()
{
  for(size_t i=0; i<3; i++) locArr[i].shrink_to_fit();
  for(size_t i=0; i<3; i++) velArr[i].shrink_to_fit();
  for(size_t i=0; i<2; i++) indArr[i].shrink_to_fit();
  wgtArr.shrink_to_fit();
}


template<std::size_t D>
size_t ParticleContainer<D>::size() 
{ 
  assert(locArr[0].size() == Nprtcls);
  assert(locArr[1].size() == Nprtcls);
  assert(locArr[2].size() == Nprtcls);

  assert(velArr[0].size() == Nprtcls);
  assert(velArr[1].size() == Nprtcls);
  assert(velArr[2].size() == Nprtcls);

  assert(indArr[0].size() == Nprtcls);
  assert(indArr[1].size() == Nprtcls);

  assert(wgtArr.size() == Nprtcls);

  return Nprtcls; 
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
    std::vector<real_prtcl> prtcl_loc,
    std::vector<real_prtcl> prtcl_vel,
    real_prtcl prtcl_wgt)
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
    std::vector<real_prtcl> prtcl_loc,
    std::vector<real_prtcl> prtcl_vel,
    real_prtcl prtcl_wgt,
    int _id, int _proc)
{
  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  indArr[0].push_back(_id);
  indArr[1].push_back(_proc);

  Nprtcls++;
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
  real_prtcl* locn[3];
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

    if ( (i != 0) || (j != 0) || (k != 0) ) 
      to_other_tiles.insert( std::make_pair( std::make_tuple(i,j,k), n) );
  }
}



template<>
void ParticleContainer<3>::check_outgoing_particles(
    std::array<double,3>& mins,
    std::array<double,3>& maxs)
{
  to_other_tiles.clear();

  // unpack limits
  double 
    xmin = mins[0],
    ymin = mins[1],
    zmin = mins[2],

    xmax = maxs[0],
    ymax = maxs[1],
    zmax = maxs[2];

  int lenx = static_cast<int>( xmax - xmin );
  int leny = static_cast<int>( ymax - ymin );
  int lenz = static_cast<int>( zmax - zmin );

  int i0, j0, k0;

  // shortcut for particle locations
  real_prtcl* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  int i,j,k; // relative indices
  for(size_t n=0; n<size(); n++) {
    i = 0;
    j = 0;
    k = 0;

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
      to_other_tiles.insert( std::make_pair( std::make_tuple(i,j,k), n) );
  }
}

template<std::size_t D>
void ParticleContainer<D>::delete_particles(std::vector<int> to_be_deleted) 
{
  std::sort(to_be_deleted.begin(), to_be_deleted.end(), std::greater<int>() );

  real_prtcl* locn[3];
  for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  real_prtcl* veln[3];
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

  return;
}


template<std::size_t D>
void ParticleContainer<D>::delete_transferred_particles()
{
  std::vector<int> to_be_deleted;

  // get transferred 
  for(auto& elem : to_other_tiles) to_be_deleted.push_back( elem.second );

  delete_particles(to_be_deleted);
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

  // particle overflow from tiles is done in shortest precision
  // to avoid rounding off errors and particles left in a limbo
  // between tiles.
  real_prtcl locx, locy, locz, velx, vely, velz, wgt;
  int id, proc;

  int i;
  for (auto&& elem : neigh.to_other_tiles) {
    if(std::get<0>(elem.first) == 0 && std::get<1>(elem.first) == 0) continue; 

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile

    if (std::get<0>(elem.first) == -dirs[0] &&
        std::get<1>(elem.first) == -dirs[1] ) {

      i = elem.second;

      locx = wrap( neigh.loc(0, i), static_cast<real_prtcl>(global_mins[0]), static_cast<real_prtcl>(global_maxs[0]) );
      locy = wrap( neigh.loc(1, i), static_cast<real_prtcl>(global_mins[1]), static_cast<real_prtcl>(global_maxs[1]) );
      locz = wrap( neigh.loc(2, i), static_cast<real_prtcl>(global_mins[2]), static_cast<real_prtcl>(global_maxs[2]) );

      velx = neigh.vel(0, i);
      vely = neigh.vel(1, i);
      velz = neigh.vel(2, i);

      wgt  = neigh.wgt(i);

      id   = neigh.id(0,i);
      proc = neigh.id(1,i);

      add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgt, id, proc);
    }
  }

  return;
}


template<>
void ParticleContainer<3>::transfer_and_wrap_particles( 
    ParticleContainer& neigh,
    std::array<int,3>    dirs, 
    std::array<double,3>& global_mins, 
    std::array<double,3>& global_maxs
    )
{

  // particle overflow from tiles is done in shortest precision
  // to avoid rounding off errors and particles left in a limbo
  // between tiles.
  real_prtcl locx, locy, locz, velx, vely, velz, wgt;
  int id, proc;

  int i;
  for (auto&& elem : neigh.to_other_tiles) {
      
    if(std::get<0>(elem.first) == 0 && 
       std::get<1>(elem.first) == 0 &&
       std::get<2>(elem.first) == 0) continue; 

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile

    if (std::get<0>(elem.first) == -dirs[0] &&
        std::get<1>(elem.first) == -dirs[1] &&
        std::get<2>(elem.first) == -dirs[2] ) {

      i = elem.second;

      locx = wrap( neigh.loc(0, i), static_cast<real_prtcl>(global_mins[0]), static_cast<real_prtcl>(global_maxs[0]) );
      locy = wrap( neigh.loc(1, i), static_cast<real_prtcl>(global_mins[1]), static_cast<real_prtcl>(global_maxs[1]) );
      locz = wrap( neigh.loc(2, i), static_cast<real_prtcl>(global_mins[2]), static_cast<real_prtcl>(global_maxs[2]) );

      velx = neigh.vel(0, i);
      vely = neigh.vel(1, i);
      velz = neigh.vel(2, i);

      wgt  = neigh.wgt(i);

      id   = neigh.id(0,i);
      proc = neigh.id(1,i);

      add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgt, id, proc);
    }
  }

  return;
}

template<std::size_t D>
void ParticleContainer<D>::pack_all_particles()
{
  outgoing_particles.clear();
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  int np = size() + 1;

  outgoing_particles.reserve(optimal_message_size);
  if(np-optimal_message_size > 0) {
    outgoing_extra_particles.reserve( np-optimal_message_size );
  }

  // first particle is always the message info
  outgoing_particles.emplace_back(np);

  // next, pack all other particles
  int i=1;
  for(size_t ind=0; ind < size(); ind++) {
    if(i < optimal_message_size) {
      outgoing_particles.emplace_back( 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) );
    } else {
      outgoing_extra_particles.emplace_back( 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) );
    }
    i++;
  }

  //outgoing_extra_particles.shrink_to_fit();
}


template<std::size_t D>
void ParticleContainer<D>::pack_outgoing_particles()
{
  outgoing_particles.clear();
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  int np = to_other_tiles.size() + 1;

  outgoing_particles.reserve(optimal_message_size);
  if (np-optimal_message_size > 0) {
    outgoing_extra_particles.reserve( np-optimal_message_size);
  }

  // first particle is always the message info
  outgoing_particles.emplace_back(np);

  // next, pack all other particles
  int i=1, ind;
  for (auto&& elem : to_other_tiles) {
    ind = elem.second;

    if(i < optimal_message_size) {
      outgoing_particles.emplace_back( 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) );
    } else {
      outgoing_extra_particles.emplace_back( 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind), 
        id(0, ind), id(1, ind) );
    }

    i++;
  }

  //outgoing_extra_particles.shrink_to_fit();

  // TODO: set next message size dynamically according to history
  //optimal_message_size = np;

  return;
}


template<std::size_t D>
void ParticleContainer<D>::unpack_incoming_particles()
{
  real_prtcl locx, locy, locz, velx, vely, velz, wgts;
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

  return;
}


template<std::size_t D>
void ParticleContainer<D>::set_keygen_state(int __key, int __rank)
{
  _key  = __key;
  _rank = __rank;
}


} // end ns pic


template class pic::ParticleContainer<1>;
template class pic::ParticleContainer<2>;
template class pic::ParticleContainer<3>;
