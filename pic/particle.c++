
#include "particle.h"
#include "../tools/wrap.h"

#include <algorithm>
#include <map>


namespace pic {

inline Particle::Particle(
    double x, double y, double z,
    double ux, double uy, double uz, 
    double wgt)
{
  data[0] = x;
  data[1] = y;
  data[2] = z;
  data[3] = ux;
  data[4] = uy;
  data[5] = uz;
  data[6] = wgt;
}


ParticleContainer::ParticleContainer()
{ 
  locArr.resize(3);
  velArr.resize(3);
};


void ParticleContainer::reserve(size_t N) {

  // always reserve at least 1 element to ensure proper array initialization
  if (N <= 0) N = 1;

  for(size_t i=0; i<3; i++) locArr[i].reserve(N);
  for(size_t i=0; i<3; i++) velArr[i].reserve(N);
  wgtArr.reserve(N);
    
  // reserve 1d N x D array for particle-specific fields
  Epart.reserve(N*3);
  Bpart.reserve(N*3);
}

void ParticleContainer::resize(size_t N)
{
  for(size_t i=0; i<3; i++) locArr[i].resize(N);
  for(size_t i=0; i<3; i++) velArr[i].resize(N);
  wgtArr.resize(N);
  Nprtcls = N;
}


size_t ParticleContainer::size() 
{ 
  // FIXME: these fail
  //std::cout << "prtcl container size()" << Nprtcls << " vs " 
  //  << locArr[0].size() << "/"
  //  << locArr[1].size() << "/"
  //  << locArr[2].size() << "/"
  //  << velArr[0].size() << "/"
  //  << velArr[1].size() << "/"
  //  << velArr[2].size() << "/"
  //  << wgtArr.size() << "\n";

  assert(locArr[0].size() == Nprtcls);
  assert(locArr[1].size() == Nprtcls);
  assert(locArr[2].size() == Nprtcls);

  assert(velArr[0].size() == Nprtcls);
  assert(velArr[1].size() == Nprtcls);
  assert(velArr[2].size() == Nprtcls);

  assert(wgtArr.size() == Nprtcls);

  //return locArr[0].size();
  return Nprtcls; // FIXME: this is the correct way to return
}


void ParticleContainer::add_particle (
    std::vector<double> prtcl_loc,
    std::vector<double> prtcl_vel,
    double prtcl_wgt)
{

  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  Nprtcls++;
}


void ParticleContainer::check_outgoing_particles(
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

  // shortcut for particle locations
  double* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  double x0, y0, z0;

  int i,j,k; // relative indices
  for(size_t n=0; n<size(); n++) {
    i = 0;
    j = 0;
    k = 0;

    x0 = locn[0][n];
    y0 = locn[1][n];
    z0 = locn[2][n];

    if(x0 <  xmin) i--; // left wrap
    if(x0 >= xmax) i++; // right wrap

    if(y0 <  ymin) j--; // bottom wrap
    if(y0 >= ymax) j++; // top wrap

    if(z0 <  zmin) k--; // back
    if(z0 >= zmax) k++; // front

    // FIXME: hack to make this work with 2D 
    if ((i == 0) && (j == 0)) continue; 

    if ( (i != 0) || (j != 0) || (k != 0) ) 
      to_other_tiles.insert( std::make_pair( std::make_tuple(i,j,k), n) );
  }
}


void ParticleContainer::delete_transferred_particles()
{
  std::vector<int> to_be_deleted;

  // get and reverse sort deleted particles
  for(auto& elem : to_other_tiles) to_be_deleted.push_back( elem.second );
  std::sort(to_be_deleted.begin(), to_be_deleted.end(), std::greater<int>() );

  double* locn[3];
  for(int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  double* veln[3];
  for(int i=0; i<3; i++) veln[i] = &( vel(i,0) );

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
  }

  // resize if needed and take care of the size
  last = last < 0 ? 0 : last;
  if ((last != (int)size()) && (size() > 0)) resize(last);

  return;
}


void ParticleContainer::transfer_and_wrap_particles( 
    ParticleContainer& neigh,
    std::array<int,3>    dirs, 
    std::array<double,3>& global_mins, 
    std::array<double,3>& global_maxs
    )
{
  double locx, locy, locz, velx, vely, velz, wgt;
  int ind;

  for (auto&& elem : neigh.to_other_tiles) {
      
    //TODO: collapsed z-dimension due to 2D corgi tiles
    if (std::get<0>(elem.first) == 0 &&
        std::get<1>(elem.first) == 0 ) 
    { continue; }

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile
    if (std::get<0>(elem.first) == -dirs[0] &&
        std::get<1>(elem.first) == -dirs[1] ) {

      ind = elem.second;

      locx = wrap( neigh.loc(0, ind), global_mins[0], global_maxs[0] );
      locy = wrap( neigh.loc(1, ind), global_mins[1], global_maxs[1] );
      locz = wrap( neigh.loc(2, ind), global_mins[2], global_maxs[2] );

      velx = neigh.vel(0, ind);
      vely = neigh.vel(1, ind);
      velz = neigh.vel(2, ind);

      wgt  = neigh.wgt(ind);

      add_particle({locx,locy,locz}, {velx,vely,velz}, wgt);
    }
  }

  return;
}

void ParticleContainer::pack_outgoing_particles()
{
  outgoing_particles.clear();
  outgoing_extra_particles.clear();
    
  // +1 for info particle
  int np = to_other_tiles.size() + 1;
  InfoParticle infoprtcl(np);

  //if (np>1) {
  //  std::cout << "Packing Np:" << np << " and extra is: " << np-optimal_message_size << "\n";
  //}

  outgoing_particles.reserve(optimal_message_size);
  if (np-optimal_message_size > 0) {
    //std::cout << "EXTRA send with " << np-optimal_message_size << "\n";
    outgoing_extra_particles.reserve( np-optimal_message_size);
  }

  // first particle is always the message info
  outgoing_particles.push_back(infoprtcl);

  // next, pack all other particles
  int i=1, ind;
  for (auto&& elem : to_other_tiles) {
    ind = elem.second;

    if(i < optimal_message_size) {
      outgoing_particles.emplace_back( 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind));
    } else {
      outgoing_extra_particles.emplace_back( 
        loc(0, ind), loc(1, ind), loc(2, ind), 
        vel(0, ind), vel(1, ind), vel(2, ind), 
        wgt(ind));
    }

    i++;
  }

  //std::cout << " outg arr size:" << outgoing_particles.size()
  //          << " outgE arr size: " << outgoing_extra_particles.size()
  //          << "\n";

  // TODO: set next message size dynamically according to history
  //optimal_message_size = np;

  return;
}


void ParticleContainer::unpack_incoming_particles()
{
  double locx, locy, locz, velx, vely, velz, wgt;

  // get real number of incoming particles
  InfoParticle msginfo(incoming_particles[0]);
  int number_of_incoming_particles = msginfo.size();

  int number_of_primary_particles = 
    number_of_incoming_particles > optimal_message_size 
    ? optimal_message_size : number_of_incoming_particles;

  int number_of_secondary_particles = incoming_extra_particles.size();

  //if (number_of_incoming_particles > 1) {
  //  std::cout << "unpacking Np:" << number_of_incoming_particles
  //    << " /+ " << number_of_primary_particles << "/"
  //    << number_of_secondary_particles << "\n";
  //}

  // skipping 1st info particle
  for(int i=1; i<number_of_primary_particles; i++){
    locx = incoming_particles[i].x();
    locy = incoming_particles[i].y();
    locz = incoming_particles[i].z();

    velx = incoming_particles[i].ux();
    vely = incoming_particles[i].uy();
    velz = incoming_particles[i].uz();
    wgt  = incoming_particles[i].wgt();

    add_particle({locx,locy,locz}, {velx,vely,velz}, wgt);
  }

  for(int i=0; i<number_of_secondary_particles; i++){
    locx = incoming_extra_particles[i].x();
    locy = incoming_extra_particles[i].y();
    locz = incoming_extra_particles[i].z();

    velx = incoming_extra_particles[i].ux();
    vely = incoming_extra_particles[i].uy();
    velz = incoming_extra_particles[i].uz();
    wgt  = incoming_extra_particles[i].wgt();

    add_particle({locx,locy,locz}, {velx,vely,velz}, wgt);
  }

  return;
}





} // end ns pic
