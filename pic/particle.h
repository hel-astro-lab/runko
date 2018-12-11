#pragma once

#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <cassert>
#include <iostream>

#include "../definitions.h"



namespace pic {


/// Particle class for easier communication
class Particle
{
public:

  std::array<double,7> data;

  Particle() {};

  Particle(double x,  double y,  double z,
           double ux, double uy, double uz,
           double wgt);

  inline double& x()   { return data[0]; };
  inline double& y()   { return data[1]; };
  inline double& z()   { return data[2]; };
  inline double& ux()  { return data[3]; };
  inline double& uy()  { return data[4]; };
  inline double& uz()  { return data[5]; };
  inline double& wgt() { return data[6]; };
};


/// Special handling of particle MPI message info 
// via this auxiliary helper class
class InfoParticle : public Particle
{
public:

  InfoParticle(size_t np) 
    : Particle(
        static_cast<double>(np), 0,0,
        0,0,0,
        0) {}

  InfoParticle(Particle& prtcl) {
    data[0] = prtcl.x();
  }

  size_t size() {return static_cast<size_t>(Particle::x());}

private:
  using Particle::x;
  using Particle::y;
  using Particle::z;
  using Particle::ux;
  using Particle::uy;
  using Particle::uz;
  using Particle::wgt;

};




/*! \brief Container of particles inside the tile
*
* Container to hold plasma particles 
*
* includes:
*   pos/loc
*   vel
*   wgt
*   qm
*
*/
class ParticleContainer {

  protected:

  size_t Nprtcls;

  std::vector< std::vector<double> > locArr;
  std::vector< std::vector<double> > velArr;
  std::vector<double> wgtArr;

  public:
    
  /// packed outgoing particles
  std::vector<Particle> outgoing_particles;
  std::vector<Particle> outgoing_extra_particles;

  /// pack particles that are marked as outflowing
  void pack_outgoing_particles();

  /// packed incoming particles
  std::vector<Particle> incoming_particles;
  std::vector<Particle> incoming_extra_particles;

  /// unpack incoming particles into internal vectors
  void unpack_incoming_particles();

  /// dynamic message size that traces the optimal
  // message length (i.e., number of particles) hand 
  // in hand with the corresponding receiver side.
  int optimal_message_size = 10;

  //! particle specific electric field components
  std::vector<double> Epart;

  //! particle specific magnetic field components
  std::vector<double> Bpart;

  //! multimap of particles going to other tiles
  typedef std::multimap<std::tuple<int,int,int>, int> mapType;
  mapType to_other_tiles;


  // size of the internal mesh
  size_t Nx;
  size_t Ny;
  size_t Nz;

  double q = 1.0; // normalization factor


  /// Constructor 
  ParticleContainer();

  // default virtual dtor
  virtual ~ParticleContainer() = default;


  //--------------------------------------------------
    
  /// reserve memory for particles
  virtual void reserve(size_t N);

  // resize everything
  virtual void resize(size_t N);

  /// size of the container (in terms of particles)
  size_t size();


  //--------------------------------------------------
  // locations
  virtual inline double loc( size_t idim, size_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  virtual inline double& loc( size_t idim, size_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  virtual inline std::vector<double> loc(size_t idim) const 
  {
    return locArr[idim];
  }


  //--------------------------------------------------
  // velocities
  virtual inline double vel( size_t idim, size_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  virtual inline double& vel( size_t idim, size_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  virtual inline std::vector<double> vel(size_t idim) const 
  {
    return velArr[idim];
  }

  //--------------------------------------------------
  // weights
  virtual inline double wgt( size_t iprtcl ) const
  {
    return wgtArr[iprtcl];
  }

  virtual inline double& wgt( size_t iprtcl )       
  {
    return wgtArr[iprtcl];
  }

  virtual inline std::vector<double>& wgt() 
  {
    return wgtArr;
  }

  // particle creation
  virtual void add_particle (
      std::vector<double> prtcl_loc,
      std::vector<double> prtcl_vel,
      double prtcl_wgt);

  // --------------------------------------------------
  // particle boundary checks

  /// check and mark particles exceeding given limits
  void check_outgoing_particles(
      std::array<double,3>&,
      std::array<double,3>& );


  /// delete particles that went beyond boundaries, i.e.,
  // ended up in to_other_tiles box
  void delete_transferred_particles();


  /// transfer particles between blocks
  void transfer_and_wrap_particles(
      ParticleContainer&, 
      std::array<int,3>,
      std::array<double,3>&,
      std::array<double,3>&);

};





} // end of namespace pic
