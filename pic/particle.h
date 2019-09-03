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

  /// actual particle data
  std::array<double,7> data = {-1,-2,-3,-4,-5,-7};

  /// particle id
  int _id = 0;

  /// mpi rank separator
  int _proc = 0;

  Particle() {};

  /// standard ctor
  Particle(double x,  double y,  double z,
           double ux, double uy, double uz,
           double wgt,
           int __ind, int __proc
           );

  /// special ctor for info prtcl
  Particle(size_t number_of_particles);

  inline double& x()   { return data[0]; };
  inline double& y()   { return data[1]; };
  inline double& z()   { return data[2]; };
  inline double& ux()  { return data[3]; };
  inline double& uy()  { return data[4]; };
  inline double& uz()  { return data[5]; };
  inline double& wgt() { return data[6]; };

  inline int& id()   { return _id;   };
  inline int& proc() { return _proc; };

  virtual ~Particle() = default;

  /// special method for info particle that re-uses x mem location
  size_t number_of_particles();

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

  private:

  /// mpi rank
  int _rank = 0;

  /// running key generator seed
  int _key = 0;

  /// unique key generator
  std::pair<int,int> keygen();


  protected:

  size_t Nprtcls = 0;

  std::vector< std::vector<double> > locArr;
  std::vector< std::vector<double> > velArr;
  std::vector< std::vector<int> > indArr;
  std::vector<double> wgtArr;

  public:
    
  /// packed outgoing particles
  std::vector<Particle> outgoing_particles;
  std::vector<Particle> outgoing_extra_particles;

  /// pack all particles in the container
  void pack_all_particles();

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
  int optimal_message_size = 3000;

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

  virtual inline std::vector<double> wgt() const
  {
    return wgtArr;
  }


  //--------------------------------------------------
  // id
  virtual inline int id( size_t idim, size_t iprtcl ) const
  {
    return indArr[idim][iprtcl];
  }

  virtual inline int& id( size_t idim, size_t iprtcl )       
  {
    return indArr[idim][iprtcl];
  }

  virtual inline std::vector<int> id(size_t idim) const 
  {
    return indArr[idim];
  }


  // particle creation
  virtual void add_particle (
      std::vector<double> prtcl_loc,
      std::vector<double> prtcl_vel,
      double prtcl_wgt);

  // particle creation
  virtual void add_identified_particle (
      std::vector<double> prtcl_loc,
      std::vector<double> prtcl_vel,
      double prtcl_wgt, 
      int _ind, int _proc);



  // --------------------------------------------------
  // particle boundary checks

  /// check and mark particles exceeding given limits
  void check_outgoing_particles(
      std::array<double,3>&,
      std::array<double,3>& );


  /// delete particles that went beyond boundaries, i.e.,
  // ended up in to_other_tiles box
  void delete_transferred_particles();

  /// process through an index list and delete particles in it
  void delete_particles(std::vector<int> l);

  /// transfer particles between blocks
  void transfer_and_wrap_particles(
      ParticleContainer&, 
      std::array<int,3>,
      std::array<double,3>&,
      std::array<double,3>&);


  /// set keygenerator state
  void set_keygen_state(int __key, int __rank);

};




} // end of namespace pic
