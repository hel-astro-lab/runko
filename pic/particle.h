#pragma once

#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <cassert>
#include <iostream>

#include "../definitions.h"

#include "../tools/iter/dynArray.h"



namespace pic {


/// Particle class for easier communication
class Particle
{
public:

  /// actual particle data
  std::array<real_prtcl, 7> data = {-1,-2,-3,-4,-5,-7};

  /// particle id
  int _id = 0;

  /// mpi rank separator
  int _proc = 0;

  Particle() = default;

  /// standard ctor
  Particle(real_prtcl x,  real_prtcl y,  real_prtcl z,
           real_prtcl ux, real_prtcl uy, real_prtcl uz,
           real_prtcl wgt,
           int __ind, int __proc
           );

  /// special ctor for info prtcl
  Particle(size_t number_of_particles);

  inline real_prtcl& x()   { return data[0]; };
  inline real_prtcl& y()   { return data[1]; };
  inline real_prtcl& z()   { return data[2]; };
  inline real_prtcl& ux()  { return data[3]; };
  inline real_prtcl& uy()  { return data[4]; };
  inline real_prtcl& uz()  { return data[5]; };
  inline real_prtcl& wgt() { return data[6]; };

  inline int& id()   { return _id;   };
  inline int& proc() { return _proc; };

  virtual ~Particle() = default;

  /// special method for info particle that re-uses x mem location
  size_t number_of_particles();

};


struct to_other_tiles_struct{
  //
  int i;
  int j;
  int k;
  int n;
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
template<std::size_t D>
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

  std::array<DevVec<real_prtcl>, 3 > locArr;
  std::array<DevVec<real_prtcl>, 3 > velArr;
  std::array<DevVec<int>, 2 > indArr;
  std::vector<real_prtcl> wgtArr;

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
  DevVec<real_prtcl> Epart;

  //! particle specific magnetic field components
  DevVec<real_prtcl> Bpart;

  //! multimap of particles going to other tiles
  using mapType = DevVec<to_other_tiles_struct>;
  mapType to_other_tiles;

  // normalization factor
  double q = 1.0; 

  /// Constructor 
  ParticleContainer();

  // default virtual dtor
  virtual ~ParticleContainer() = default;


  //--------------------------------------------------
    
  /// reserve memory for particles
  virtual void reserve(size_t N);

  // resize everything
  virtual void resize(size_t N);

  // "shrink to fit" all internal main containers
  virtual void shrink_to_fit();

  /// size of the container (in terms of particles)
  size_t size();


  //--------------------------------------------------
  // locations
  virtual inline real_prtcl loc( size_t idim, size_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  virtual inline real_prtcl& loc( size_t idim, size_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  virtual inline std::vector<real_prtcl> loc(size_t idim) const 
  {
    return locArr[idim].toVector();
  }

/*
  virtual inline std::vector<real_prtcl>& loc(size_t idim)
  {
    return locArr[idim];
  }
*/
  //--------------------------------------------------
  // velocities
  virtual inline real_prtcl vel( size_t idim, size_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  virtual inline real_prtcl& vel( size_t idim, size_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  virtual inline std::vector<real_prtcl> vel(size_t idim) const 
  {
    return velArr[idim].toVector();
  }
/*
  virtual inline std::vector<real_prtcl>& vel(size_t idim)
  {
    return velArr[idim];
  }
*/
  //--------------------------------------------------
  // weights
  virtual inline real_prtcl wgt( size_t iprtcl ) const
  {
    return wgtArr[iprtcl];
  }

  virtual inline real_prtcl& wgt( size_t iprtcl )       
  {
    return wgtArr[iprtcl];
  }

  virtual inline std::vector<real_prtcl> wgt() const
  {
    return wgtArr;
  }

/*
  virtual inline std::vector<real_prtcl>& wgt()
  {
    return wgtArr;
  }
*/
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
    return indArr[idim].toVector();
  }

/*
  virtual inline std::vector<int>& id(size_t idim)
  {
    return indArr[idim];
  }
*/
  // particle creation
  virtual void add_particle (
      std::vector<real_prtcl> prtcl_loc,
      std::vector<real_prtcl> prtcl_vel,
      real_prtcl prtcl_wgt);

  // particle creation
  virtual void add_identified_particle (
      std::vector<real_prtcl> prtcl_loc,
      std::vector<real_prtcl> prtcl_vel,
      real_prtcl prtcl_wgt, 
      int _id, int _proc);



  // --------------------------------------------------
  // particle boundary checks

  /// check and mark particles exceeding given limits
  //template <size_t D> 
  void check_outgoing_particles(
      std::array<double,3>&,
      std::array<double,3>& );


  /// delete particles that went beyond boundaries, i.e.,
  // ended up in to_other_tiles box
  void delete_transferred_particles();

  /// process through an index list and delete particles in it
  void delete_particles(std::vector<int> to_be_deleted);

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
