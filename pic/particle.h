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
  std::array<float_tp, 7> data = {-1,-2,-3,-4,-5,-7};

  /// particle id
  int _id = 0;

  /// mpi rank separator
  int _proc = 0;

  Particle() {};

  /// standard ctor
  Particle(float_tp x,  float_tp y,  float_tp z,
           float_tp ux, float_tp uy, float_tp uz,
           float_tp wgt,
           int __ind, int __proc
           );

  /// special ctor for info prtcl
  Particle(size_t number_of_particles);

  inline float_tp& x()   { return data[0]; };
  inline float_tp& y()   { return data[1]; };
  inline float_tp& z()   { return data[2]; };
  inline float_tp& ux()  { return data[3]; };
  inline float_tp& uy()  { return data[4]; };
  inline float_tp& uz()  { return data[5]; };
  inline float_tp& wgt() { return data[6]; };

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

  std::vector< std::vector<float_tp> > locArr;
  std::vector< std::vector<float_tp> > velArr;
  std::vector< std::vector<int> > indArr;
  std::vector<float_tp> wgtArr;

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
  std::vector<float_tp> Epart;

  //! particle specific magnetic field components
  std::vector<float_tp> Bpart;

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

  // "shrink to fit" all internal main containers
  virtual void shrink_to_fit();

  /// size of the container (in terms of particles)
  size_t size();


  //--------------------------------------------------
  // locations
  virtual inline float_tp loc( size_t idim, size_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  virtual inline float_tp& loc( size_t idim, size_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  virtual inline std::vector<float_tp> loc(size_t idim) const 
  {
    return locArr[idim];
  }

  virtual inline std::vector<float_tp>& loc(size_t idim)
  {
    return locArr[idim];
  }

  //--------------------------------------------------
  // velocities
  virtual inline float_tp vel( size_t idim, size_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  virtual inline float_tp& vel( size_t idim, size_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  virtual inline std::vector<float_tp> vel(size_t idim) const 
  {
    return velArr[idim];
  }

  virtual inline std::vector<float_tp>& vel(size_t idim)
  {
    return velArr[idim];
  }

  //--------------------------------------------------
  // weights
  virtual inline float_tp wgt( size_t iprtcl ) const
  {
    return wgtArr[iprtcl];
  }

  virtual inline float_tp& wgt( size_t iprtcl )       
  {
    return wgtArr[iprtcl];
  }

  virtual inline std::vector<float_tp> wgt() const
  {
    return wgtArr;
  }

  virtual inline std::vector<float_tp>& wgt()
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

  virtual inline std::vector<int>& id(size_t idim)
  {
    return indArr[idim];
  }

  // particle creation
  virtual void add_particle (
      std::vector<float_tp> prtcl_loc,
      std::vector<float_tp> prtcl_vel,
      float_tp prtcl_wgt);

  // particle creation
  virtual void add_identified_particle (
      std::vector<float_tp> prtcl_loc,
      std::vector<float_tp> prtcl_vel,
      float_tp prtcl_wgt, 
      int _ind, int _proc);



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
