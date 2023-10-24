#pragma once

#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <cassert>
#include <iostream>
#include <string>

#include "../definitions.h"

#include "../tools/iter/dynArray.h"
#include "../tools/iter/allocator.h"
#include "../tools/iter/managed_alloc.h"
#include "../tools/sort.h"



namespace pic {


/// Particle class for easier communication
class Particle
{
public:

  /// actual particle data
  std::array<float_p, 7> data = {-1,-2,-3,-4,-5,-7};

  /// particle id
  int _id = 0;

  /// mpi rank separator
  int _proc = 0;

  Particle() = default;

  /// standard ctor
  Particle(float_p x,  float_p y,  float_p z,
           float_p ux, float_p uy, float_p uz,
           float_p wgt,
           int __ind, int __proc
           );

  /// special ctor for info prtcl
  Particle(int number_of_particles);

  inline float_p& x()   { return data[0]; };
  inline float_p& y()   { return data[1]; };
  inline float_p& z()   { return data[2]; };
  inline float_p& ux()  { return data[3]; };
  inline float_p& uy()  { return data[4]; };
  inline float_p& uz()  { return data[5]; };
  inline float_p& wgt() { return data[6]; };

  inline int& id()   { return _id;   };
  inline int& proc() { return _proc; };

  virtual ~Particle() = default;

  /// special method for info particle that re-uses x loc mem slot
  int number_of_particles();

};

struct to_other_tiles_struct{
  //
  int i;
  int j;
  int k;
  size_t n;
};

/*! \brief Container of particles inside the tile
*
* Container to hold plasma particles . Includes: pos/loc vel wgt qm.
*
*/
template<std::size_t D>
class ParticleContainer{

  private:

  /// mpi rank
  int _rank = 0;

  /// running key generator seed
  int _key = 0;

  /// unique key generator
  std::pair<int,int> keygen();

  protected:

  std::array<ManVec<float_p>, 3 > locArr; // x y z location
  std::array<ManVec<float_p>, 3 > velArr; // vx vy vz velocities
  std::array<ManVec<int>, 2 > indArr;     // cpu,id index
  ManVec<float_p> wgtArr;                 // weight

  public:

  int Nprtcls = 0;
  int cid = 0; // container identification id (minted with tile id)

  // these arrays are required for QED interactions
  ManVec<float_p> wgtCumArr;              // cumulative weights; kept 0 if not needed
  ManVec<float_p> eneArr;                 // particle energies

    
  /// packed outgoing particles
  ManVec<Particle> outgoing_particles;
  ManVec<Particle> outgoing_extra_particles;

  /// pack all particles in the container
  void pack_all_particles();

  /// pack particles that are marked as outflowing
  void pack_outgoing_particles();

  /// packed incoming particles
  ManVec<Particle> incoming_particles;
  ManVec<Particle> incoming_extra_particles;

#ifdef GPU
  // incoming indexes, to optimize transfer_and_wrap_particles for GPUs
  //int incoming_count;
  ManVec<int> incoming_particleIndexes;
  ManVec<int> particleIndexesA;
  ManVec<int> particleIndexesB;
  int pCount;

  void     *d_temp_storage = NULL;
  size_t   temp_storage_bytes = 0;

  // incoming indexes, to optimize transfer_and_wrap_particles for GPUs
  ManVec<int> outgoing_particleIndexes;
#endif

  /// number of particles flowing out from the tile
  int outgoing_count;

  /// unpack incoming particles into internal vectors
  void unpack_incoming_particles();

  // size of MPI particle buffers
  const int first_message_size = 4096; 
  // NOTE maximum prtcl size during first iteration is 2*first_msg; then resized

  //! particle specific electric field components
  ManVec<float_p> Epart;

  //! particle specific magnetic field components
  ManVec<float_p> Bpart;

  //! multimap of particles going to other tiles
  using mapType = ManVec<to_other_tiles_struct>;
  mapType to_other_tiles;

  // particle charge 
  double q = 1.0; 

  // particle mass
  double m = 1.0; 

  // type identifier; default generic type
  std::string type = "gen";

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
  //DEVCALLABLE size_t size() const { return Nprtcls; }
  //DEVCALLABLE size_t size() const { return locArr[0].size(); } // FIXME defaul

  DEVCALLABLE size_t size() const { 

    bool ts[9] = {0,0,0,0,0,0,0,0,0};

    ts[0] = locArr[0].size() == locArr[1].size();
    ts[1] = locArr[0].size() == locArr[2].size();
    ts[2] = locArr[0].size() == velArr[0].size();
    ts[3] = locArr[0].size() == velArr[1].size();
    ts[4] = locArr[0].size() == velArr[2].size();
    ts[5] = locArr[0].size() == indArr[0].size();
    ts[6] = locArr[0].size() == indArr[1].size();
    ts[7] = locArr[0].size() == wgtArr.size();
    ts[8] = locArr[0].size() == Nprtcls;

    //bool ts2 = 0;
    //for(size_t i = 0; i<9; i++) ts2 += !ts[i];

    if( !ts[0] ||
        !ts[1] ||
        !ts[2] ||
        !ts[3] ||
        !ts[4] ||
        !ts[5] ||
        !ts[6] ||
        !ts[7] ||
        !ts[8] ){
      std::cerr << "ERROR: particle number mismatch\n";
      std::cerr << locArr[1].size() << std::endl;;
      std::cerr << locArr[2].size() << std::endl;;
      std::cerr << velArr[0].size() << std::endl;;
      std::cerr << velArr[1].size() << std::endl;;
      std::cerr << velArr[2].size() << std::endl;;
      std::cerr << indArr[0].size() << std::endl;;
      std::cerr << indArr[1].size() << std::endl;;
      std::cerr << wgtArr.size()    << std::endl;;
      std::cerr << Nprtcls          << std::endl;;
      assert(false);
    }

    return Nprtcls; //locArr[0].size(); 
  }

  //--------------------------------------------------
  // locations
  DEVCALLABLE
  inline float_p loc( size_t idim, size_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  DEVCALLABLE
  inline float_p& loc( size_t idim, size_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  inline std::vector<float_p> loc(size_t idim) 
  {
    std::vector<float_p> ret;
    for(const auto& e: locArr[idim])
      ret.push_back(e);
    return ret;//locArr[idim];
  }

  //--------------------------------------------------
  // velocities
  DEVCALLABLE
  inline float_p vel( size_t idim, size_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  DEVCALLABLE
  inline float_p& vel( size_t idim, size_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  inline std::vector<float_p> vel(size_t idim) 
  {
    //return velArr[idim];
    std::vector<float_p> ret;
    for(const auto& e: velArr[idim])
      ret.push_back(e);
    return ret;
  }
/*
  virtual inline std::vector<float_p>& vel(size_t idim)
  {
    return velArr[idim];
  }
*/
  //--------------------------------------------------
  // weights
  DEVCALLABLE
  inline float_p wgt( size_t iprtcl ) const
  {
    return wgtArr[iprtcl];
  }

  DEVCALLABLE
  inline float_p& wgt( size_t iprtcl )       
  {
    return wgtArr[iprtcl];
  }

  inline std::vector<float_p> wgt()
  {
    //return wgtArr;
    std::vector<float_p> ret;
    for(const auto& e: wgtArr)
      ret.push_back(e);
    return ret;
  }

/*
  virtual inline std::vector<float_p>& wgt()
  {
    return wgtArr;
  }
*/
  //--------------------------------------------------
  // id
  DEVCALLABLE
  inline int id( size_t idim, size_t iprtcl ) const
  {
    return indArr[idim][iprtcl];
  }

  DEVCALLABLE
  inline int& id( size_t idim, size_t iprtcl )       
  {
    return indArr[idim][iprtcl];
  }

  inline std::vector<int> id(size_t idim) 
  {
    //return indArr[idim];
    std::vector<int> ret;
    for(const auto& e: indArr[idim])
      ret.push_back(e);
    return ret;
  }

/*
  virtual inline std::vector<int>& id(size_t idim)
  {
    return indArr[idim];
  }
*/

  //--------------------------------------------------
  // EM fields
  DEVCALLABLE inline float_p& ex(size_t iprtcl ) { return Epart[0*size() + iprtcl]; };
  DEVCALLABLE inline float_p& ey(size_t iprtcl ) { return Epart[1*size() + iprtcl]; };
  DEVCALLABLE inline float_p& ez(size_t iprtcl ) { return Epart[2*size() + iprtcl]; };

  DEVCALLABLE inline float_p& bx(size_t iprtcl ) { return Bpart[0*size() + iprtcl]; };
  DEVCALLABLE inline float_p& by(size_t iprtcl ) { return Bpart[1*size() + iprtcl]; };
  DEVCALLABLE inline float_p& bz(size_t iprtcl ) { return Bpart[2*size() + iprtcl]; };

  DEVCALLABLE inline float_p ex(size_t iprtcl ) const {return Epart[0*size() + iprtcl]; };
  DEVCALLABLE inline float_p ey(size_t iprtcl ) const {return Epart[1*size() + iprtcl]; };
  DEVCALLABLE inline float_p ez(size_t iprtcl ) const {return Epart[2*size() + iprtcl]; };

  DEVCALLABLE inline float_p bx(size_t iprtcl ) const {return Bpart[0*size() + iprtcl]; };
  DEVCALLABLE inline float_p by(size_t iprtcl ) const {return Bpart[1*size() + iprtcl]; };
  DEVCALLABLE inline float_p bz(size_t iprtcl ) const {return Bpart[2*size() + iprtcl]; };


  // particle creation
  virtual void add_particle (
      std::vector<float_p> prtcl_loc,
      std::vector<float_p> prtcl_vel,
      float_p prtcl_wgt);

  // particle creation v2
  //virtual void add_particle2 (
  //  float_p lx, float_p ly, float_p lz, 
  //  float_p ux, float_p uy, float_p uz,
  //  float_p prtcl_wgt);

  // particle creation
  virtual void add_identified_particle (
      std::vector<float_p> prtcl_loc,
      std::vector<float_p> prtcl_vel,
      float_p prtcl_wgt, 
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
  
  // return energy of i:th particle 
  float_p get_prtcl_ene(size_t n);
  
  // sort particles in reverse (ascending) energy order
  void apply_permutation(ManVec<size_t>& indices);

  // sort particles in reverse (ascending) energy order
  void sort_in_rev_energy();

  // update internal cumulative weight arrays of particles
  void update_cumulative_arrays();
};




} // end of namespace pic
