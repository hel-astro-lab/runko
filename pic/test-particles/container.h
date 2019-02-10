#pragma once

#include "../particle.h"

namespace pic {

/*! \brief Test particle container
 *
 * Holds test particles that are a derivative
 * of normal particles but with additional properties.
 * Especially the I/O is different as we typically save
 * the full history.
 */
class TestParticleContainer :
  public pic::ParticleContainer
{

  private:

    /// global mpi rank
    int rank;

    /// comm size
    int mpi_world_size;

    /// running key generator seed
    int key = 1;

    /// unique key generator
    double keygen();


  public:

  /// Constructor 
  TestParticleContainer(); 

  // default virtual dtor
  virtual ~TestParticleContainer() = default;

  /// add particle with automatic unique key generation
  void add_test_particle (
      std::vector<double> prtcl_loc,
      std::vector<double> prtcl_vel);


  //--------------------------------------------------
  // id
  //inline int id( size_t iprtcl ) const
  //{
  //  return wgtArr[iprtcl];
  //}

  //inline int& id( size_t iprtcl )       
  //{
  //  return wgtArr[iprtcl];
  //}

  //inline std::vector<int> id() const
  //{
  //  return wgtArr;
  //}

  // overwrite weight and return unity
  // This is to make sure that test particles behave
  // correctly with depositer etc. that need weights.

  // note: can't do this because ParticleContainer::pack_xxx
  // etc use wgt(ind)



};



} // end of ns pic
