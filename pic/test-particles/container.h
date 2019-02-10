#pragma once

#include "../particle.h"


/*! \brief Test particle container
 *
 * Holds test particles that are a derivative
 * of normal particles but with additional properties.
 * Especially the I/O is different as we typically save
 * the full history.
 */
class TestParticleContainer :
  public virtual ParticleContainer
{

  /// Constructor 
  TestParticleContainer();

  // default virtual dtor
  virtual ~TestParticleContainer() = default;





  //--------------------------------------------------
  // id
  inline int id( size_t iprtcl ) const
  {
    return wgtArr[iprtcl];
  }

  inline int& id( size_t iprtcl )       
  {
    return wgtArr[iprtcl];
  }

  inline std::vector<int> id() const
  {
    return wgtArr;
  }





};








