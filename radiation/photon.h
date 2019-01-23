#pragma once

#include "../pic/particle.h"



namespace rad {

/*! \brief Container of photons inside the tile
 *
 * Container to hold photon particles, i.e., computational
 * particles with v=c.
 *
 * includes:
 *   pos/loc
 *   vel
 *   wgt
 *   qm
 *
 */
class PhotonContainer : 
  virtual public pic::ParticleContainer
{
  public:

  /// Constructor 
  PhotonContainer() : 
    pic::ParticleContainer()
  { };

  virtual ~PhotonContainer() = default;

  /// particle energy (h\nu in m_e c^2)
  std::vector< double > eneArr;


  /// initializes internal arrays
  virtual void reserve(size_t N) override
  {
    eneArr.reserve(N);
    pic::ParticleContainer::reserve(N);
  }


  /// resize 
  void resize(size_t N) override
  {
    eneArr.resize(N);
    pic::ParticleContainer::resize(N);
  }

  /// special method 
  void add_particle(
    std::vector<double> prtcl_loc,
    std::vector<double> prtcl_vel,
    double weight,
    double energy) 
  {
    eneArr.push_back(energy);
    pic::ParticleContainer::add_particle(prtcl_loc, prtcl_vel, weight);
  }

  // explicitly disallow the usage of base class member
  private:
    using pic::ParticleContainer::add_particle;

};









} // end of namespace radiation
