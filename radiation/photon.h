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
 *   Based on pic::Container<3> i.e., assuming 3D geometry
 *   by default. This is to simplify the algorithm.
 */
class PhotonContainer : 
  public pic::ParticleContainer<3>
{
  public:

  /// Constructor 
  PhotonContainer()  
    
  { };

  //virtual ~PhotonContainer() = default;

  /// particle energy (h\nu in m_e c^2)
  std::vector< float_p > eneArr;


  /// initializes internal arrays
  void reserve(size_t N) override
  {
    eneArr.reserve(N);
    pic::ParticleContainer<3>::reserve(N);
  }


  /// resize 
  void resize(size_t N) override
  {
    eneArr.resize(N);
    pic::ParticleContainer<3>::resize(N);
  }

  /// special method 
  void add_particle(
    std::vector<float_p> prtcl_loc,
    std::vector<float_p> prtcl_vel,
    float_p weight,
    float_p energy) 
  {
    eneArr.push_back(energy);
    pic::ParticleContainer<3>::add_particle(prtcl_loc, prtcl_vel, weight);
  }

  // explicitly disallow the usage of base class member
  private:
    using pic::ParticleContainer<3>::add_particle;

};









} // end of namespace radiation
