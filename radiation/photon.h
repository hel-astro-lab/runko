#pragma once

#include "../pic/particle.h"



namespace rad {

/*! \brief Block of photons inside the tile
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
class PhotonBlock : 
  virtual public pic::ParticleBlock 
{
  public:

  /// Constructor 
  PhotonBlock(size_t Nx, size_t Ny, size_t Nz) : 
    pic::ParticleBlock(Nx, Ny, Nz)
  { };

  virtual ~PhotonBlock() = default;


  /// particle weight
  std::vector< Realf > wgtArr;

  /// particle energy (h\nu in m_e c^2)
  std::vector< Realf > eneArr;



  /// initializes internal arrays
  virtual void reserve(size_t N) override
  {
    wgtArr.reserve(N);
    eneArr.reserve(N);

    pic::ParticleBlock::reserve(N);
  }


  /// resize 
  void resize(size_t N) override
  {
    // reserve size for photon specific stuff
    wgtArr.resize(N);
    eneArr.resize(N);

    // and finally to original arrays as well
    pic::ParticleBlock::resize(N);
  }

  /// special method 
  void add_particle(
    std::vector<Realf> prtcl_loc,
    std::vector<Realf> prtcl_vel,
    Realf weight,
    Realf energy) 
  {
    wgtArr.push_back(weight);
    eneArr.push_back(energy);

    pic::ParticleBlock::add_particle(prtcl_loc, prtcl_vel);
  }


  // explicitly disallow the usage of base class member
  private:
    using pic::ParticleBlock::add_particle;
    using pic::ParticleBlock::resize_em;

};









} // end of namespace radiation
