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
    ParticleBlock(Nx, Ny, Nz)
  { };

  /// particle weight
  std::vector< double > wgtArr;

  /// particle energy (h\nu in m_e c^2)
  std::vector< double > eneArr;



  /// initializes internal arrays
  virtual void reserve(size_t N) override
  {
    wgtArr.reserve(N);
    eneArr.reserve(N);

    locArr.resize(3);
    for(size_t i=0; i<3; i++) locArr[i].reserve(N);

    velArr.resize(3);
    for(size_t i=0; i<3; i++) velArr[i].reserve(N);
  }


  /// resize 
  void resize(size_t N) override
  {
    // reserve size for photon specific stuff
    wgtArr.resize(N);
    eneArr.resize(N);

    // and finally to original arrays as well
    ParticleBlock::resize(N);

    //for(size_t i=0; i<3; i++) locArr[i].resize(N);
    //for(size_t i=0; i<3; i++) velArr[i].resize(N);
    //Nprtcls = N;
  }

  /// special method 
  void add_particle(
    std::vector<Realf> prtcl_loc,
    std::vector<Realf> prtcl_ang,
    Realf weight,
    Realf energy) 
  {
    assert(prtcl_loc.size() == 3);
    assert(prtcl_ang.size() == 3);

    for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
    for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_ang[i]);

    wgtArr.push_back(weight);
    eneArr.push_back(energy);

    Nprtcls++;
  }

  //inline std::vector<Realf> wgt(size_t idim) const 
  //{
  //  return wgtArr[idim];
  //}


  // explicitly disallow the usage of base class member
  private:
    using pic::ParticleBlock::add_particle;


};









} // end of namespace radiation
