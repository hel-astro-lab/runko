#pragma once

#include <vector>
#include <map>
#include <cmath>

#include "../definitions.h"



namespace pic {

/*! \brief Block of particles inside the tile
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
class ParticleBlock {


  size_t Nprtcls;

  std::vector< std::vector<Realf> > locArr;
  std::vector< std::vector<Realf> > velArr;
  //std::vector< std::vector<Realf> > wgtArr;

  public:

  //! particle specific electric field components
  //std::vector<Realf> Epart;
  std::vector< std::vector<Realf> > Epart;

  //! particle specific magnetic field components
  //std::vector<Realf> Bpart;
  std::vector< std::vector<Realf> > Bpart;

  //! multimap of particles going to other tiles
  typedef std::multimap<std::tuple<int,int,int>, int> mapType;
  mapType to_other_tiles;


  // size of the internal mesh
  size_t Nx;
  size_t Ny;
  size_t Nz;

  Realf q = 1.0; // normalization factor


  /// Constructor 
  ParticleBlock(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz)
  { 
    locArr.resize(3);
    velArr.resize(3);
  
  };
    




  //--------------------------------------------------
    
  /// reserve memory for particles
  void reserve(uint64_t N, uint64_t D) {

    locArr.resize(3);
    for(uint64_t i=0; i<D; i++) locArr[i].reserve(N);

    velArr.resize(D);
    for(uint64_t i=0; i<D; i++) velArr[i].reserve(N);

    // reserve 1d N x D array for particle-specific fields
    // XXX is this needed since we can not push_back?
    //Epart.reserve(N*D);
    //Bpart.reserve(N*D);
  }

  // resize everything
  void resize(uint64_t N) {

    for(uint64_t i=0; i<3; i++) locArr[i].resize(N);
    for(uint64_t i=0; i<3; i++) velArr[i].resize(N);

    Nprtcls = N;
  }

  // resize arrays for fields
  void resizeEM(uint64_t N, uint64_t D) {
    if (N <= 0) N = 1;

    Epart.resize(D);
    for(uint64_t i=0; i<D; i++) Epart[i].resize(N);

    Bpart.resize(D);
    for(uint64_t i=0; i<D; i++) Bpart[i].resize(N);

  }



  /// size of the container (in terms of particles)
  size_t size() { return locArr[0].size(); }



  //--------------------------------------------------
  // locations
  inline Realf loc( uint64_t idim, uint64_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  inline Realf& loc( uint64_t idim, uint64_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  inline std::vector<Realf> loc(uint64_t idim) const 
  {
    return locArr[idim];
  }


  //--------------------------------------------------
  // velocities
  inline Realf vel( uint64_t idim, uint64_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  inline Realf& vel( uint64_t idim, uint64_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  inline std::vector<Realf> vel(uint64_t idim) const 
  {
    return velArr[idim];
  }


  // --------------------------------------------------
  // particle creation & destruction methods

  void add_particle (
      std::vector<Realf> prtcl_loc,
      std::vector<Realf> prtcl_vel) 
  {

    for (size_t i=0; i<3; i++)  
      locArr[i].push_back(prtcl_loc[i]);

    for (size_t i=0; i<3; i++)  
      velArr[i].push_back(prtcl_vel[i]);

    Nprtcls++;
  }



  // --------------------------------------------------
  //! Lorentz factor
  inline Realf gamma(unsigned int iprtcl) {
    return sqrt(1.+
         pow(vel(0,iprtcl),2)
        +pow(vel(1,iprtcl),2)
        +pow(vel(2,iprtcl),2));
  }

  //! Reciprocal of Lorentz factor
  //inline Realf invgamma(unsigned int iprtcl) {
  //    return 1./sqrt(1.+pow(vel(0,iprtcl),2)
  //                     +pow(vel(1,iprtcl),2)
  //                     +pow(vel(2,iprtcl),2));
  //}





  // --------------------------------------------------
  // Other useful mistileaneous methods
    
  //! test if prtcl is inside this block
  //bool is_local(uint64_t iprtcl);

};





} // end of namespace pic
