#pragma once

#include <vector>
#include <map>
#include <cmath>




namespace pic {

/*! \brief Block of particles inside the cell
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

  std::vector< std::vector<double> > locArr;
  std::vector< std::vector<double> > velArr;
  //std::vector< std::vector<double> > wgtArr;

  public:

  //! particle specific electric field components
  //std::vector<double> Epart;
  std::vector< std::vector<double> > Epart;

  //! particle specific magnetic field components
  //std::vector<double> Bpart;
  std::vector< std::vector<double> > Bpart;

  //! multimap of particles going to other tiles
  typedef std::multimap<std::tuple<int,int,int>, int> mapType;
  mapType to_other_tiles;


  // size of the internal mesh
  size_t Nx;
  size_t Ny;
  size_t Nz;

  double q = 1.0; // normalization factor


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

    for(uint64_t i=0; i<3; i++) Epart[i].resize(N);
    for(uint64_t i=0; i<3; i++) Bpart[i].resize(N);

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
  inline double loc( uint64_t idim, uint64_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  inline double& loc( uint64_t idim, uint64_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  inline std::vector<double> loc(uint64_t idim) const 
  {
    return locArr[idim];
  }


  //--------------------------------------------------
  // velocities
  inline double vel( uint64_t idim, uint64_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  inline double& vel( uint64_t idim, uint64_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  inline std::vector<double> vel(uint64_t idim) const 
  {
    return velArr[idim];
  }


  // --------------------------------------------------
  // particle creation & destruction methods

  void add_particle (
      std::vector<double>& prtcl_loc,
      std::vector<double>& prtcl_vel) 
  {

    for (size_t i=0; i<3; i++)  
      locArr[i].push_back(prtcl_loc[i]);

    for (size_t i=0; i<3; i++)  
      velArr[i].push_back(prtcl_vel[i]);

  }



  // --------------------------------------------------
  //! Lorentz factor
  inline double gamma(unsigned int iprtcl) {
    return sqrt(1.+
         pow(vel(0,iprtcl),2)
        +pow(vel(1,iprtcl),2)
        +pow(vel(2,iprtcl),2));
  }

  //! Reciprocal of Lorentz factor
  //inline double invgamma(unsigned int iprtcl) {
  //    return 1./sqrt(1.+pow(vel(0,iprtcl),2)
  //                     +pow(vel(1,iprtcl),2)
  //                     +pow(vel(2,iprtcl),2));
  //}





  // --------------------------------------------------
  // Other useful miscellaneous methods
    
  //! test if prtcl is inside this block
  //bool is_local(uint64_t iprtcl);

};





} // end of namespace pic
