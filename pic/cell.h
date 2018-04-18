#pragma once

#include "../definitions.h"
#include "../corgi/cell.h"

#include "../em-fields/fields.h"

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
  //std::vector< std::vector<double> > Weight;

  //! particle specific electric field components
  std::vector<double> Epart;

  //! particle specific magnetic field components
  std::vector<double> Bpart;


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
  //! Lorentz factor
    inline double gamma(unsigned int iprtcl) {
        return sqrt(1.+pow(vel(0,iprtcl),2)
                      +pow(vel(1,iprtcl),2)
                      +pow(vel(2,iprtcl),2));
    }

    //! Reciprocal of Lorentz factor
    //inline double invgamma(unsigned int iprtcl) {
    //    return 1./sqrt(1.+pow(vel(0,iprtcl),2)
    //                     +pow(vel(1,iprtcl),2)
    //                     +pow(vel(2,iprtcl),2));
    //}




  //! test if prtcl is inside this block
  //bool is_local(uint64_t iprtcl);

};




/*! \brief PiC cell
 *
 * Cell infrastructures are inherited from corgi::Cell
 * Maxwell field equation solver is inherited from fields::PlasmaCell
*/
class PicCell :
  virtual public fields::PlasmaCell, 
  virtual public corgi::Cell {

public:

  /// Size of the internal grid
  size_t NxGrid;
  size_t NyGrid;
  size_t NzGrid;

  size_t Nspecies = 2;
  size_t Nsteps   = 2;


  ParticleBlock container;


  /// constructor
  PicCell(size_t i, size_t j, 
             int o, 
             size_t NxG, size_t NyG,
             size_t NxMesh, size_t NyMesh
             ) : 
    corgi::Cell(i, j, o, NxG, NyG),
    fields::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh, 1),
    NxGrid(NxMesh),
    NyGrid(NyMesh),
    NzGrid(1)
  { }


  /// destructor
  ~PicCell() { };

  /// cell temporal and spatial scales
  Realf dt = 0.0;
  Realf dx = 0.0;



};

} // end of namespace pic
