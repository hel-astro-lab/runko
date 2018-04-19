#pragma once

#include "../definitions.h"
#include "../corgi/cell.h"

#include "particle.h"

#include "../em-fields/fields.h"

namespace pic {


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
    NxGrid(NxMesh), NyGrid(NyMesh), NzGrid(1),
    container(NxMesh, NyMesh, 1)
  { }


  /// destructor
  ~PicCell() { };

  /// cell temporal and spatial scales
  Realf dt = 0.0;
  Realf dx = 0.0;

};



} // end of namespace pic
