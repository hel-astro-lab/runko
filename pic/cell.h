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



  //ParticleBlock container;
  std::vector<ParticleBlock> containers;

  /// get i:th container
  ParticleBlock& get_container(size_t i) { return containers[i]; };

  /// set i:th container
  void set_container(ParticleBlock block) {containers.push_back(block);};

  size_t Nspecies() { return containers.size(); };


  /// constructor
  PicCell(size_t i, size_t j, 
             int o, 
             size_t NxG, size_t NyG,
             size_t NxMesh, size_t NyMesh
             ) : 
    corgi::Cell(i, j, o, NxG, NyG),
    fields::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh, 1),
    NxGrid(NxMesh), NyGrid(NyMesh), NzGrid(1)
  { }


  /// destructor
  ~PicCell() override { };

  /// cell temporal and spatial scales
  using fields::PlasmaCell::cfl;
  using fields::PlasmaCell::dt;
  using fields::PlasmaCell::dx;


};



} // end of namespace pic
