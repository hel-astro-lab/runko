#pragma once

#include "fields.h"



namespace fields {


/* \brief Damping electro-magnetic cell
 *
 */
class PlasmaCellDamped : 
  virtual public fields::PlasmaCell 
{

  public:

    PlasmaCellDamped(
      size_t i, size_t j, 
      int o,
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh, size_t NzMesh);

    ~PlasmaCellDamped() { };


    void pushE();
    void pushE2d_damped();


};












} // end on namespace
