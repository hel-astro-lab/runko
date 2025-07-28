#pragma once

#include <random>
#include <cassert>

#include "core/pic/tile.h"
#include "definitions.h"
#include "tools/vector.h"


namespace pic {

/// Magnetospheric gap
template<size_t D>
class Gap
{

private:

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float> uni_dis;

  // using raw pointer instead of smart ptrs; it does not take ownership of the object
  // so the container is not deleted when the temporary storage goes out of scope.
  using ConPtr = pic::ParticleContainer<D>* ;


public:

  // members
  double B0;     // B_bg
  double E0;     // E_rot
                      
  double gap_length; // length of the gap in dx
                       
  double x_left = 5;     // start of the box
  double x_right;        // end of the box
  double delta_left  = 2;  // smoothing parameter for left BC in units of dx
  double delta_right = 4;  // smoothing parameter for right BC in units of dx
  double Nx;     // total grid length

  int halo = 3; // extent of the halo regions not updated in left/right edges of the domain



  // constructor
  Gap() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)
  {
    assert(D==1); // only 1d implemented
  }

  // random numbers between [0, 1[
  float rand() { return uni_dis(gen); };

  // tile functions
  void insert_em(pic::Tile<D>&  tile);

  void update_b(pic::Tile<D>&  tile);

  void update_e(pic::Tile<D>&  tile);

  void update_j(pic::Tile<D>&  tile);

  void solve(pic::Tile<D>&  tile);
  
}; // end of class

} // end of namespace pic
