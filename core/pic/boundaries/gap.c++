
#include "core/pic/boundaries/gap.h"
#include "tools/vector.h"
#include "tools/signum.h"
#include "tools/staggered_grid.h"

#include <cmath> 
#include <cassert>
#include <string>
#include <map>


using std::min;
using std::max;
using std::abs;
using std::sqrt;
using std::sin;
using std::cos;

using toolbox::Vec3;
using toolbox::norm;
using toolbox::norm1d;
using toolbox::norm2d;
using toolbox::StaggeredSphericalCoordinates;
using toolbox::shape; // tanh function



template<size_t D>
void pic::Gap<D>::insert_em(
    pic::Tile<D>& tile)
{

  // Tile limits
  auto mins     = tile.mins;
  auto& gs      = tile.get_grids();
  const float c = tile.cfl; // (numerical) speed of light


  // loop indices
  const int imin = -3, imax = tile.mesh_lengths[0]+3;
  const int jmin =  0, jmax = 1;
  const int kmin =  0, kmax = 1;


  #pragma omp simd
  for(int i=imin; i<imax; i++) {

    // global grid coordinates
    float iglob = (D>=1) ? i + mins[0] : 0;
    float h = iglob - radius; // height in units of cells, h=0 is the surface
      
    //--------------------------------------------------
    // magnetic field

    gs.bx(i,0,0) = B0; // constant background field
    gs.by(i,0,0) = 0.0;
    gs.bz(i,0,0) = 0.0;


    //--------------------------------------------------
    // electric field

    // linear profile 
    float erot = h < radius_pc ? E0*( 1.0 - h/radius_pc ) : 0.0;

    const float erot1 = erot;
    const float erot2 = 0.0f;
    const float erot3 = 0.0f;

    //--------------------------------------------------
    // blending of old + new solution

    //auto s  = 1.0f;
    auto s  = 1.0f - shape( h, 0.0, delta); // height smoothing parameter

    gs.ex(i,0,0) = s*erot1 + (1.0f - s)*gs.ex(i,0,0); 
    gs.ey(i,0,0) = s*erot2 + (1.0f - s)*gs.ey(i,0,0); 
    gs.ez(i,0,0) = s*erot3 + (1.0f - s)*gs.ez(i,0,0); 
  }

}



template<size_t D>
void pic::Gap<D>::update_b(
    pic::Tile<D>& tile)
{

}



template<size_t D>
void pic::Gap<D>::update_e(
    pic::Tile<D>& tile)
{

}


// current from moving frame
template<size_t D>
void pic::Gap<D>::update_j(
  pic::Tile<D>& tile)
{

}


template<size_t D>
void pic::Gap<D>::solve(
    pic::Tile<D>& tile)
{


}






//--------------------------------------------------
// explicit template instantiation
template class pic::Gap<1>; // 1D only
