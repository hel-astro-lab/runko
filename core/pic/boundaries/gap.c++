
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

  // loop indices
  const int imin = -3, imax = tile.mesh_lengths[0]+3;

  #pragma omp simd
  for(int i=imin; i<imax; i++) {

    // global grid coordinates
    float iglob = (D>=1) ? i + mins[0] : 0;
    float h = iglob - x_left; // height in units of cells, h=0 is the surface
      
    //--------------------------------------------------
    // magnetic field

    gs.bx(i,0,0) = B0; // constant background field
    gs.by(i,0,0) = 0.0;
    gs.bz(i,0,0) = 0.0;

    //--------------------------------------------------
    // electric field

    // linear profile 
    float erot = h < gap_length ? E0*( 1.0 - h/gap_length ) : 0.0;

    const float erot1 = erot;
    const float erot2 = 0.0f;
    const float erot3 = 0.0f;

    //--------------------------------------------------
    // blending of old + new solution
    auto s  = 1.0f - shape( h, 0.0, delta_left); // height smoothing parameter

    gs.ex(i,0,0) = s*erot1 + (1.0f - s)*gs.ex(i,0,0); 
    gs.ey(i,0,0) = s*erot2 + (1.0f - s)*gs.ey(i,0,0); 
    gs.ez(i,0,0) = s*erot3 + (1.0f - s)*gs.ez(i,0,0); 
  }

  return;
}


//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
template<size_t D>
void pic::Gap<D>::update_b(
    pic::Tile<D>& tile)
{
    
  // Tile limits
  auto mins     = tile.mins;
  auto maxs     = tile.maxs;
  auto& gs      = tile.get_grids();

  // loop indices
  const int imin = -halo, imax = tile.mesh_lengths[0]+halo;

  // find top and bottom tiles and only operate on them
  bool top   = false;
  bool bot   = false;
  if( mins[0] < x_left  +3*delta_left  + tile.mesh_lengths[0]) bot   = true; 
  if( maxs[0] > x_right -3*delta_right - tile.mesh_lengths[0]) top   = true; 


  // bottom boundary
  if( bot ) {
    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      float h = iglob - x_left; // height in units of cells, h=0 is the surface
        
      //--------------------------------------------------
      // magnetic field

      float bx = B0; // constant background field
      float by = 0.0;
      float bz = 0.0;

      //--------------------------------------------------
      // blending of old + new solution
      auto s  = shape( h, 0.0, delta_left); // height smoothing parameter

      gs.bx(i,0,0) = s*bx + (1.0f - s)*gs.bx(i,0,0); 
      gs.by(i,0,0) = s*by + (1.0f - s)*gs.by(i,0,0); 
      gs.bz(i,0,0) = s*bz + (1.0f - s)*gs.bz(i,0,0); 
    }
  }

  //-------------------------------------------------- 
  // top boundary
  if( top ) {

    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
        
      //--------------------------------------------------
      // magnetic field

      float bx = B0; // constant background field
      float by = 0.0;
      float bz = 0.0;

      //--------------------------------------------------
      // blending of old + new solution
      auto s  = 1.0f - shape( iglob, x_right, delta_right); // height smoothing parameter

      gs.bx(i,0,0) = s*bx + (1.0f - s)*gs.bx(i,0,0); 
      gs.by(i,0,0) = s*by + (1.0f - s)*gs.by(i,0,0); 
      gs.bz(i,0,0) = s*bz + (1.0f - s)*gs.bz(i,0,0); 
    }
  }


  //--------------------------------------------------
  // hard-coded left/star BC
  if( mins[0] < 1 ) {
    for(int i=imin; i<=halo; i++) {
      gs.bx(i,0,0) = B0; 
      gs.by(i,0,0) = 0.0; 
      gs.bz(i,0,0) = 0.0; 
    }
  }

  //--------------------------------------------------
  // hard-coded right/const BC
  if( maxs[0] > Nx-1 ) {
    for(int i=imax-2*halo; i<=imax; i++) {
      gs.bx(i,0,0) = B0; 
      gs.by(i,0,0) = 0.0; 
      gs.bz(i,0,0) = 0.0; 
    }
  }

  return;
}




//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
template<size_t D>
void pic::Gap<D>::update_e(
    pic::Tile<D>& tile)
{
  // Tile limits
  auto mins     = tile.mins;
  auto maxs     = tile.maxs;
  auto& gs      = tile.get_grids();

  // loop indices
  const int imin = -halo, imax = tile.mesh_lengths[0]+halo;

  // find top and bottom tiles and only operate on them
  bool top   = false;
  bool bot   = false;
  if( mins[0] < x_left  +3*delta_left  + tile.mesh_lengths[0]) bot = true; 
  if( maxs[0] > x_right -3*delta_right - tile.mesh_lengths[0]) top = true; 


  // bottom boundary
  //if( bot ) {
  //  #pragma omp simd
  //  for(int i=imin; i<imax; i++) {

  //    // global grid coordinates
  //    float iglob = (D>=1) ? i + mins[0] : 0;
  //    float h = iglob - x_left; // height in units of cells, h=0 is the surface

  //    // linear profile 
  //    float erot = h < gap_length ? E0*( 1.0 - h/gap_length ) : 0.0;

  //    const float erot1 = erot;
  //    const float erot2 = 0.0f;
  //    const float erot3 = 0.0f;

  //    //--------------------------------------------------
  //    // blending of old + new solution
  //    auto s = 1.0f - shape( h, 0.0, delta_left); // height smoothing parameter

  //    gs.ex(i,0,0) = s*erot1 + (1.0f - s)*gs.ex(i,0,0); 
  //    gs.ey(i,0,0) = s*erot2 + (1.0f - s)*gs.ey(i,0,0); 
  //    gs.ez(i,0,0) = s*erot3 + (1.0f - s)*gs.ez(i,0,0); 
  //  }
  //}

  //-------------------------------------------------- 
  // top boundary
  if( top ) {
    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      auto s  = 1.0f - shape( iglob, x_right, delta_right); // height smoothing parameter
                                                              
      gs.ex(i,0,0) = s*(0.0f) + (1.0f - s)*gs.ex(i,0,0); 
      gs.ey(i,0,0) = s*(0.0f) + (1.0f - s)*gs.ey(i,0,0); 
      gs.ez(i,0,0) = s*(0.0f) + (1.0f - s)*gs.ez(i,0,0); 
    }
  }


  //--------------------------------------------------
  // hard-coded left (star) BC
  if( mins[0] < 1 ) {
    for(int i=imin; i<=halo; i++) {
      gs.ex(i,0,0) = 0.0;
      gs.ey(i,0,0) = 0.0; 
      gs.ez(i,0,0) = 0.0; 
    }
  }

  //--------------------------------------------------
  // hard-coded right/vacuum BC
  if( maxs[0] > Nx-1 ) {
    for(int i=imax-2*halo; i<=imax; i++) {
      gs.ex(i,0,0) = 0.0; 
      gs.ey(i,0,0) = 0.0; 
      gs.ez(i,0,0) = 0.0; 
    }
  }

  return;
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
