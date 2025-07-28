
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

  // Tile limits
  auto mins = tile.mins;
  auto& gs  = tile.get_grids();
  const float c = tile.cfl; // Delta t
  const float v = E0/B0; // v_rot
 
  // loop indices
  const int imin = 0, imax = tile.mesh_lengths[0]; // NOTE: no halos

  // set current
  #pragma omp simd
  for(int i=imin; i<imax; i++) {

    // global grid coordinates
    float ig = (D>=1) ? i + mins[0] : 0;
    float jx=0.0f, jy=0.0f, jz=0.0f;

    //--------------------------------------------------
    // current in 1D moving frame v = v \hat{y}
      
    // j_m = 
    //x:  - v d( E_y)
    //y:  -d( v E_x) + d( v^2 B_z)
    //z:  0

    // jx
    float dx_ey_at_x = gs.ey(i,0,0) - gs.ey(i-1,0,0); // term2: partial_x(E_y) at i,j,k
    //float dx_ey_at_x = gs.ey(i+1,j,k) - gs.ey(i,j,k); // term2: partial_x(E_y) at i,j,k // BAD: oscillates
    //jx = -v*dx_ey_at_x; // IGNORED term2

    // jy
    float dx_v_ex_at_y  = v*gs.ex(i,0,0) - v*gs.ex(i-1,0,0); // term1: partial_y(v E_x) at i,j,k
    //float dx_v_ex_at_y  = v*gs.ex(i+1,j,k) - v*gs.ex(i,j,k); // term1: partial_y(v E_x) at i,j,k // BAD: oscillates
	  
    float dx_v2_bz_at_y = v*v*gs.bz(i,0,0) - v*v*gs.bz(i-1,0,0);   // term3: partial_x( v^2 B_z )
    //float dx_v2_bz_at_y = v*v*gs.bz(i+1,j,k) - v*v*gs.bz(i,j,k); // term3: partial_x (v^2 B_z) // UNTESTED
    //jy = -dx_v_ex_at_y + dx_v2_bz_at_y; // IGNORED term3
    jy = -dx_v_ex_at_y; // ONLY term1 active

    //jz
    jz = 0.0f;

    //--------------------------------------------------
    // ver 1; smooth tanh profile
    auto s_l = 1.0f - shape( ig, x_left,  delta_left); 
    auto s_r =        shape( ig, x_right, delta_right); 
    auto s = s_l*s_r;

    //--------------------------------------------------
    // add to the current
    gs.jx(i,0,0) += s*jx*c;
    gs.jy(i,0,0) += s*jy*c;
    gs.jz(i,0,0) += s*jz*c;
  }


}


template<size_t D>
void pic::Gap<D>::solve(
    pic::Tile<D>& tile)
{


}






//--------------------------------------------------
// explicit template instantiation
template class pic::Gap<1>; // 1D only
