
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


//-------------------------------------------------- 
//  magnetic  field
template<size_t D>
float pic::Gap<D>::B(float h) 
{
  if(      b_profile_mode == 0 ) { return B0; }
  else if( b_profile_mode == 1 ) { return B0*std::max(0.0f, 1.0f - h/gap_length); } // linear slope
  else if( b_profile_mode == 2 ) { return B0*std::max(0.0f, std::min(1.0f, 1.0f - h/Nx)); } // long linear (full box)
  else if( b_profile_mode == 3 ) { return B0*shape(h + x_left, x_right, delta_right); } // const + smooth damping
  else if( b_profile_mode == 4 ) { return B0*std::max(0.0f, std::min(1.0f, 1.0f - (h - gap_length)/gap_length)); } // const + linear 

  assert(false);
}


//-------------------------------------------------- 
//  electric field
template<size_t D>
float pic::Gap<D>::E(float h) 
{
  if(h <= 0) return 0.0f; // set always to zero inside star to ensure proper interior boundaries

  //--------------------------------------------------
  if(      e_profile_mode == 0 ) { return std::max(0.0f, E0*std::min(1.0f, (1.0f - h/gap_length) ) ); } // Ruderman-style setups
  else if( e_profile_mode == 1 ) { return std::max(0.0f, E0*std::min(1.0f, (1.0f - h/gap_length) ) ) -E0; } // SCLF 
  else if( e_profile_mode == 2 ) { return 0.0f; } // zero everywhere
  else if( e_profile_mode == 3 ) { 
    // v1: linear gap
    float erot = h < gap_length ? -E0*( h/gap_length ) : 0.0; // value outside
    erot = h > gap_length ? -E0 : erot; // const value higher up
    return erot;
  }
  else if( e_profile_mode == 4 ) { 
    // v2: co-rotation charge density is taken to depend on dimensionless height h = (x/H) as:
    //   eta_co/eta_co0 = 1 + A*h
    //
    // then, the electric field is 
    //    E(x) = 4\pi \int (eta_co0 - eta_co) dx
    //    E(x) = -E_rot * 0.5*A*h^2
    //const float eta_co_A = 0.07f; //  realistic value
    const float eta_co_A = 0.8f; // numerical value
    const float e_max = 0.5*E0*eta_co_A;
    float erot = h < gap_length ? e_max*std::pow( h/gap_length, 2 ) : 0.0f;
    return erot;
  }

  assert(false);
}


// is bottom (left) tile?
template<size_t D>
bool pic::Gap<D>::is_bot(pic::Tile<D>& tile) 
{
  return tile.mins[0] < x_left + (3*delta_left + tile.mesh_lengths[0]);
}

// is top (right) tile?
template<size_t D>
bool pic::Gap<D>::is_top(pic::Tile<D>& tile) 
{
  return tile.maxs[0] > x_right - (3*delta_right + tile.mesh_lengths[0]);
}


//template<size_t D>
//float pic::Gap<D>::blend_right(float ig, float val, float tar) 
//{
//  // height smoothing parameter
//  auto s  = shape( ig, x_right, delta_right); 
//  // shape of s(x): ----\____ rightmost BC
//
//
//  auto dv = tar - val; // \delta F; how much value deviates from target
//  return tar - s*dv;
//  return tar - s*tar + s*val = tar*(1-s) + s*val
//
//  auto s  = 1.0f - shape( ig, x_right, delta_right); 
//  // shape of s(x): ___ / ---- rightmost BC
//  return s*tar + (1.0f - s)*val;
//}


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
    gs.bx(i,0,0) = B(h);
    gs.by(i,0,0) = 0.0;
    gs.bz(i,0,0) = 0.0;

    //--------------------------------------------------
    // electric field
    float erot = E(h); // E(x)
    float erot_inside = E(0.0f); // field inside star

    //--------------------------------------------------
    // blending of inside/outside solutions; removes sharp jumps between interior and exterior
    auto s  = shape( h, 0.0, delta_left); // height smoothing parameter

    gs.ex(i,0,0) = s*erot_inside + (1.0f-s)*erot; 
    gs.ey(i,0,0) = 0.0f; 
    gs.ez(i,0,0) = 0.0f; 
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
  const int imin = 0, imax = tile.mesh_lengths[0];
  //const int imin = -halo, imax = tile.mesh_lengths[0]+halo;

  // find top and bottom tiles and only operate on them
  bool bot = is_bot(tile);
  bool top = is_top(tile);

  // bottom boundary
  if( bot ) {
    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      float h = iglob - x_left; // height in units of cells, h=0 is the surface
        
      //--------------------------------------------------
      // magnetic field
      float bx = B(0.0f); // constant background field
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
  //if( true ) { // NOTE: some profiles need to be always on since they bleeds to left on every step
  //if( top ) {
  //  #pragma omp simd
  //  for(int i=imin; i<imax; i++) {

  //    // global grid coordinates
  //    float iglob = (D>=1) ? i + mins[0] : 0;
  //      
  //    //--------------------------------------------------
  //    // magnetic field
  //    float bx = B(Nx); 
  //    float by = 0.0;
  //    float bz = 0.0;

  //    //--------------------------------------------------
  //    // blending of old + new solution
  //    auto s  = 1.0f - shape( iglob, x_right, delta_right); // height smoothing parameter

  //    gs.bx(i,0,0) = s*bx + (1.0f - s)*gs.bx(i,0,0); 
  //    gs.by(i,0,0) = s*by + (1.0f - s)*gs.by(i,0,0); 
  //    gs.bz(i,0,0) = s*bz + (1.0f - s)*gs.bz(i,0,0); 
  //  }
  //}


  //--------------------------------------------------
  // hard-coded left/star BC
  if( mins[0] < 1 ) {
    #pragma omp simd
    for(int i=-halo; i<=halo; i++) {
      gs.bx(i,0,0) = B(0.0f); 
      gs.by(i,0,0) = 0.0f; 
      gs.bz(i,0,0) = 0.0f; 
    }
  }

  //--------------------------------------------------
  // hard-coded right/const BC
  //if( maxs[0] > Nx-1 ) {
  //  #pragma omp simd
  //  for(int i=tile.mesh_lengths[0]-halo; i<=tile.mesh_lengths[0]+halo; i++) {
  //    gs.bx(i,0,0) = B(0.0f); 
  //    gs.by(i,0,0) = 0.0; 
  //    gs.bz(i,0,0) = 0.0; 
  //  }
  //}

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
  //const int imin = 0, imax = tile.mesh_lengths[0];
  const int imin = -halo, imax = tile.mesh_lengths[0]+halo;

  // find top and bottom tiles and only operate on them
  bool bot = is_bot(tile);
  bool top = is_top(tile);


  // bottom boundary
  if( bot ) {
    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      float h = iglob - x_left; // height in units of cells, h=0 is the surface
                                  
      // blending of boundary condition + active solution

      // sharp cutoff
      //auto s = h < 0 ? 1.0f : 0.0f; // here h=0 is the first "real" cell with dynamic E field

      //// height smoothing parameter
      auto s = shape( h, 0.0, delta_left); 

      gs.ex(i,0,0) = s*E(0.0f) + (1.0f - s)*gs.ex(i,0,0); 
      gs.ey(i,0,0) = s*0.0f    + (1.0f - s)*gs.ey(i,0,0); 
      gs.ez(i,0,0) = s*0.0f    + (1.0f - s)*gs.ez(i,0,0); 
    }
  }

  //-------------------------------------------------- 
  // top boundary; smooth w/ tanh profile
  if( top ) {
    #pragma omp simd
    for(int i=imin; i<imax; i++) {
      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      auto s  = 1.0f - shape( iglob, x_right, delta_right); // height smoothing parameter
                                                              
      gs.ex(i,0,0) = s*E(0.0f)+ (1.0f - s)*gs.ex(i,0,0); 
      gs.ey(i,0,0) = s*(0.0f) + (1.0f - s)*gs.ey(i,0,0); 
      gs.ez(i,0,0) = s*(0.0f) + (1.0f - s)*gs.ez(i,0,0); 
    }
  }

  // copy values from last point
  //if( maxs[0] >= Nx ) {
  //  //--------------------------------------------------
  //  // copy the last "real" value to the rest of the grid
  //  int i_r = x_right - mins[0]; 
  //  assert(i_r >= 0); // check that we are on the right tile
  //  //float ex_bc = gs.ex(i_r,0,0); // value at the edge
  //    
  //  for(int i=i_r+1; i<=tile.mesh_lengths[0]+halo; i++) {
  //    gs.ex(i,0,0) = 0.0f; //ex_bc; 
  //  }
  //}


  //--------------------------------------------------
  // hard-coded left (star) BC
  if( mins[0] < 1 ) {
    #pragma omp simd
    for(int i=-halo; i<=halo; i++) {
      gs.ex(i,0,0) = E(0.0f);
      gs.ey(i,0,0) = 0.0f; 
      gs.ez(i,0,0) = 0.0f; 
    }
  }


  //--------------------------------------------------
  // hard-coded right/vacuum BC
  if( maxs[0] >= Nx ) {
    #pragma omp simd
    for(int i=tile.mesh_lengths[0]; i<=tile.mesh_lengths[0]+halo; i++) {
      gs.ex(i,0,0) = E(0.0f); 
      gs.ey(i,0,0) = 0.0f; 
      gs.ez(i,0,0) = 0.0f; 
    }
  }

  return;
}


// current from moving frame
template<size_t D>
void pic::Gap<D>::add_jrot(
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
    // ver 1; smooth tanh profile inside the magnetosphere
    //const auto h_l = ig + std::max(1.0, 2.0*delta_left); // shift by +2 delta
    //const auto h_r = ig - std::max(1.0, 2.0*delta_left); // shift by -2 delta

    //auto s_l = 1.0f - shape( h, x_left,  delta_left); 
    //auto s_r = 1.0f; //shape( h_r, x_right, delta_right); 
    //auto s = s_l*s_r;
      
    auto s = 1.0f; // no boundary smoothing

    //--------------------------------------------------
    // add to the current
    gs.jx(i,0,0) += s*jx*c;
    gs.jy(i,0,0) += s*jy*c;
    gs.jz(i,0,0) += s*jz*c;
  }

  return;
}


template<size_t D>
void pic::Gap<D>::add_jext(
  pic::Tile<D>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto& gs  = tile.get_grids();
  const float c = tile.cfl; // Delta t
 
  // loop indices
  const int imin = 0, imax = tile.mesh_lengths[0]; // NOTE: no halos
  //const int imin = -3, imax = tile.mesh_lengths[0]+3;

  // set current
  #pragma omp simd
  for(int i=imin; i<imax; i++) {

    // global grid coordinates
    float ig = (D>=1) ? i + mins[0] : 0;

    // external current
    float jx = j_ext*( B(ig - x_left)/B0 ); 
    //float jy = 0.0f;
    //float jz = 0.0f;

    //--------------------------------------------------
    // various kinds of suppressions at the grid edges
    //
    // suppress current at boundaries; double tanh profile
    //auto s_l = 1.0f - shape( ig, x_left,  delta_left); 
    //auto s_r =        shape( ig, x_right, delta_right); 

    // sharp cutoff at halos
    //float s_l = 0.0f ? ig < halo : 1.0f; // sharp cutoff
    //float s_r = 0.0f ? ig >= Nx - halo : 1.0f; // sharp cutoff
      
    //auto s_l = 1.0f - shape( ig, halo,    2); 
    //auto s_r =        shape( ig, Nx-halo, 2); 

    //--------------------------------------------------
    float s_l = 0.0f ? ig < x_left : 1.0f; // sharp cutoff
    float s_r = 1.0f; //? ig > x_right: 1.0f; // no cutoff
    auto s = s_l*s_r;
      
    //-------------------------------------------------- 
    //auto s = 1.0f; // NOTE: no smoothing

    gs.jx(i,0,0) -= jx*s; // external current
  }

  return;
}



template<size_t D>
void pic::Gap<D>::update_j(
  pic::Tile<D>& tile)
{
  // NOTE: this function is only called when filtering to set grid edge boundaries

  // tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;
  auto& gs  = tile.get_grids();
 
  // loop indices
  //const int imin = -3, imax = tile.mesh_lengths[0]+3;
  const int imin = 0, imax = tile.mesh_lengths[0];

  // find top and bottom tiles and only operate on them
  bool bot = is_bot(tile);
  bool top = is_top(tile);

  // operate only on roughly correct tiles 
  if(!(top || bot)) return;

  // set current
  //#pragma omp simd
  //for(int i=imin; i<imax; i++) {

  //  // global grid coordinates
  //  float ig = (D>=1) ? i + mins[0] : 0;

  //  // suppress current at boundaries; double tanh profile
  //  //auto s_l = 1.0f - shape( ig, x_left,  delta_left); 
  //  //auto s_r =        shape( ig, x_right, delta_right);  // tanh profile
  //    
  //  // sharp cutoff at halos
  //  float s_l = 0.0f ? ig <  halo : 1.0f; // sharp cutoff
  //  float s_r = 0.0f ? ig >= Nx   : 1.0f; // sharp cutoff

  //  // cutoff outside internal tile boundaries
  //  //float s_l2 = 0.0f ? i < 0   : 1.0f; // sharp cutoff
  //  //float s_r2 = 0.0f ? i >= Nx : 1.0f; // sharp cutoff
  //                                               
  //  // combine
  //  auto s  = s_l*s_r; // s now looks like 0 -> 1 -> 0
  //  //     s *= s_l2*s_r2; // s now looks like 0 -> 1 -> 0

  //  // suppression of current at the boundaries
  //  gs.jx(i,0,0) *= s;
  //  gs.jy(i,0,0) *= s;
  //  gs.jz(i,0,0) *= s;
  //}


  //--------------------------------------------------
  // hard-coded left (star) BC
  //if( bot ) {
  if( mins[0] < 1 ) {
    int i_l = x_left - 5; // location of the last point where current can be deposited
    assert(i_l >= 0); // check that we are on the right tile
    float jx_bc = gs.jx(i_l,0,0); // value at the edge

    //std::cout << " j_upd L: " 
    //  << " i-1: " << gs.jx(i_l-1,0,0)
    //  << " i  : " << gs.jx(i_l  ,0,0)
    //  << " i+1: " << gs.jx(i_l+1,0,0)
    //  << "\n";

    for(int i=-halo; i<i_l; i++) {
      gs.ex(i,0,0) = jx_bc;
      gs.ey(i,0,0) = 0.0f; 
      gs.ez(i,0,0) = 0.0f; 
    }
  }

  //--------------------------------------------------
  // hard-coded right/vacuum BC
  //if( top ) {
  if( maxs[0] >= Nx ) {
    int i_r = x_right - mins[0]; 
    assert(i_r >= 0); // check that we are on the right tile

    float jx_bc = gs.jx(i_r,0,0); // value at the edge
                            
    //std::cout << " j_upd R: " 
    //  << " i-1: " << gs.jx(i_r-1,0,0)
    //  << " i  : " << gs.jx(i_r  ,0,0)
    //  << " i+1: " << gs.jx(i_r+1,0,0)
    //  << "\n";

    for(int i=i_r+1; i<=tile.mesh_lengths[0]+halo; i++) {
      gs.jx(i,0,0) = jx_bc; 
      gs.jy(i,0,0) = 0.0f; 
      gs.jz(i,0,0) = 0.0f; 
    }
  }


  return;
}



template<size_t D>
void pic::Gap<D>::inject_prtcls(
    pic::Tile<D>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

                                                     
  // find top and bottom tiles and only operate on them
  bool bot = is_bot(tile);
  bool top = is_top(tile);

  // operate only on bottom tiles
  if(!bot) return;

  auto& gs = tile.get_grids(); // get fields grid

  //--------------------------------------------------
  // shortcut for containers 
  std::map<std::string, ConPtr> cons;
  for(auto&& con : tile.containers) cons.emplace(con.type, &con );

  // get charge (assume q_- = q_+)
  const float q = cons["e-"]->q;
  const float c = tile.cfl;
  const float n_co  = abs(E0/q)/gap_length; // co-rotation (Goldreich-Julian) density
  const float n_ext = abs(j_ext/q)/c/c; // maximum density to screen j_ext*dt


  //--------------------------------------------------
  // loop over grid points
  const int imin = 0, imax = tile.mesh_lengths[0]; // NOTE: no halos
                                                     
  for(int i=imin; i<imax; i++) {
    float iglob = (D>=1) ? i + mins[0] : 0;
    float h = iglob - x_left; // height in units of cells, h=0 is the surface

    bool inside_injection_layer = false;
      
    // v0: wide injection region
    //const float x_min_inj = x_left - std::max(1.0, 1.0*delta_left); // left boundary of atmosphere
    //const float x_max_inj = x_left + std::max(1.0, 1.0*delta_left); // right boundary of atmosphere
    //const float inj_width = x_max_inj - x_min_inj;  // width of the region in cells
    //if( (x_min_inj <= iglob) && (iglob <= x_max_inj) ) inside_injection_layer = true;

    // v1: narrow injection region, just behind the surface 
    if( h == -3.0f ) inside_injection_layer = true;
    const float inj_width = 1.0f;
      

    //std::cout << "inj: i" << iglob << " min:" << x_min_inj << " max:" << x_max_inj << " w:" << inj_width << " ins:" << inside_injection_layer << "\n";

    if( inside_injection_layer ) {

      // v1 (constant rate)
      //float ninj = inj_rate_pairs; // constant number of injections

      // v2 (screening of E_x)
      auto epar = gs.ex(i,0,0); // E_par
      //float ninj = inj_rate_pairs*abs(epar/E0)*n_co;

      // v3 (screening of j_ext)
      //float ninj = inj_rate_pairs*n_ext*abs(30.0f*epar/E0);
      float ninj = inj_rate_pairs*n_ext/inj_width;

      // add ninj pairs with MC injection; results on average in ninj injections
      float ncop = 0.0f; // number of pairs added
      float z1 = rand();

      //--------------------------------------------------
      //add particles
      while( ninj > z1 + ncop ) {

        // sample velocity from thermal distribution
        // using Box-Muller method to draw thermal velocities; valid for v <~ 0.2c
        float rr1 = rand(); //zeta4[ncop];
        float vr = sqrt( -2.0f*log(rr1))*temp_pairs;

        // Using now a random 3D distribution instead:
        float zeta      = 2.0f*PI*rand();
        float mu        = -1.0f + 2.0f*rand();
        float sin_theta = sqrt(1.0f-pow(mu,2));

        // construct velocities
        auto ux1 = vr*sin_theta*cos(zeta);
        auto uy1 = vr*sin_theta*sin(zeta);
        auto uz1 = vr*mu;

        //ux1 += 0.1f; // add upwards motion to help atmospheric current form

        // inject location is set randomly inside the cell
        // added to the first half of the domain between x_p \in [0, -0.5] 
        float dx = 0.5*rand();  
        //float dx = rand(); 
        //float dx = 0.0f;
                             
        cons["e-"]->add_particle( {{iglob + dx, 0.0f, 0.0f }}, {{ ux1, uy1, uz1 }}, wep);
        cons["p" ]->add_particle( {{iglob + dx, 0.0f, 0.0f }}, {{ ux1, uy1, uz1 }}, wep);

        ncop += 1;
      }
    } // end of inside_injection_layer



    //-------------------------------------------------- 
    // atmospheric photons
    
    if( inside_injection_layer ) {
      float ninj = inj_rate_phots; 

      float ncop = 0.0f; // number of phots added
      float z1 = rand();

      //--------------------------------------------------
      // add photons
      while( ninj > z1 + ncop ) {

        // draw random isotropic 3d vector
        float xia = rand();
        float xib = rand();
        float vx = 2.0f*xia -1.0f;
        float vy = 2.0f*sqrt(xia*(1.0f-xia))*cos(2.0f*PI*xib);
        float vz = 2.0f*sqrt(xia*(1.0f-xia))*sin(2.0f*PI*xib);

        //--------------------------------------------------
        // draw energy sample from a black body distribution
        float xi1 = rand();
        float xi2 = rand();
        float xi3 = rand();
        float xi4 = rand();
    
        float xi, jj, fsum;
        if( 1.202f*xi1 < 1.0f ){
            xi = 1.0f;
        } else {
            jj = 1.0f;
            fsum = std::pow(jj, -3);
            while( 1.202f*xi1 > fsum + std::pow(jj + 1.0f, -3) )
            {
                jj   += 1.0f;
                fsum += std::pow(jj, -3);
            }
            xi = jj + 1.0f;
        }
        float xinj = -temp_phots*std::log( xi2*xi3*xi4 )/xi;

        //--------------------------------------------------
        auto ux = xinj*vx;
        auto uy = xinj*vy;
        auto uz = xinj*vz;

        float dx = rand(); // inject location is set randomly inside the cell
                             
        cons["ph"]->add_particle( {{iglob + dx, 0.0f, 0.0f }}, {{ ux, uy, uz }}, wph);

        ncop += 1;
      }

    
    } // end of inside_injection_layer
  } // end of for loop over grid points


  return;
}


template<size_t D>
void pic::Gap<D>::delete_prtcls(
    pic::Tile<D>& tile)
{

  // more aggressive boundaries for particle removal
  bool bot = tile.mins[0] < x_left;
  bool top = tile.maxs[0] > Nx-1.0f; // rightmost tile

  if(!(top || bot)) return;

  //-------------------------------------------------- 
  // NOTE: we only enter here if boundary tile

  //-------------------------------------------------- 
  // remove outflowing particles
  for(auto&& con : tile.containers) {
    for(size_t n=0; n<con.size(); n++) {
                                                                                      
      // left BC
      if(  con.loc(0,n) < x_left - 3.0f )                          con.info(n) = -1; // inside
      //if( (con.loc(0,n) < x_left - 2.5f ) && con.vel(0,n) < 0.0f ) con.info(n) = -1; // in-flowing
                                                                                    
      // right BC
      //if(  con.loc(0,n) > Nx -  1.0f ) con.info(n) = -1; // outside
      //if( (con.loc(0,n) > Nx - 10.0f ) && (con.vel(0,n) < 0.0f) )    con.info(n) = -1; // reflected
                                                                                       
    }
  }

  // remove outflowing prtcls; 
  // here the transfer storage is used as a tmp container for storing the indices
  for(auto&& con : tile.containers) {
    con.delete_transferred_particles(); 
  }

  return;
}


template<size_t D>
void pic::Gap<D>::drive_prtcls(
    pic::Tile<D>& tile)
{

  const float c = tile.cfl; // c = dt

  //-------------------------------------------------- 
  // add (negative) acceleration from star's gravity
  for(auto&& con : tile.containers) {

     // gravitational force is F_g = G M m/r^2 = g m
     // surface gravity is given as g = G M /r^2 
     // hence gravitational acceleration is a_grav = F_grav/m = G M/r^2 = g_0 (R/r)^2
    const float& m = con.m; // prtcl mass in units of m_e
    const float  g_m = grav_const*m;

    for(size_t n=0; n<con.size(); n++) {
                                                                                      
      auto& x = con.loc(0,n); // x location
      const float r = std::max(1.0f, 1.0f + (x - x_left)/gap_length); // r = R_star + height

      con.vel(0,n) -= c*g_m/r/r;
    }
  }

  return;
}





//--------------------------------------------------
// explicit template instantiation
template class pic::Gap<1>; // 1D only
