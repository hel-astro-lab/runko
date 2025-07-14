#include "core/emf/boundaries/conductor.h"
#include "tools/vector.h"

#include <cmath> 
#include <cassert>

using std::min;
using std::max;
using std::abs;
using std::sqrt;
using toolbox::Vec3;
using toolbox::norm;
using toolbox::norm1d;
using toolbox::norm2d;
using emf::StaggeredSphericalCoordinates;


// Decaying dipole-like field in 1D
template<>
Vec3<float> emf::Conductor<1>::dipole(Vec3<float>& xvec)
{
  Vec3<float> ret( 2.0/(pow(xvec(0), 3) + EPS), 0.0, 0.0); // real dipole
  //Vec3<float> ret( 2.0, 0.0, 0.0); // const
  return ret;
}

// General dipole formula in 2D cartesian coordinates
template<>
Vec3<float> emf::Conductor<2>::dipole(Vec3<float>& xvec)
{

  //non-rotating magnetic moment vector components
  float p1 = sin(chi_mu), p2 = cos(chi_mu); //, p3 = 0.0;   // 2D orientation

  // final rotating magnetic moment with phase included; rotates in xz plane
  // TODO rotation turned off for 2D; i.e., no phase dependency
  float mux = p1; //*cos(phase) - p3*sin(phase);
  float muy = p2;
  //float muz = p3; //*sin(phase) + p3*cos(phase);

  float rad = std::sqrt(xvec(0)*xvec(0) + xvec(1)*xvec(1));

  // mu . r
  float mudotr = mux*xvec(0) + muy*xvec(1);

  Vec3<float> ret(
        3.0*xvec(0)*mudotr/( pow(rad,5) + EPS) - mux/(pow(rad,3) + EPS), //x
        3.0*xvec(1)*mudotr/( pow(rad,5) + EPS) - muy/(pow(rad,3) + EPS), //y
        0.0); 

  return ret;
}


// General dipole formula in 3D cartesian coordinates
template<>
Vec3<float> emf::Conductor<3>::dipole(Vec3<float>& xvec)
{

  //non-rotating magnetic moment vector components
  float p1 = sin(chi_mu), p2 = 0.0, p3 = cos(chi_mu);

  // final rotating magnetic moment with phase included
  float mux = p1*cos(phase_mu) - p2*sin(phase_mu);
  float muy = p1*sin(phase_mu) + p2*cos(phase_mu);
  float muz = p3;

  float rad = norm(xvec); 

  // mu . r
  float mudotr = mux*xvec(0) + muy*xvec(1) + muz*xvec(2);

  Vec3<float> ret(
   3.0*xvec(0)*mudotr/( pow(rad,5) + EPS) - mux/(pow(rad,3) + EPS), //x
   3.0*xvec(1)*mudotr/( pow(rad,5) + EPS) - muy/(pow(rad,3) + EPS), //y
   3.0*xvec(2)*mudotr/( pow(rad,5) + EPS) - muz/(pow(rad,3) + EPS)  //z
   );

  return ret;
}



template<size_t D>
void emf::Conductor<D>::insert_em(
    emf::Tile<D>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;

  auto& gs = tile.get_grids();

  const float c = tile.cfl; // (numerical) speed of light
                              
  //--------------------------------------------------
  // angular velocity
  float Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check
                                  
  Vec3<float> Om;
  if(D == 1) Om.set(Omega, 0.0, 0.0); // Omega unit vector along x-axis
  if(D == 2) Om.set(0.0, Omega, 0.0); // Omega unit vector along y-axis
  if(D == 3) Om.set( sin(chi_om)*cos(phase_om)*Omega, sin(chi_om)*sin(phase_om)*Omega, cos(chi_om)*Omega ); 
  //if(D == 3) Om.set( sin(chi_om)*Omega, 0.0,   cos(chi_om)*Omega ); // general Omega unit vector 


  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  // loop indices
  int imin, imax, jmin, jmax, kmin, kmax;
  if(D == 1){
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin =  0, jmax = 1;
    kmin =  0, kmax = 1;
  } else if (D == 2) {
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin = -3, jmax = tile.mesh_lengths[1]+3;
    kmin =  0, kmax = 1;
  } else if (D == 3) {
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin = -3, jmax = tile.mesh_lengths[1]+3;
    kmin = -3, kmax = tile.mesh_lengths[2]+3;
  }


  for(int k=kmin; k<kmax; k++) 
  for(int j=jmin; j<jmax; j++) 
  for(int i=imin; i<imax; i++) {

    // global grid coordinates
    float iglob = (D>=1) ? i + mins[0] : 0;
    float jglob = (D>=2) ? j + mins[1] : 0;
    float kglob = (D>=3) ? k + mins[2] : 0;

    //--------------------------------------------------
    // magnetic field
    
    //--------------------------------------------------
    // bx
    auto r1  = coord.bx().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
    auto bxd = B0*dipole(r1); // diple field

    // by
    auto r2  = coord.by().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
    auto byd = B0*dipole(r2); // diple field

    // bz
    auto r3  = coord.bz().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
    auto bzd = B0*dipole(r3); // diple field

    gs.bx(i,j,k) = bxd(0);
    gs.by(i,j,k) = byd(1);
    gs.bz(i,j,k) = bzd(2);


    //--------------------------------------------------
    // special mode to set constant background field in 1D
    if( (D == 1) && set_const_b ) {

      // get value of the dipole at r=R
      auto r1_left = coord.bx().vec(0, 0, 0, D); // cartesian position vector in "star's coordinates"
      bxd = B0*dipole(r1_left); // dipole field at r=R

      gs.bx(i,j,k) = bxd(0);
      gs.by(i,j,k) = 0.0;
      gs.bz(i,j,k) = 0.0;
    }


    //--------------------------------------------------
    // electric field

    //--------------------------------------------------
    // x coord staggering for Ex
    //xr = (D>=1) ? coord.ex().x(iglob) : 0;
    //yr = (D>=2) ? coord.ex().y(jglob) : 0;
    //zr = (D>=3) ? coord.ex().z(kglob) : 0;
    //r.set(xr,yr,zr);

    //Bd(0) = B0*dipole(xr,yr,zr,0);
    //Bd(1) = B0*dipole(xr,yr,zr,1);
    //Bd(2) = B0*dipole(xr,yr,zr,2);

    //auto vrot1 = cross(Om, r);
    //auto erot1 = -1.0f*cross(vrot1, Bd);
    //exd = erot1(0); // x component

    ////--------------------------------------------------
    //// y coord staggering for Ey
    //xr = (D>=1) ? coord.ey().x(iglob) : 0;
    //yr = (D>=2) ? coord.ey().y(jglob) : 0;
    //zr = (D>=3) ? coord.ey().z(kglob) : 0;
    //r.set(xr,yr,zr);

    //Bd(0) = B0*dipole(xr,yr,zr,0);
    //Bd(1) = B0*dipole(xr,yr,zr,1);
    //Bd(2) = B0*dipole(xr,yr,zr,2);

    //auto vrot2 = cross(Om, r);
    //auto erot2 = -1.0f*cross(vrot2, Bd);
    //eyd = erot2(1); // y component

    ////--------------------------------------------------
    //// z coord staggering for Ez
    //xr = (D>=1) ? coord.ez().x(iglob) : 0;
    //yr = (D>=2) ? coord.ez().y(jglob) : 0;
    //zr = (D>=3) ? coord.ez().z(kglob) : 0;
    //r.set(xr,yr,zr);
    //        
    //Bd(0) = B0*dipole(xr,yr,zr,0);
    //Bd(1) = B0*dipole(xr,yr,zr,1);
    //Bd(2) = B0*dipole(xr,yr,zr,2);
    //        
    //auto vrot3 = cross(Om, r);
    //auto erot3 = -1.0f*cross(vrot3, Bd);
    //ezd = erot3(2); // z component

    //gs.ex(i,j,k) = exd;
    //gs.ey(i,j,k) = eyd;
    //gs.ez(i,j,k) = ezd;

  }

  //--------------------------------------------------
  // start magnetosphere from non-rotating configuration, i.e., set E_rot 
  if( (D == 1) ) {

    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      const float iglob = (D>=1) ? i + mins[0] : 0;
                                   
      // TODO define rotating electric field inside star
      auto r  = coord.ex().vec(iglob, 0.0f, 0.0f, D); // cartesian position vector in "star's coordinates"
      auto bd = B0*dipole(r); // diple field
      auto h  = abs(r(0)); // cylindrical coordinate system height

      //--------------------------------------------------
      // special mode to set constant background field in 1D
      if( (D == 1) && set_const_b ) { // get value of the dipole at r=R
        auto r1_left = coord.bx().vec(0, 0, 0, D); // cartesian position vector in "star's coordinates"
        bd = B0*dipole(r1_left); // dipole field at r=R
      }

                             
      auto s  = 1.0f - shape( h, radius, delta); // height smoothing parameter
      //auto s  = h < radius + 4 ? 1.0f : 0.0f; // step function; field inside star h < r_*
      //auto s  = h > radius + 4 ? 1.0f : 0.0f; // step function; field outside star h > r_*

      //const auto rcyl1 = radius_pc;
      //const float offs = delta_pc; // expanded polar cap
      //sx      *= shape(rcyl1, radius_pc + offs, delta_pc); // no damping from polar cap edges in 1D
                                   
      float vrot = Om(0)*radius_pc/c; //r1(0); // Omega x r_pc
      float erot = 1.0f*vrot*bd(0); //-v x B

      // linear voltage
      erot = iglob < radius_pc ? erot*( 1.0 - iglob/radius_pc ) : 0.0;

      const float erot1 = erot;
      const float erot2 = 0.0f;
      const float erot3 = 0.0f;

      //--------------------------------------------------
      // blending of old + new solution
      gs.ex(i,0,0) = s*erot1 + (1.0f - s)*gs.ex(i,0,0); 
      gs.ey(i,0,0) = s*erot2 + (1.0f - s)*gs.ey(i,0,0); 
      gs.ez(i,0,0) = s*erot3 + (1.0f - s)*gs.ez(i,0,0); 
    }
  }

}



template<size_t D>
void emf::Conductor<D>::update_b(
    emf::Tile<D>& tile)
{

  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  // loop indices
  int imin, imax, jmin, jmax, kmin, kmax;
  if(D == 1){
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin =  0, jmax = 1;
    kmin =  0, kmax = 1;
  } else if (D == 2) {
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin = -3, jmax = tile.mesh_lengths[1]+3;
    kmin =  0, kmax = 1;
  } else if (D == 3) {
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin = -3, jmax = tile.mesh_lengths[1]+3;
    kmin = -3, kmax = tile.mesh_lengths[2]+3;
  }


  auto& gs = tile.get_grids();

  //-------------------------------------------------- 
  // find if this is a corner tile
  bool left  = false;
  bool right = false;
  bool top   = false;
  bool bot   = false;
  bool front = false;
  bool back  = false;

  Vec3<float> x1, x2, x3, x4, x5, x6, x7, x8;
  if(D == 1){
    if( mins[0] < 1 )    bot   = true; 
    if( maxs[0] > Nx-1 ) top   = true; 

    x1 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x2 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
    x3 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x4 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x5 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
    x6 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
    x7 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x8 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
  } else if(D == 2){
    if( mins[1] < 1 )    bot   = true; 
    if( mins[0] < 1 )    left  = true; 
    if( maxs[1] > Ny-1 ) top   = true; 
    if( maxs[0] > Nx-1 ) right = true; 

    // all 8 of tile corners are compared to radius to find if we are inside the spherical region
    x1 = coord.mid().vec(mins[0], mins[1], 0.0, D); 
    x2 = coord.mid().vec(maxs[0], mins[1], 0.0, D); 
    x3 = coord.mid().vec(mins[0], maxs[1], 0.0, D); 
    x4 = coord.mid().vec(mins[0], mins[1], 0.0, D); 
    x5 = coord.mid().vec(maxs[0], maxs[1], 0.0, D); 
    x6 = coord.mid().vec(maxs[0], mins[1], 0.0, D); 
    x7 = coord.mid().vec(mins[0], maxs[1], 0.0, D); 
    x8 = coord.mid().vec(maxs[0], maxs[1], 0.0, D); 
  } else if (D == 3) {
    if( mins[0] < 1 )    left  = true; 
    if( maxs[0] > Nx-1 ) right = true; 
    if( mins[1] < 1 )    front = true; 
    if( maxs[1] > Ny-1 ) back  = true; 
    if( mins[2] < 1 )    bot   = true; 
    if( maxs[2] > Nz-1 ) top   = true; 
      
    // all 8 of tile corners are compared to radius to find if we are inside the spherical region
    x1 = coord.mid().vec(mins[0], mins[1], mins[2], D); 
    x2 = coord.mid().vec(maxs[0], mins[1], mins[2], D); 
    x3 = coord.mid().vec(mins[0], maxs[1], mins[2], D); 
    x4 = coord.mid().vec(mins[0], mins[1], maxs[2], D); 
    x5 = coord.mid().vec(maxs[0], maxs[1], mins[2], D); 
    x6 = coord.mid().vec(maxs[0], mins[1], maxs[2], D); 
    x7 = coord.mid().vec(mins[0], maxs[1], maxs[2], D); 
    x8 = coord.mid().vec(maxs[0], maxs[1], maxs[2], D); 
  }

  const int H = 2; // halo region size for nulling of boundaries

    
  float tile_len = 0.0;
  if( D == 1 ) tile_len = tile.mesh_lengths[0];
  if( D == 2 ) tile_len = tile.mesh_lengths[1];
  if( D == 3 ) tile_len = tile.mesh_lengths[2];

  //--------------------------------------------------
  // inside star

  // TODO this is not a bulletproof method to find if the star is inside tile;
  //      if a small star is in the middle of the tile it can be missed...
  bool inside_star = 
    norm(x1) < 1.1*radius ||
    norm(x2) < 1.1*radius ||
    norm(x3) < 1.1*radius ||
    norm(x4) < 1.1*radius ||
    norm(x5) < 1.1*radius ||
    norm(x6) < 1.1*radius ||
    norm(x7) < 1.1*radius ||
    norm(x8) < 1.1*radius;

  const float b_offset = 0.0; // height offset of b smoothing

  if( inside_star) {

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      float jglob = (D>=2) ? j + mins[1] : 0;
      float kglob = (D>=3) ? k + mins[2] : 0;

      //--------------------------------------------------
      // bx
      auto r1  = coord.bx().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
      auto bxd = B0*dipole(r1); // diple field
      auto h1  = abs(r1(D-1)); // cylindrical coordinate system height
      //auto sx  = shape( norm(r1), radius + b_offset, delta); // radial smoothing parameter
      auto sx  = shape( h1, radius + b_offset, delta); // radial smoothing parameter

      // by
      auto r2  = coord.by().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
      auto byd = B0*dipole(r2); // diple field
      auto h2  = abs(r2(D-1)); // cylindrical coordinate system height
      //auto sy  = shape( norm(r2), radius + b_offset, delta); // radial smoothing parameter
      auto sy  = shape( h2, radius + b_offset, delta); // radial smoothing parameter

      // bz
      auto r3  = coord.bz().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
      auto bzd = B0*dipole(r3); // diple field
      auto h3  = abs(r3(D-1)); // cylindrical coordinate system height
      //auto sz  = shape( norm(r3), radius + b_offset, delta); // radial smoothing parameter
      auto sz  = shape( h3, radius + b_offset, delta); // radial smoothing parameter

      //--------------------------------------------------
      // blending of old + new solution
      gs.bx(i,j,k) = sx*bxd(0) + (1.0f - sx)*gs.bx(i,j,k); 
      gs.by(i,j,k) = sy*byd(1) + (1.0f - sy)*gs.by(i,j,k); 
      gs.bz(i,j,k) = sz*bzd(2) + (1.0f - sz)*gs.bz(i,j,k); 
    }}}
  }

  //--------------------------------------------------
  // damping on a cylindrical region around the pcap
    
  // norm2d gives the length of the vector x and y components (ignoring z)
  // this gives us the cylindrical coordinate

  bool inside_cyl_bcs = false;
  float rbox = 0.0;

  // N/A for 1D; ignored

  if(D == 2) {
    rbox = 0.5*Nx - 0.5*tile.mesh_lengths[0]; // half of box - half tile
                               
    inside_cyl_bcs = 
      norm1d(x1) > rbox || norm1d(x2) > rbox || norm1d(x3) > rbox || norm1d(x4) > rbox ||
      norm1d(x5) > rbox || norm1d(x6) > rbox || norm1d(x7) > rbox || norm1d(x8) > rbox;
  } else if(D == 3) {
    rbox = 0.5*Nx-H-1; // half of box size in x direction (not incl halos)
                               
    inside_cyl_bcs = 
      norm2d(x1) > rbox || norm2d(x2) > rbox || norm2d(x3) > rbox || norm2d(x4) > rbox ||
      norm2d(x5) > rbox || norm2d(x6) > rbox || norm2d(x7) > rbox || norm2d(x8) > rbox;
  }

  if( (D>=2) && inside_cyl_bcs ) {

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

      // global grid coordinates
      float iglob = (D>=1) ? i + mins[0] : 0;
      float jglob = (D>=2) ? j + mins[1] : 0;
      float kglob = (D>=3) ? k + mins[2] : 0;


      //--------------------------------------------------
      // bx
      auto r1    = coord.bx().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
      auto bxd   = B0*dipole(r1); // diple field
      auto rcyl1 = (D == 2) ? norm1d(r1) : norm2d(r1); // cylindrical radius

      // by
      auto r2    = coord.by().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
      auto byd   = B0*dipole(r2); // diple field
      auto rcyl2 = (D == 2) ? norm1d(r2) : norm2d(r2); // cylindrical radius

      // bz
      auto r3    = coord.bz().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
      auto bzd   = B0*dipole(r3); // diple field
      auto rcyl3 = (D == 2) ? norm1d(r3) : norm2d(r3); // cylindrical radius

      //--------------------------------------------------
      // ver 1; tanh profile
      auto delta_ext = 0.25f*tile_len; // 1/4 of tile size
      auto sx = shape(rcyl1, rbox, delta_ext); 
      auto sy = shape(rcyl2, rbox, delta_ext);
      auto sz = shape(rcyl3, rbox, delta_ext); 

      //ver 4; exp profile
      //auto sx = 1.0f - min(1.0, exp( (rcyl1-rbox)/(3.0*tile_len) ) );
      //auto sy = 1.0f - min(1.0, exp( (rcyl2-rbox)/(3.0*tile_len) ) );
      //auto sz = 1.0f - min(1.0, exp( (rcyl3-rbox)/(3.0*tile_len) ) );

      //--------------------------------------------------
      // damp to dipole solution
      gs.bx(i,j,k) = sx*gs.bx(i,j,k) + (1.0f-sx)*bxd(0);
      gs.by(i,j,k) = sy*gs.by(i,j,k) + (1.0f-sy)*byd(1);
      gs.bz(i,j,k) = sz*gs.bz(i,j,k) + (1.0f-sz)*bzd(2);
    }}}
  }


  //--------------------------------------------------
  // top absorping layer

  if(top) {

    //float radius_ext = (D == 2) ? Ny - 0.5*tile_len : Nz - 0.5*tile_len;
      
    float radius_ext = 0.0;
    if(D == 1) radius_ext = Nx - 0.5*tile_len;
    if(D == 2) radius_ext = Ny - 0.5*tile_len;
    if(D == 3) radius_ext = Nz - 0.5*tile_len;

    //float radius_ext = (D == 2) ? Ny - H : Nz - H; // exp profile
    float delta_ext  = 0.25f*tile_len; // 1/4 of tile size
                                           
    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          const float iglob = (D>=1) ? i + mins[0] : 0;
          const float jglob = (D>=2) ? j + mins[1] : 0;
          const float kglob = (D>=3) ? k + mins[2] : 0;

          //auto h = (D == 2 ) ? jglob : kglob; // height
          const auto h = (D==1) ? iglob : (D==2) ? jglob : kglob; // height

          //--------------------------------------------------
          // bx
          auto r1  = coord.bx().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bxd = B0*dipole(r1); // diple field

          // by
          auto r2  = coord.by().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto byd = B0*dipole(r2); // diple field

          // bz
          auto r3  = coord.bz().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bzd = B0*dipole(r3); // diple field

          //--------------------------------------------------
          // ver 1; tanh profile
          auto s = shape(h, radius_ext, delta_ext); // tanh

          // ver2 linear profile
          //float radius_ext = (D == 2) ? Ny - tile_len : Nz - tile_len;
          //float delta_ext = 1.0*tile_len; // full tile size
          //s = 1.0f - std::max(0.0f, std::min(1.0f, (h-radius_ext)/delta_ext) ); //RELU

          //ver 3; mimic lambda profile from PML
          //float lam = pow( (h - radius_ext)/(1.0 - radius_ext), 3);
          //s = min(1.0f, lam);

          //ver 4; exp profile
          //auto s = 1.0f - min(1.0, exp( (h-radius_ext)/(3.0*tile_len) ) );

          //--------------------------------------------------
          // damp to dipole solution
          gs.bx(i,j,k) = s*gs.bx(i,j,k) + (1.0f-s)*bxd(0);
          gs.by(i,j,k) = s*gs.by(i,j,k) + (1.0f-s)*byd(1);
          gs.bz(i,j,k) = s*gs.bz(i,j,k) + (1.0f-s)*bzd(2);
    }}}
  }


  //--------------------------------------------------
  // bottom absorping layer

  if(bot) {

    float radius_ext = 0.5*tile_len; // exp profile
    float delta_ext  = 0.25f*tile_len; // 1/4 of tile size
                                           
    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          const float iglob = (D>=1) ? i + mins[0] : 0;
          const float jglob = (D>=2) ? j + mins[1] : 0;
          const float kglob = (D>=3) ? k + mins[2] : 0;

          auto rvec = coord.mid().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          const auto h = (D==1) ? iglob : (D==2) ? jglob : kglob; // height
          const auto rcyl = (D==1) ? 0.0f : (D==2) ? norm1d(rvec) : norm2d(rvec); // cylindrical radius

          //--------------------------------------------------
          // bx
          auto r1  = coord.bx().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bxd = B0*dipole(r1); // diple field

          // by
          auto r2  = coord.by().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto byd = B0*dipole(r2); // diple field

          // bz
          auto r3  = coord.bz().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bzd = B0*dipole(r3); // diple field

          //--------------------------------------------------
          // ver 1; tanh profile
          //auto s = 1.0f - shape(h, radius_ext, delta_ext); // tanh
          //if( rcyl < 1.5*radius_pc ) s = 1.0f; // act normal inside star

          // tanh profile outside the cylindrical region
          auto s = 1.0f ? rcyl < 1.5*radius_pc : 1.0f - shape(h, radius_ext, delta_ext);
                                                    
          //--------------------------------------------------
          // damp to dipole solution
          gs.bx(i,j,k) = s*gs.bx(i,j,k) + (1.0f-s)*bxd(0);
          gs.by(i,j,k) = s*gs.by(i,j,k) + (1.0f-s)*byd(1);
          gs.bz(i,j,k) = s*gs.bz(i,j,k) + (1.0f-s)*bzd(2);
    }}}

  }





  //--------------------------------------------------
  // box boundaries

  if( left  ||
      right || 
      front ||
      back  ||
      bot   ||
      top
    ){

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          float iglob = (D>=1) ? i + mins[0] : 0;
          float jglob = (D>=2) ? j + mins[1] : 0;
          float kglob = (D>=3) ? k + mins[2] : 0;

          //--------------------------------------------------
          // bx
          auto r1  = coord.bx().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bxd = B0*dipole(r1); // diple field

          // by
          auto r2  = coord.by().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto byd = B0*dipole(r2); // diple field

          // bz
          auto r3  = coord.bz().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bzd = B0*dipole(r3); // diple field

          //--------------------------------------------------
          // sides
          bool inside_bot   = (D==1) ? iglob < H    : (D==2) ? jglob < H    : kglob < H; // y or z direction flip 
          bool inside_top   = (D==1) ? iglob > Nx-1 : (D==2) ? jglob > Ny-1 : kglob > Nz-1; // y or z direction flip 

          bool inside_left  = left  && iglob < H; 
          bool inside_right = right && iglob > Nx-H-1; 

          bool inside_front = front && jglob < H; 
          bool inside_back  = back  && jglob > Ny-H-1; 

          // cartesian boundaries
          bool inside_box_bcs = inside_left  ||
                                inside_right ||
                                inside_front ||
                                inside_back  ||
                                inside_top   ||
                                inside_bot;

          //--------------------------------------------------
          // vector friendly application of boolean boundary
          auto s = static_cast<float>(inside_box_bcs);

          gs.bx(i,j,k) = s*bxd(0) + (1.0f-s)*gs.bx(i,j,k);
          gs.by(i,j,k) = s*byd(1) + (1.0f-s)*gs.by(i,j,k);
          gs.bz(i,j,k) = s*bzd(2) + (1.0f-s)*gs.bz(i,j,k);
    }}}
  }

}



template<size_t D>
void emf::Conductor<D>::update_e(
    emf::Tile<D>& tile)
{

  // angular velocity
  Vec3<float> Om; 
  float Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check

  const float c = tile.cfl; // (numerical) speed of light


  if(D == 1) Om.set(Omega, 0.0, 0.0); // Omega unit vector along y-axis
  if(D == 2) Om.set(0.0, Omega, 0.0); // Omega unit vector along y-axis
  if(D == 3) Om.set( sin(chi_om)*cos(phase_om)*Omega, sin(chi_om)*sin(phase_om)*Omega, cos(chi_om)*Omega ); 

  //if(D == 3) Om.set( sin(chi_om)*Omega, 0.0,   cos(chi_om)*Omega ); // general Omega unit vector 
    
  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& gs = tile.get_grids();

  // loop indices
  int imin, imax, jmin, jmax, kmin, kmax;
  if(D == 1){
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin =  0, jmax = 1;
    kmin =  0, kmax = 1;
  } else if (D == 2) {
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin = -3, jmax = tile.mesh_lengths[1]+3;
    kmin =  0, kmax = 1;
  } else if (D == 3) {
    imin = -3, imax = tile.mesh_lengths[0]+3;
    jmin = -3, jmax = tile.mesh_lengths[1]+3;
    kmin = -3, kmax = tile.mesh_lengths[2]+3;
  }

  //-------------------------------------------------- 
  // null sides to prevent periodic bc conditions
  bool left  = false;
  bool right = false;
  bool top   = false;
  bool bot   = false;
  bool front = false;
  bool back  = false;

  Vec3<float> x1, x2, x3, x4, x5, x6, x7, x8;
  if(D == 1){
    if( mins[0] < 1 )    bot   = true; 
    if( maxs[0] > Nx-1 ) top   = true; 

    x1 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x2 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
    x3 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x4 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x5 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
    x6 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
    x7 = coord.mid().vec(mins[0], 0.0, 0.0, D); 
    x8 = coord.mid().vec(maxs[0], 0.0, 0.0, D); 
  } else if(D == 2){
    if( mins[1] < 1 )    bot   = true; 
    if( mins[0] < 1 )    left  = true; 
    if( maxs[1] > Ny-1 ) top   = true; 
    if( maxs[0] > Nx-1 ) right = true; 

    // all 8 of tile corners are compared to radius to find if we are inside the spherical region
    x1 = coord.mid().vec(mins[0], mins[1], 0.0, D); 
    x2 = coord.mid().vec(maxs[0], mins[1], 0.0, D); 
    x3 = coord.mid().vec(mins[0], maxs[1], 0.0, D); 
    x4 = coord.mid().vec(mins[0], mins[1], 0.0, D); 
    x5 = coord.mid().vec(maxs[0], maxs[1], 0.0, D); 
    x6 = coord.mid().vec(maxs[0], mins[1], 0.0, D); 
    x7 = coord.mid().vec(mins[0], maxs[1], 0.0, D); 
    x8 = coord.mid().vec(maxs[0], maxs[1], 0.0, D); 
  } else if (D == 3) {
    if( mins[0] < 1 )    left  = true; 
    if( maxs[0] > Nx-1 ) right = true; 
    if( mins[1] < 1 )    front = true; 
    if( maxs[1] > Ny-1 ) back  = true; 
    if( mins[2] < 1 )    bot   = true; 
    if( maxs[2] > Nz-1 ) top   = true; 

    // all 8 of tile corners are compared to radius to find if we are inside the spherical region
    x1 = coord.mid().vec(mins[0], mins[1], mins[2], D); 
    x2 = coord.mid().vec(maxs[0], mins[1], mins[2], D); 
    x3 = coord.mid().vec(mins[0], maxs[1], mins[2], D); 
    x4 = coord.mid().vec(mins[0], mins[1], maxs[2], D); 
    x5 = coord.mid().vec(maxs[0], maxs[1], mins[2], D); 
    x6 = coord.mid().vec(maxs[0], mins[1], maxs[2], D); 
    x7 = coord.mid().vec(mins[0], maxs[1], maxs[2], D); 
    x8 = coord.mid().vec(maxs[0], maxs[1], maxs[2], D); 
  }

  const int H = 2;
    
  float tile_len = 0.0;
  if( D == 1 ) tile_len = tile.mesh_lengths[0];
  if( D == 2 ) tile_len = tile.mesh_lengths[1];
  if( D == 3 ) tile_len = tile.mesh_lengths[2];

  // tile lengths
  const int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  const int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  const int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;

  //--------------------------------------------------
  // inside star

  //--------------------------------------------------
  // additionally; null epar in the closed field line region

  // for dipole fields r/Rbc = sin^2(theta) where r and theta are in spherical coordinates 
  // tolerance of 0.95 here so that we start killing the epar earlier
  auto sint = radius_pc/radius; // sin\theta = R_pc/R_star
  auto Rbc  = radius/sint/sint;  


    

  // NOTE: do not use full tile boundaries w/ halos for epar removal
  if( D >= 2 ) { // epar removal only for D=2 and 3

  #pragma omp parallel for
  for(int k=0; k<nz_tile; k++) {
    for(int j=0; j<ny_tile; j++) {
      #pragma omp simd
      for(int i=0; i<nx_tile; i++) {
      
        // global grid coordinates
        const float iglob = (D>=1) ? i + mins[0] : 0;
        const float jglob = (D>=2) ? j + mins[1] : 0;
        const float kglob = (D>=3) ? k + mins[2] : 0;

        //--------------------------------------------------
        // spherical coordinates; TODO ignoring staggering
        auto rvec = coord.mid().vec(iglob, jglob, kglob, D); //cartesian radius vector in star's coords
        const auto rad = norm(rvec); // length of radius

        const auto rcyl = (D==1) ? abs(rvec(0)) : (D==2) ? norm1d(rvec) : norm2d(rvec); // cylindrical radius
                                                                                          
        //auto sint = rcycl/rad; // \sin\theta
        //auto eta = rad/Rbc; // dimensionless dipole coordinate radius

        //bool inside_closed_field_region = eta < 1.0*sint*sint;

        //--------------------------------------------------
        // vector friendly application of boolean boundary
        //auto s = static_cast<float>(inside_closed_field_region);

        //--------------------------------------------------
        // smoothing function
        //s = shape( rad, 1.0*Rbc*sint*sint, 50.0); // location of the open-closed field line bc

        // condition for closed field line region is
        //eta = sint*sint
        //rad/Rbc = (rcycl/rad)^2 = rcyl^2/rad^2
        //rcycl^2 = (rad^3/Rbc)
        //rcycl > sqrt( rad^3/Rbc ) = rad*sqrt(rad/Rbc)
        // Then, we can use a smoothing function similar to that in pcap radius:
        auto s = 1.0f - shape(rcyl, rad*sqrt(rad/Rbc), delta_pc);  

        //--------------------------------------------------
        // additional height dependent damping inside the atmosphere; overwrites other smoothing
        //bool inside_atmos  = norm(rvec) < 1.0*radius + 2;

        const auto h = (D==1) ? iglob : (D==2) ? jglob : kglob; // cylindrical coordinate system height
        bool inside_atmos  = h < 1.0*radius + 2;
          
        if( inside_atmos ) s = 1.0f - shape(rcyl, radius_pc, delta_pc);  
                                                         
        //if(inside_atmos) s = shape(h, radius, delta); // radial smoothing parameter
        //if(inside_atmos) s = 1.0f - shape(norm(rvec), radius, delta); // radial smoothing parameter


        //--------------------------------------------------
        // epar; ver 1
        auto exi = gs.ex(i,j,k);
        auto eyi = gs.ey(i,j,k);
        auto ezi = gs.ez(i,j,k);

        auto bxi = gs.bx(i,j,k);
        auto byi = gs.by(i,j,k);
        auto bzi = gs.bz(i,j,k);
        auto bn  = std::sqrt( bxi*bxi + byi*byi + bzi*bzi ) + EPS;

        // E_\parallel
        auto epar = (exi*bxi + eyi*byi + ezi*bzi)/bn;

        //--------------------------------------------------
        // take out eparallel component from electric field
        auto exnew = exi - epar*bxi/bn;
        auto eynew = eyi - epar*byi/bn;
        auto eznew = ezi - epar*bzi/bn;

        // blend solution in with a smoothing function
        gs.ex(i,j,k) = s*exnew + (1.0f - s)*exi;
        gs.ey(i,j,k) = s*eynew + (1.0f - s)*eyi;
        gs.ez(i,j,k) = s*eznew + (1.0f - s)*ezi;

        //--------------------------------------------------
        // epar; ver 2 with interpolation
        //float exi, eyi, ezi, bxi, byi, bzi, b2, eparb;

        //const size_t iy = D >= 2 ? gs.ex.indx(0,1,0) - gs.ex.indx(0,0,0) : 0;
        //const size_t iz = D >= 3 ? gs.ex.indx(0,0,1) - gs.ex.indx(0,0,0) : 0;
        //const size_t ind = gs.ex.indx(i,j,k);

        //// ex
        //exi = gs.ex(ind);
        //eyi = 0.25 *(gs.ey(ind)    + gs.ey(ind+1)     + gs.ey(ind-iy)      + gs.ey(ind+1-iy));
        //ezi = 0.25 *(gs.ez(ind)    + gs.ez(ind+1)     + gs.ez(ind-iz)      + gs.ez(ind+1-iz));

        //bxi = 0.125*(gs.bx(ind)    + gs.bx(ind-iy)    + gs.bx(ind+1-iy)    + gs.bx(ind+1) +
        //             gs.bx(ind-iz) + gs.bx(ind-iy-iz) + gs.bx(ind+1-iy-iz) + gs.bx(ind+1-iz));
        //byi = 0.5 * (gs.by(ind)    + gs.by(ind-iz));
        //bzi = 0.5 * (gs.bz(ind)    + gs.bz(ind-iy));

        //b2  = bxi*bxi + byi*byi + bzi*bzi + EPS;
        //eparb = exi*bxi + eyi*byi + ezi*bzi;
        //auto exnew = exi-eparb*bxi/b2;

        //gs.ex(i,j,k) = s*exnew + (1.0f - s)*exi;

        ////--------------------------------------------------
        //// ey
        //exi = 0.25 * (gs.ex(ind)    + gs.ex(ind-1)    + gs.ex(ind+iy)      + gs.ex(ind-1+iy));
        //eyi =         gs.ey(ind);
        //ezi = 0.25 * (gs.ez(ind)    + gs.ez(ind+iy)   + gs.ez(ind-iz)      + gs.ez(ind+iy-iz));

        //bxi = 0.5 * ( gs.bx(ind)    + gs.bx(ind-iz));
        //byi = 0.125*( gs.by(ind)    + gs.by(ind-1)    + gs.by(ind-1+iy)    + gs.by(ind+iy) +
        //              gs.by(ind-iz) + gs.by(ind-1-iz) + gs.by(ind-1+iy-iz) + gs.by(ind+iy-iz));
        //bzi = 0.5 * ( gs.bz(ind)    + gs.bz(ind-1));

        //b2  = bxi*bxi + byi*byi + bzi*bzi + EPS;
        //eparb = exi*bxi + eyi*byi + ezi*bzi;
        //auto eynew = eyi-eparb*byi/b2;
        //gs.ey(i,j,k) = s*eynew + (1.0f - s)*eyi;

        ////--------------------------------------------------
        ////ez
        //exi = 0.25 * (gs.ex(ind)    + gs.ex(ind-1)   +  gs.ex(ind+iz)      + gs.ex(ind-1+iz ));
        //eyi = 0.25 * (gs.ey(ind)    + gs.ey(ind-iy)  +  gs.ey(ind+iz)      + gs.ey(ind-iy+iz));
        //ezi =         gs.ez(ind);

        //bxi = 0.5 * ( gs.bx(ind)    + gs.bx(ind-iy));
        //byi = 0.5 * ( gs.by(ind)    + gs.by(ind-1));
        //bzi = 0.125*( gs.bz(ind)    + gs.bz(ind-1)    + gs.bz(ind-1-iy)    + gs.bz(ind-iy) + 
        //              gs.bz(ind+iz) + gs.bz(ind-1+iz) + gs.bz(ind-1-iy+iz) + gs.bz(ind-iy+iz));

        //b2  = bxi*bxi + byi*byi + bzi*bzi + EPS;
        //eparb = exi*bxi + eyi*byi + ezi*bzi;
        //auto eznew = ezi - eparb*bzi/b2;
        //gs.ez(i,j,k) = s*eznew + (1.0f - s)*ezi;

  }}}

  } // end of D>=2



  //--------------------------------------------------
  // TODO this is not a bulletproof method to find if the star is inside tile;
  //      if a small star is in the middle of the tile it can be missed...
  bool inside_star = 
    norm(x1) < 1.1*radius ||
    norm(x2) < 1.1*radius ||
    norm(x3) < 1.1*radius ||
    norm(x4) < 1.1*radius ||
    norm(x5) < 1.1*radius ||
    norm(x6) < 1.1*radius ||
    norm(x7) < 1.1*radius ||
    norm(x8) < 1.1*radius;

  if( inside_star && (D == 1) ) {

    #pragma omp simd
    for(int i=imin; i<imax; i++) {

      // global grid coordinates
      const float iglob = (D>=1) ? i + mins[0] : 0;
                                   
      // rotating electric field inside star
      auto r  = coord.ex().vec(iglob, 0.0f, 0.0f, D); // cartesian position vector in "star's coordinates"
      auto bd = B0*dipole(r); // diple field
      auto h  = abs(r(0)); // cylindrical coordinate system height
                             
      auto s  = shape( h, radius, delta); // height smoothing parameter
      //auto s  = h < radius + 4 ? 1.0f : 0.0f; // step function; field inside star h < r_*
      //auto s  = h > radius + 4 ? 1.0f : 0.0f; // step function; field outside star h > r_*
                                   
      float vrot = Om(0)*radius_pc/c; //r1(0); // Omega x r_pc
      float erot = 1.0f*vrot*bd(0); //-v x B

      const float erot1 = 0.0f; // erot; // B_\parallel direction
      const float erot2 = 0.0f; 
      const float erot3 = 0.0f;

      //--------------------------------------------------
      // blending of old + new solution
      gs.ex(i,0,0) = s*erot1 + (1.0f - s)*gs.ex(i,0,0); 
      gs.ey(i,0,0) = s*erot2 + (1.0f - s)*gs.ey(i,0,0); 
      gs.ez(i,0,0) = s*erot3 + (1.0f - s)*gs.ez(i,0,0); 
    }


  } else if( inside_star && (D>=2) ) {

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          const float iglob = (D>=1) ? i + mins[0] : 0;
          const float jglob = (D>=2) ? j + mins[1] : 0;
          const float kglob = (D>=3) ? k + mins[2] : 0;

          const float offs = delta_pc; // expanded polar cap

          //--------------------------------------------------
          // ex
          auto r1  = coord.ex().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bxd = B0*dipole(r1); // diple field
          auto h1  = abs(r1(D-1)); // cylindrical coordinate system height
          //auto sx  = shape( norm(r1), radius, delta); // radial smoothing parameter
          auto sx  = shape( h1, radius, delta); // height smoothing parameter
                                                             
          auto rcyl1 = (D==1) ? 0.0f : (D==2) ? norm1d(r1) : norm2d(r1); // cylindrical radius
          sx        *= shape(rcyl1, radius_pc + offs, delta_pc); // damp off edges of polar cap
                                                             
          auto vrot1 = cross(Om, r1); // Omega x r
          auto erot1 = -1.0f*cross(vrot1, bxd); //-v x B


          //--------------------------------------------------
          // ey
          auto r2  = coord.ey().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto byd = B0*dipole(r2); // diple field
          auto h2  = abs(r2(D-1)); // cylindrical coordinate system height

          auto sy  = shape( h2, radius, delta); // height smoothing parameter

          auto rcyl2 = (D==1) ? 0.0f : (D==2) ? norm1d(r2) : norm2d(r2); // cylindrical radius
          sy      *= shape(rcyl2, radius_pc + offs, delta_pc); // damp off edges of polar cap

          auto vrot2 = cross(Om, r2); // Omega x r
          auto erot2 = -1.0f*cross(vrot2, byd); //-v x B


          //--------------------------------------------------
          // ez
          auto r3  = coord.ez().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto bzd = B0*dipole(r3); // diple field
          auto h3  = abs(r3(D-1)); // cylindrical coordinate system height
          //auto sz  = shape( norm(r3), radius, delta); // radial smoothing parameter
          auto sz  = shape( h3, radius, delta); // height smoothing parameter
                                                             
          auto rcyl3 = (D==1) ? 0.0f : (D==2) ? norm1d(r3) : norm2d(r3); // cylindrical radius
          sz        *= shape(rcyl3, radius_pc + offs, delta_pc); // damp off edges of polar cap
                                                             
          auto vrot3 = cross(Om, r3); // Omega x r
          auto erot3 = -1.0f*cross(vrot3, bzd); //-v x B


          //--------------------------------------------------
          // blending of old + new solution
          gs.ex(i,j,k) = sx*erot1(0) + (1.0f - sx)*gs.ex(i,j,k); 
          gs.ey(i,j,k) = sy*erot2(1) + (1.0f - sy)*gs.ey(i,j,k); 
          gs.ez(i,j,k) = sz*erot3(2) + (1.0f - sz)*gs.ez(i,j,k); 

    }}}
  }

  //--------------------------------------------------
  // damping on a cylindrical region around the pcap
    
  // norm2d gives the length of the vector x and y components (ignoring z)
  // this gives us the cylindrical coordinate
    
                              
  bool inside_cyl_bcs = false;
  float rbox = 0.0;

  if(D == 2) {
    rbox = 0.5*Nx - 0.5*tile_len; // half of box - half tile
                                        
    inside_cyl_bcs = 
      norm1d(x1) > rbox || norm1d(x2) > rbox || norm1d(x3) > rbox || norm1d(x4) > rbox ||
      norm1d(x5) > rbox || norm1d(x6) > rbox || norm1d(x7) > rbox || norm1d(x8) > rbox;
  } else if(D == 3) {
    rbox = 0.5*Nx-H-1; // half of box size in x direction (not incl halos)

    inside_cyl_bcs = 
      norm2d(x1) > rbox || norm2d(x2) > rbox || norm2d(x3) > rbox || norm2d(x4) > rbox ||
      norm2d(x5) > rbox || norm2d(x6) > rbox || norm2d(x7) > rbox || norm2d(x8) > rbox;
  }


  if( (D>=2) && inside_cyl_bcs ) {

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          float iglob = (D>=1) ? i + mins[0] : 0;
          float jglob = (D>=2) ? j + mins[1] : 0;
          float kglob = (D>=3) ? k + mins[2] : 0;

          //--------------------------------------------------
          auto rvec = coord.mid().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          auto rcyl = (D == 2) ? norm1d(rvec) : norm2d(rvec); // cylindrical radius

          //--------------------------------------------------
          // ver 1; tanh profile
          auto delta_ext = 0.25f*tile_len; // 1/4 of tile size
          auto s = shape(rcyl, rbox, delta_ext); 

          // ver2 linear profile
          //float radius_ext = (D == 2) ? Ny - tile_len : Nz - tile_len;
          //float delta_ext = 1.0*tile_len; // full tile size
          //s = 1.0f - std::max(0.0f, std::min(1.0f, (h-radius_ext)/delta_ext) ); //RELU

          //ver 3; mimic lambda profile from PML
          //float lam = pow( (h - radius_ext)/(1.0 - radius_ext), 3);
          //s = min(1.0f, lam);

          //ver 4; exp profile
          //float radius_ext = rbox;
          //s = 1.0f - min(1.0, exp( (rcyl-radius_ext)/(3.0*tile_len) ) );

          //--------------------------------------------------
          // damp to vacuum
          gs.ex(i,j,k) = s*gs.ex(i,j,k) + (1.0f-s)*0.0f;
          gs.ey(i,j,k) = s*gs.ey(i,j,k) + (1.0f-s)*0.0f;
          gs.ez(i,j,k) = s*gs.ez(i,j,k) + (1.0f-s)*0.0f;
    }}}
  }

  //--------------------------------------------------
  // inside top

  if(top) {

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          const float iglob = (D>=1) ? i + mins[0] : 0;
          const float jglob = (D>=2) ? j + mins[1] : 0;
          const float kglob = (D>=3) ? k + mins[2] : 0;

          const auto h = (D==1) ? iglob : (D==2) ? jglob : kglob; // height

          //--------------------------------------------------
          // ver 1; tanh profile
          const float radius_ext = (D==1) ? Nx - 2.5*tile_len : (D==2) ? Ny - 0.5*tile_len : Nz - 0.5*tile_len;
          const float delta_ext = 0.25*tile_len; // 1/4 of tile size
          const auto s = shape(h, radius_ext, delta_ext); // tanh

          // ver2 linear profile
          //float radius_ext = (D == 2) ? Ny - tile_len : Nz - tile_len;
          //float delta_ext = 1.0*tile_len; // full tile size
          //s = 1.0f - std::max(0.0f, std::min(1.0f, (h-radius_ext)/delta_ext) ); //RELU

          //ver 3; mimic lambda profile from PML
          //float lam = pow( (h - radius_ext)/(1.0 - radius_ext), 3);
          //s = min(1.0f, lam);

          //ver 4; exp profile
          //float radius_ext = (D == 2) ? Ny - H : Nz - H;
          //auto s = 1.0f - min(1.0, exp( (h-radius_ext)/(3.0*tile_len) ) );


          //--------------------------------------------------
          // damp to vacuum
          gs.ex(i,j,k) = s*gs.ex(i,j,k) + (1.0f-s)*0.0f;
          gs.ey(i,j,k) = s*gs.ey(i,j,k) + (1.0f-s)*0.0f;
          gs.ez(i,j,k) = s*gs.ez(i,j,k) + (1.0f-s)*0.0f;

    }}}
  }


  //--------------------------------------------------
  // inside bottom / outside star
  if( bot && (D>=2) ) { // only for D=2,3

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          const float iglob = (D>=1) ? i + mins[0] : 0;
          const float jglob = (D>=2) ? j + mins[1] : 0;
          const float kglob = (D>=3) ? k + mins[2] : 0;

          auto rvec = coord.mid().vec(iglob, jglob, kglob, D); // cartesian position vector in "star's coordinates"
          const auto h = (D==1) ? iglob : (D==2) ? jglob : kglob; // height
          const auto rcyl = (D==1) ? 0.0f : (D==2) ? norm1d(rvec) : norm2d(rvec); // cylindrical radius

          //--------------------------------------------------
          // ver 1; tanh profile
          const float radius_ext = 0.5*tile_len;
          const float delta_ext = 0.25*tile_len; // 1/4 of tile size
          auto s = 1.0 - shape(h, radius_ext, delta_ext); // tanh

          if( rcyl < 1.5*radius_pc ) s = 1.0f; // act normal inside star

          //--------------------------------------------------
          // damp to vacuum
          gs.ex(i,j,k) = s*gs.ex(i,j,k) + (1.0f-s)*0.0f;
          gs.ey(i,j,k) = s*gs.ey(i,j,k) + (1.0f-s)*0.0f;
          gs.ez(i,j,k) = s*gs.ez(i,j,k) + (1.0f-s)*0.0f;
    }}}

  }



  //--------------------------------------------------
  // box boundaries

  // TODO optimize for 1D where -3,-2,-1 and nx-3,-2,-1 values can be set manually much faster

  if( left  ||
      right || 
      front ||
      back  ||
      bot   ||
      top
    ){

    #pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {

          // global grid coordinates
          float iglob = (D>=1) ? i + mins[0] : 0;
          float jglob = (D>=2) ? j + mins[1] : 0;
          float kglob = (D>=3) ? k + mins[2] : 0;

          //--------------------------------------------------
          // sides
          bool inside_bot   = (D==1) ? iglob < H    : (D==2) ? jglob < H    : kglob < H; // y or z direction flip 
          bool inside_top   = (D==1) ? iglob > Nx-1 : (D==2) ? jglob > Ny-1 : kglob > Nz-1; // y or z direction flip 

          bool inside_left  = left  && iglob < H; 
          bool inside_right = right && iglob > Nx-H-1; 

          bool inside_front = front && jglob < H; 
          bool inside_back  = back  && jglob > Ny-H-1; 

          // cartesian boundaries
          bool inside_box_bcs = inside_left  ||
                                inside_right ||
                                inside_front ||
                                inside_back  ||
                                inside_top   ||
                                inside_bot;

          //--------------------------------------------------
          // vector friendly application of boolean boundary
          auto s = static_cast<float>(inside_box_bcs);

          //--------------------------------------------------
          // damp to vacuum
          gs.ex(i,j,k) = s*0.0f + (1.0f - s)*gs.ex(i,j,k);
          gs.ey(i,j,k) = s*0.0f + (1.0f - s)*gs.ey(i,j,k);
          gs.ez(i,j,k) = s*0.0f + (1.0f - s)*gs.ez(i,j,k);

    }}}
  }

}



// current from moving frame
template<size_t D>
void emf::Conductor<D>::update_j(
  emf::Tile<D>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;
  auto& gs = tile.get_grids();


  // tile lengths
  const int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  const int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  const int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;
  //const int H = 2;

  // cutoff radius after which jm is not applied. 
  // NOTE: this needs to be less than update_e top cutoff; hence we set it to 0.25*tile (and E cutoff is 0.5*tile)
  float radius_ext = 0.0;
  if(D == 1) radius_ext = Nx - 0.25*nx_tile;
  if(D == 2) radius_ext = Ny - 0.25*ny_tile;
  if(D == 3) radius_ext = Nz - 0.25*nz_tile;

  float delta_ext  = 0.10f*nx_tile; // sharp enough profile so that at radius_ext + 2*delta_ext = 0

  //--------------------------------------------------
  const float Omega = 2.0*PI/period;
  const float v = Omega*radius_pc; // v_rot = Omega \times r_pc

  const float c = tile.cfl; // \Delta t

  // set current
  #pragma omp parallel for
  for(int k=0; k<nz_tile; k++) {
    for(int j=0; j<ny_tile; j++) {
      #pragma omp simd
      for(int i=0; i<nx_tile; i++) {

        // global grid coordinates
        float ig = (D>=1) ? i + mins[0] : 0;
        float jg = (D>=2) ? j + mins[1] : 0;
        float kg = (D>=3) ? k + mins[2] : 0;

        float jx=0.0f, jy=0.0f, jz=0.0f;

        //--------------------------------------------------
        // current in 1D moving frame v = v \hat{y}
        if( D == 1 ) {
          // j_m = 
          //x:  - v d( E_y)
          //y:  -d( v E_x) + d( v^2 B_z)
          //z:  0

          // jx
          float dx_ey_at_x = gs.ey(i,j,k) - gs.ey(i-1,j,k); // term2: partial_x(E_y) at i,j,k
          //float dx_ey_at_x = gs.ey(i+1,j,k) - gs.ey(i,j,k); // term2: partial_x(E_y) at i,j,k // BAD: oscillates
          //jx = -v*dx_ey_at_x; // IGNORED term2

          // jy
          float dx_v_ex_at_y  = v*gs.ex(i,j,k) - v*gs.ex(i-1,j,k); // term1: partial_y(v E_x) at i,j,k
          //float dx_v_ex_at_y  = v*gs.ex(i+1,j,k) - v*gs.ex(i,j,k); // term1: partial_y(v E_x) at i,j,k // BAD: oscillates
	    
          float dx_v2_bz_at_y = v*v*gs.bz(i,j,k) - v*v*gs.bz(i-1,j,k);   // term3: partial_x( v^2 B_z )
          //float dx_v2_bz_at_y = v*v*gs.bz(i+1,j,k) - v*v*gs.bz(i,j,k); // term3: partial_x (v^2 B_z) // UNTESTED
          //jy = -dx_v_ex_at_y + dx_v2_bz_at_y; // IGNORED term3
          jy = -dx_v_ex_at_y; // ONLY term1 active

          //jz
          jz = 0.0f;
        }
          
        // ver 1; tanh profile
        const auto h = (D==1) ? ig : (D==2) ? jg : kg; // height
        auto s = shape(h, radius_ext, delta_ext); // tanh: suppress rotational currents at the right edge of the box

        //--------------------------------------------------
        // add to the current
        gs.jx(i,j,k) += s*jx*c;
        gs.jy(i,j,k) += s*jy*c;
        gs.jz(i,j,k) += s*jz*c;
      }
    }
  }



}

//--------------------------------------------------
// explicit template instantiation
template class emf::Conductor<1>; // 1D
template class emf::Conductor<2>; // 2D
template class emf::Conductor<3>; // 3D
