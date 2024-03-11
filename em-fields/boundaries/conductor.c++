#include "conductor.h"
#include "../../tools/vector.h"

#include <cmath> 
#include <cassert>

using std::min;
using std::max;
using std::abs;
using std::sqrt;
using toolbox::Vec3;


// General dipole formula in 2D cartesian coordinates
template<>
double fields::Conductor<2>::dipole(
        double x,
        double y,
        double z,
        int dim) {

  //non-rotating magnetic moment vector components
  double p1 = sin(chi_mu), p2 = cos(chi_mu), p3 = 0.0;   // 2D orientation

  // final rotating magnetic moment with phase included; rotates in xz plane
  // TODO rotation turned off for 2D; i.e., no phase dependency
  double mux = p1; //*cos(phase) - p3*sin(phase);
  double muy = p2;
  double muz = p3; //*sin(phase) + p3*cos(phase);

  double rad = std::sqrt(x*x + y*y);

  // mu . r
  double mudotr = mux*x + muy*y;

   if     (dim == 0) return 3.0*x*mudotr/( pow(rad,5) + EPS) - mux/(pow(rad,3) + EPS); //x
   else if(dim == 1) return 3.0*y*mudotr/( pow(rad,5) + EPS) - muy/(pow(rad,3) + EPS); //y
   else if(dim == 2) return 0.0; 
   return 0.0;
}


// General dipole formula in 3D cartesian coordinates
template<>
double fields::Conductor<3>::dipole(
        double x,
        double y,
        double z,
        int dim) {

  //non-rotating magnetic moment vector components
  double p1 = sin(chi_mu), p2 = 0.0, p3 = cos(chi_mu);

  // final rotating magnetic moment with phase included
  double mux = p1*cos(phase_mu) - p2*sin(phase_mu);
  double muy = p1*sin(phase_mu) + p2*cos(phase_mu);
  double muz = p3;

  double rad = std::sqrt(x*x + y*y + z*z);

  // mu . r
  double mudotr = mux*x + muy*y + muz*z;

   if     (dim == 0) return 3.0*x*mudotr/( pow(rad,5) + EPS) - mux/(pow(rad,3) + EPS); //x
   else if(dim == 1) return 3.0*y*mudotr/( pow(rad,5) + EPS) - muy/(pow(rad,3) + EPS); //y
   else if(dim == 2) return 3.0*z*mudotr/( pow(rad,5) + EPS) - muz/(pow(rad,3) + EPS); //z
   return 0.0;
}


template<size_t D>
void fields::Conductor<D>::insert_em(
    fields::Tile<D>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  float_m bxd, byd, bzd, exd, eyd, ezd;
  float_m iglob, jglob, kglob;
  float_m xr,yr,zr;
  float_m vx,vy;
  Vec3<float_m> r, Bd; // tmp variables

  // angular velocity
  float_m Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check
                                  
  Vec3<float_m> Om;
  if(D == 2) Om.set(0.0,                Omega, 0.0); // Omega unit vector along y-axis
  //if(D == 3) Om.set( sin(chi_om)*Omega, 0.0,   cos(chi_om)*Omega ); // general Omega unit vector 
  if(D == 3) Om.set( sin(chi_om)*cos(phase_om)*Omega, sin(chi_om)*sin(phase_om)*Omega, cos(chi_om)*Omega ); 


  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);


  int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;
  
  for(int k=-3; k<nz_tile+3; k++) 
  for(int j=-3; j<ny_tile+3; j++) 
  for(int i=-3; i<nx_tile+3; i++) {

    // global grid coordinates
    iglob = (D>=1) ? i + mins[0] : 0;
    jglob = (D>=2) ? j + mins[1] : 0;
    kglob = (D>=3) ? k + mins[2] : 0;

    //--------------------------------------------------
    // magnetic field
    
    // x coord staggering for Bx
    xr = (D>=1) ? coord.bx().x(iglob) : 0;
    yr = (D>=2) ? coord.bx().y(jglob) : 0;
    zr = (D>=3) ? coord.bx().z(kglob) : 0;
    bxd = B0*dipole(xr,yr,zr,0);

    // y coord staggering for By
    xr = (D>=1) ? coord.by().x(iglob) : 0;
    yr = (D>=2) ? coord.by().y(jglob) : 0;
    zr = (D>=3) ? coord.by().z(kglob) : 0;
    byd = B0*dipole(xr,yr,zr,1);

    // z coord staggering for Bz
    xr = (D>=1) ? coord.bz().x(iglob) : 0;
    yr = (D>=2) ? coord.bz().y(jglob) : 0;
    zr = (D>=3) ? coord.bz().z(kglob) : 0;
    bzd = B0*dipole(xr,yr,zr,2);

    yee.bx(i,j,k) = bxd;
    yee.by(i,j,k) = byd;
    yee.bz(i,j,k) = bzd;

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

    //yee.ex(i,j,k) = exd;
    //yee.ey(i,j,k) = eyd;
    //yee.ez(i,j,k) = ezd;

  }
}



template<size_t D>
void fields::Conductor<D>::update_b(
    fields::Tile<D>& tile)
{

  float_m bxd, byd, bzd;
  float_m bxi, byi, bzi;
  float_m bxnew, bynew, bznew;

  float_m iglob, jglob, kglob;

  float_m xr0,yr0,zr0;
  float_m xr,yr,zr;

  float_m s,r;

  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  // smoothing scales for different components
  //
  // NOTE: keeping here for historical purposes; not used
  //
  //float_m delta = 1.0;  //from conductor.h
  //
  //float_m delta_erad   = 1.0*delta; 
  //float_m delta_eperp  = 0.5*delta;
  //float_m delta_brad   = 1.0*delta;
  //float_m delta_bperp  = 1.0*delta;

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  //-------------------------------------------------- 
  // null sides to prevent periodic bc conditions
  bool left  = false;
  bool right = false;
  bool top   = false;
  bool bot   = false;
  bool front = false;
  bool back  = false;

  if(D == 2){
    if( mins[1] < 1 )    bot   = true; 
    if( mins[0] < 1 )    left  = true; 
    if( maxs[1] > Ny-1 ) top   = true; 
    if( maxs[0] > Nx-1 ) right = true; 
  } else if (D == 3) {
    if( mins[0] < 1 )    left  = true; 
    if( maxs[0] > Nx-1 ) right = true; 
    if( mins[1] < 1 )    front = true; 
    if( maxs[1] > Ny-1 ) back  = true; 
    if( mins[2] < 1 )    bot   = true; 
    if( maxs[2] > Nz-1 ) top   = true; 
  }

  const int H = 2;

  //--------------------------------------------------
  // additionally, define quantities for the closed field line region
  float_m sint = radius_pc/radius; // sin\theta = R_pc/R_star
  float_m Rbc  = radius/sint/sint;  
  float_m rad, eta;

  //--------------------------------------------------
  // loop over grid
  int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;

  for(int k=-3; k<nz_tile+3; k++) 
  for(int j=-3; j<ny_tile+3; j++) 
  for(int i=-3; i<nx_tile+3; i++) {
      
    // global grid coordinates
    iglob = (D>=1) ? i + mins[0] : 0;
    jglob = (D>=2) ? j + mins[1] : 0;
    kglob = (D>=3) ? k + mins[2] : 0;

    // spherical coordinates 
    xr0 = (D>=1) ? coord.mid().x(iglob) : 0;
    yr0 = (D>=2) ? coord.mid().y(jglob) : 0;
    zr0 = (D>=3) ? coord.mid().z(kglob) : 0;

    //--------------------------------------------------
    // regional checks
    
    // check if we are inside star
    bool inside_star = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.1*radius;        // curved surface
    //bool inside_star  = (D==3) ? abs(zr0) <= 1.1*radius : abs(yr0) <= 1.1*radius + 0.0; // flat surface

    //--------------------------------------------------
    // closed field line zone
    //rad  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0);
    //sint = D == 2 ? abs(xr0)/rad : sqrt(xr0*xr0 + yr0*yr0)/rad; // sin\theta
    //eta  = rad/Rbc; // dimensionless dipole coordinate radius
    //bool closed_field_zone = eta < 0.9*sint*sint; // some tolerance for the boundaries
                                                    
    //std::cout << "update_b: z" << closed_field_zone 
    //          << " rad:" << rad 
    //          << " sint " << sint 
    //          << " eta " << eta 
    //          << "\n";

    //--------------------------------------------------
    // sides
    bool inside_bot   = (D == 2) ? jglob < H    : kglob < H; // y or z direction flip 
    bool inside_top   = (D == 2) ? jglob > Ny-1 : kglob > Nz-1; // y or z direction flip 

    bool inside_left  = left  && iglob < H; 
    bool inside_right = right && iglob > Nx-H-1; 

    bool inside_front = front && jglob < H; 
    bool inside_back  = back  && jglob > Ny-H-1; 

    // cartesian boundaries
    bool inside_box_bcs = inside_left  ||
                          inside_right ||
                          inside_front ||
                          inside_back;

    //--------------------------------------------------
    // cylindrical region
    float_m rcyl = (D == 2) ? abs(xr0) : sqrt(xr0*xr0 + yr0*yr0); // cylindrical radius
    float_m rbox = 0.5*Nx-H-1; // half of box size in x direction incl halos
    bool inside_cyl_bcs = rcyl > rbox;

    //--------------------------------------------------

    // operate inside the region only
    if( inside_star  ||
        inside_bot   ||
        top          ||  // selecting larger chunk for top since we damp outgoing waves there
        inside_box_bcs ||
        inside_cyl_bcs
      ) {

      //--------------------------------------------------
      // Bx: x coord staggering for Bx
      xr = (D>=1) ? coord.bx().x(iglob) : 0;
      yr = (D>=2) ? coord.bx().y(jglob) : 0;
      zr = (D>=3) ? coord.bx().z(kglob) : 0;
      r = std::sqrt( xr*xr + yr*yr + zr*zr ); // curved
      //r = (D == 2) ? abs(yr) : abs(zr); // flat 

      // alternative syntax  using internal toolbox::Vec3
      //xvec = coord.bx().vec(iglob, jglob, kglob); // cartesian position vector in "star's coordinates"
      //r = norm(xvec);


      // interior diple field
      bxd = B0*dipole(xr,yr,zr,0);
      //byd = B0*dipole(xr,yr,zr,1);
      //bzd = B0*dipole(xr,yr,zr,2);

      bxi = yee.bx(i,j,k);
      //byi = yee.by(i,j,k); // interpolate to this location
      //bzi = yee.bz(i,j,k); // interpolate to this location

      s = shape(r, radius, delta);
      bxnew = s*bxd + (1.0f-s)*bxi;

      // Bx radial  component 
      // NOTE i've sketched a possible code here to split B_r and B_\perp update
      //      Seems like it is not needed, though.
      //
      //s = shape(r, radius-offs_brad, delta_brad);
      //bxrad =     s*(bxd - (bxd*xr + byd*yr + bzd*zr)*xr/r/r)
      //        (1-s)*(bxi - (bxi*xr + byi*yr + bzi*zr)*xr/r/r);
      //
      //s = shape(r, radius-offs_bperp, delta_bperp);
      // Bx perp  component 
      //bxperp =      s*(bxd*xr + byd*yr + bzd*zr)*xr/r/r 
      //        + (1-s)*(bxi*xr + byi*yr + bzi*zr)*xr/r/r;


      //--------------------------------------------------
      // By: y coord staggering for By
      xr = (D>=1) ? coord.by().x(iglob) : 0;
      yr = (D>=2) ? coord.by().y(jglob) : 0;
      zr = (D>=3) ? coord.by().z(kglob) : 0;
      r = std::sqrt( xr*xr + yr*yr +zr*zr ); // curved
      //r = (D == 2) ? abs(yr) : abs(zr); // flat 

      // interior diple field
      byd = B0*dipole(xr,yr,zr,1);
      byi = yee.by(i,j,k);

      s = shape(r, radius, delta);
      bynew = s*byd  + (1.0f-s)*byi;

      //--------------------------------------------------
      // Bz: z coord staggering for Bz
      xr = (D>=1) ? coord.bz().x(iglob) : 0;
      yr = (D>=2) ? coord.bz().y(jglob) : 0;
      zr = (D>=3) ? coord.bz().z(kglob) : 0;
      r = std::sqrt( xr*xr + yr*yr +zr*zr ); // curved
      //r = (D == 2) ? abs(yr) : abs(zr); // flat 

      // interior diple field
      bzd = B0*dipole(xr,yr,zr,2);
      bzi = yee.bz(i,j,k);

      s = shape(r, radius, delta);
      bznew = s*bzd + (1.0-s)*bzi;

      //--------------------------------------------------
      if(!top) {

        yee.bx(i,j,k) = bxnew;
        yee.by(i,j,k) = bynew;
        yee.bz(i,j,k) = bznew;

      } else if(top) {
        // manual damping of outgoing waves
        // poor-man's version of PML absorbing boundary conditions

        float_m tile_len = (D == 2 ) ? tile.mesh_lengths[1] : tile.mesh_lengths[2];
        float_m h        = (D == 2 ) ? jglob : kglob; // height
                                                        
        // ver 1; tanh
        float_m radius_ext = (D == 2) ? Ny - 0.5*tile_len : Nz - 0.5*tile_len;
        float_m delta_ext = 0.25*tile_len; // 1/4 of tile size
        s = shape(h, radius_ext, delta_ext); // tanh

        // ver2
        //float_m radius_ext = (D == 2) ? Ny - tile_len : Nz - tile_len;
        //float_m delta_ext = 1.0*tile_len; // full tile size
        //s = 1.0f - std::max(0.0f, std::min(1.0f, (h-radius_ext)/delta_ext) ); //RELU

        //ver 3; mimic lambda profile from PML
        //float_m lam = pow( (h - radius_ext)/(1.0 - radius_ext), 3);
        //s = min(1.0f, lam);

        //ver 4; exp profile
        //float_m radius_ext = (D == 2) ? Ny - 3 : Nz - 3;
        //s = 1.0f - min(1.0, exp( (h-radius_ext)/(3.0*tile_len) ) );

        // damp to dipole solution
        yee.bx(i,j,k) = s*yee.bx(i,j,k) + (1.0f-s)*bxnew;
        yee.by(i,j,k) = s*yee.by(i,j,k) + (1.0f-s)*bynew;
        yee.bz(i,j,k) = s*yee.bz(i,j,k) + (1.0f-s)*bznew;
      }
    }
  }
}



template<size_t D>
void fields::Conductor<D>::update_e(
    fields::Tile<D>& tile)
{
  float_m exd, eyd, ezd;
  float_m bxd, byd, bzd;
  float_m exi, eyi, ezi;
  float_m exnew, eynew, eznew;

  float_m iglob, jglob, kglob;

  float_m xr0,yr0,zr0;
  float_m xr,yr,zr;

  float_m vx,vy;
  float_m s, rcyl, h;


  // angular velocity
  Vec3<float_m> r, Bd, Om; // tmp variables
                                  
  float_m Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check

  if(D == 2) Om.set(0.0,                Omega, 0.0); // Omega unit vector along y-axis
  //if(D == 3) Om.set( sin(chi_om)*Omega, 0.0,   cos(chi_om)*Omega ); // general Omega unit vector 
  if(D == 3) Om.set( sin(chi_om)*cos(phase_om)*Omega, sin(chi_om)*sin(phase_om)*Omega, cos(chi_om)*Omega ); 

  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  // smoothing scales for Brad/Bperp components
  // NOTE: keeping here for historical purposes; not used
  //
  //float_m delta = 2.0;  // from conductor.h
  //
  //float_m delta_erad   = 1.0*delta; 
  //float_m delta_eperp  = 0.5*delta;
  //float_m delta_brad   = 1.0*delta;
  //float_m delta_bperp  = 1.0*delta;

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  //--------------------------------------------------
  // values to calculate location of closed/open field line region
  //
  // for dipole fields r/Rbc = sin^2(theta) where r and theta are in spherical coordinates 
  // tolerance of 0.95 here so that we start killing the epar earlier
  float_m sint = radius_pc/radius; // sin\theta = R_pc/R_star
  float_m Rbc  = radius/sint/sint;  
  float_m rad, bn, bxi, byi, bzi, epar;


  //-------------------------------------------------- 
  // null sides to prevent periodic bc conditions
  bool left  = false;
  bool right = false;
  bool top   = false;
  bool bot   = false;
  bool front = false;
  bool back  = false;

  if(D == 2){
    if( mins[1] < 1 )    bot   = true; 
    if( mins[0] < 1 )    left  = true; 
    if( maxs[1] > Ny-1 ) top   = true; 
    if( maxs[0] > Nx-1 ) right = true; 
  } else if (D == 3) {
    if( mins[0] < 1 )    left  = true; 
    if( maxs[0] > Nx-1 ) right = true; 
    if( mins[1] < 1 )    front = true; 
    if( maxs[1] > Ny-1 ) back  = true; 
    if( mins[2] < 1 )    bot   = true; 
    if( maxs[2] > Nz-1 ) top   = true; 
  }

  const int H = 2;

  //--------------------------------------------------
  // loop over grid
  int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;

  for(int k=-3; k<nz_tile+3; k++) 
  for(int j=-3; j<ny_tile+3; j++) 
  for(int i=-3; i<nx_tile+3; i++) {

    // global grid coordinates
    iglob = (D>=1) ? i + mins[0] : 0;
    jglob = (D>=2) ? j + mins[1] : 0;
    kglob = (D>=3) ? k + mins[2] : 0;

    // spherical coordinates
    xr0 = (D>=1) ? coord.mid().x(iglob) : 0;
    yr0 = (D>=2) ? coord.mid().y(jglob) : 0;
    zr0 = (D>=3) ? coord.mid().z(kglob) : 0;

    // check if we are inside star
    //bool inside_star = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.05*radius;
    //bool inside_star = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= radius + 2.0*delta_pc + 1.0;
      
    bool inside_star  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.1*radius; // curved
    //bool inside_star  = (D==3) ? abs(zr0) <= 1.1*radius : abs(yr0) <= 1.0*radius; //flat

    //--------------------------------------------------
    // closed/open field line zone
    rad  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0);
    sint = D == 2 ? abs(xr0)/rad : sqrt(xr0*xr0 + yr0*yr0)/rad; // sin\theta
    float_m eta = rad/Rbc; // dimensionless dipole coordinate radius
    bool inside_closed_field_region = eta < 1.0*sint*sint;

    //--------------------------------------------------
    // sides
    bool inside_bot   = (D == 2) ? jglob < H    : kglob < H; // y or z direction flip 
    bool inside_top   = (D == 2) ? jglob > Ny-2 : kglob > Nz-2; // y or z direction flip 

    bool inside_left  = left  && iglob < H; 
    bool inside_right = right && iglob > Nx-H-1; 

    bool inside_front = front && jglob < H; 
    bool inside_back  = back  && jglob > Ny-H-1; 

    //--------------------------------------------------
    // cartesian boundaries
    bool inside_box_bcs = inside_left  ||
                          inside_right ||
                          inside_front ||
                          inside_back;

    //--------------------------------------------------
    //  cylindrical boundaries
    rcyl = (D == 2) ? abs(xr0) : sqrt(xr0*xr0 + yr0*yr0); // cylindrical radius
    float_m rbox = 0.5*Nx -H-1 ; // half of box size in x direction incl halos
                                                            
    bool inside_cyl_bcs = rcyl > rbox;
    //--------------------------------------------------


    //--------------------------------------------------
    if( inside_star ||
        top         
        //inside_closed_field_region
        ) {

      //-------------------------------------------------- 
      // ex
      xr = coord.ex().x(iglob);
      yr = coord.ex().y(jglob);
      zr = coord.ex().z(kglob);
      r.set(xr,yr,zr); // spherical radius
      rcyl = (D == 2) ? abs(xr) : sqrt(xr*xr + yr*yr); // cylindrical radius
      h    = (D == 2) ? abs(yr) : abs(zr); // height from center

      Bd(0) = B0*dipole(xr,yr,zr,0);
      Bd(1) = B0*dipole(xr,yr,zr,1);
      Bd(2) = B0*dipole(xr,yr,zr,2);

      auto vrot1 = cross(Om, r);
      auto erot1 = -1.0f*cross(vrot1, Bd);
      exd = erot1(0); // x component

      exi = yee.ex(i,j,k);

      s =  shape(norm(r), radius, delta); //curved
      //s =  shape(h, radius, delta); // flat
      s *= shape(rcyl, radius_pc, delta_pc); // damp off edges of polar cap

      exnew = s*exd  + (1.0f-s)*exi;


      //-------------------------------------------------- 
      // ey
      xr = coord.ey().x(iglob);
      yr = coord.ey().y(jglob);
      zr = coord.ey().z(kglob);
      r.set(xr,yr,zr);
      rcyl = (D == 2) ? abs(xr) : sqrt(xr*xr + yr*yr); 
      h    = (D == 2) ? abs(yr) : abs(zr); // height from center

      Bd(0) = B0*dipole(xr,yr,zr,0);
      Bd(1) = B0*dipole(xr,yr,zr,1);
      Bd(2) = B0*dipole(xr,yr,zr,2);

      auto vrot2 = cross(Om, r);
      auto erot2 = -1.0f*cross(vrot2, Bd);
      eyd = erot2(1); // y component

      eyi = yee.ey(i,j,k);

      s =  shape(norm(r), radius, delta); //curved
      //s =  shape(h, radius, delta); // flat // flat
      s *= shape(rcyl, radius_pc, delta_pc); // damp off edges of polar cap

      eynew = s*eyd  + (1.0f-s)*eyi;

      //-------------------------------------------------- 
      // ez
      xr = coord.ez().x(iglob);
      yr = coord.ez().y(jglob);
      zr = coord.ez().z(kglob);
      r.set(xr,yr,zr);
      rcyl = (D == 2) ? abs(xr) : sqrt(xr*xr + yr*yr); 
      h    = (D == 2) ? abs(yr) : abs(zr); // height from center
              
      Bd(0) = B0*dipole(xr,yr,zr,0);
      Bd(1) = B0*dipole(xr,yr,zr,1);
      Bd(2) = B0*dipole(xr,yr,zr,2);
              
      auto vrot3 = cross(Om, r);
      auto erot3 = -1.0f*cross(vrot3, Bd);
      ezd = erot3(2); // z component

      ezi = yee.ez(i,j,k);

      s =  shape(norm(r), radius, delta); // curved
      //s =  shape(h, radius, delta); // flat // flat
      s *= shape(rcyl, radius_pc, delta_pc); // damp off edges of polar cap

      eznew = s*ezd  + (1.0f-s)*ezi;

      //std::cout << "  ijk" << i << "," << j << "," << k << " r:" << norm(r) << " rad:" << radius << " d:" << delta << " s:" << s << "\n";

      //--------------------------------------------------
      if(!top) {

        yee.ex(i,j,k) = exnew;
        yee.ey(i,j,k) = eynew;
        yee.ez(i,j,k) = eznew;

      } else if(top) {
        // manual damping of outgoing waves 
        // poor-man's version of PML absorbing boundary conditions

        float_m tile_len = (D == 2 ) ? tile.mesh_lengths[1] : tile.mesh_lengths[2];
        float_m h        = (D == 2 ) ? jglob : kglob; // height
                                                  
        // ver 1; tanh
        float_m radius_ext = (D == 2) ? Ny - 0.5*tile_len : Nz - 0.5*tile_len;
        float_m delta_ext = 0.25*tile_len; // 1/4 of tile size
        s = shape(h, radius_ext, delta_ext); // tanh

        // ver2; linear
        //float_m radius_ext = (D == 2) ? Ny - tile_len : Nz - tile_len;
        //float_m delta_ext = 1.0*tile_len; // full tile size
        //s = 1.0f - std::max(0.0f, std::min(1.0f, (h-radius_ext)/delta_ext) ); //RELU

        //ver 3; mimic lambda profile from PML
        //float_m lam = pow( (h - radius_ext)/(1.0 - radius_ext), 3);
        //s = min(1.0f, lam);

        //ver 4; exp profile
        //float_m radius_ext = (D == 2) ? Ny - 3 : Nz - 3;
        //s = 1.0f - min(1.0, exp( (h-radius_ext)/(3.0*tile_len) ) );

        // damp to vacuum
        //yee.ex(i,j,k) = s*yee.ex(i,j,k) + (1.0f-s)*0.0;
        //yee.ey(i,j,k) = s*yee.ey(i,j,k) + (1.0f-s)*0.0;
        //yee.ez(i,j,k) = s*yee.ez(i,j,k) + (1.0f-s)*0.0;
          
        // damp to co-rotation // TODO debug
        //yee.ex(i,j,k) = s*yee.ex(i,j,k) + (1.0f-s)*exd;
        //yee.ey(i,j,k) = s*yee.ey(i,j,k) + (1.0f-s)*eyd;
        //yee.ez(i,j,k) = s*yee.ez(i,j,k) + (1.0f-s)*ezd;

        // null currents; not deposited yet so no effect
        yee.jx(i,j,k) = s*yee.jx(i,j,k); 
        yee.jy(i,j,k) = s*yee.jy(i,j,k); 
        yee.jz(i,j,k) = s*yee.jz(i,j,k); 
      }

    }


    // boundaries
    if( inside_bot     ||
        inside_top     ||
        inside_box_bcs ||
        inside_cyl_bcs 
      ) {

      //-------------------------------------------------- 
      // remove e-field from box boundaries
      yee.ex(i,j,k) = 0.0;
      yee.ey(i,j,k) = 0.0;
      yee.ez(i,j,k) = 0.0;

    }

    //  // TODO DEBUG set e-field to co-rotation value 

    //  //-------------------------------------------------- 
    //  // ex
    //  xr = coord.ex().x(iglob);
    //  yr = coord.ex().y(jglob);
    //  zr = coord.ex().z(kglob);
    //  r.set(xr,yr,zr); // spherical radius

    //  Bd(0) = B0*dipole(xr,yr,zr,0);
    //  Bd(1) = B0*dipole(xr,yr,zr,1);
    //  Bd(2) = B0*dipole(xr,yr,zr,2);

    //  auto vrot1 = cross(Om, r);
    //  auto erot1 = -1.0f*cross(vrot1, Bd);
    //  exd = erot1(0); // x component

    //  //-------------------------------------------------- 
    //  // ey
    //  xr = coord.ey().x(iglob);
    //  yr = coord.ey().y(jglob);
    //  zr = coord.ey().z(kglob);
    //  r.set(xr,yr,zr);

    //  Bd(0) = B0*dipole(xr,yr,zr,0);
    //  Bd(1) = B0*dipole(xr,yr,zr,1);
    //  Bd(2) = B0*dipole(xr,yr,zr,2);

    //  auto vrot2 = cross(Om, r);
    //  auto erot2 = -1.0f*cross(vrot2, Bd);
    //  eyd = erot2(1); // y component

    //  //-------------------------------------------------- 
    //  // ez
    //  xr = coord.ez().x(iglob);
    //  yr = coord.ez().y(jglob);
    //  zr = coord.ez().z(kglob);
    //  r.set(xr,yr,zr);
    //          
    //  Bd(0) = B0*dipole(xr,yr,zr,0);
    //  Bd(1) = B0*dipole(xr,yr,zr,1);
    //  Bd(2) = B0*dipole(xr,yr,zr,2);
    //          
    //  auto vrot3 = cross(Om, r);
    //  auto erot3 = -1.0f*cross(vrot3, Bd);
    //  ezd = erot3(2); // z component
    //                    
    //  yee.ex(i,j,k) = exd;
    //  yee.ey(i,j,k) = eyd;
    //  yee.ez(i,j,k) = ezd;

    //}
  }


  //--------------------------------------------------
  // additionally; null epar in the closed field line region

  for(int k=-3; k<nz_tile+3; k++) 
  for(int j=-3; j<ny_tile+3; j++) 
  for(int i=-3; i<nx_tile+3; i++) {
      
    // global grid coordinates
    iglob = (D>=1) ? i + mins[0] : 0;
    jglob = (D>=2) ? j + mins[1] : 0;
    kglob = (D>=3) ? k + mins[2] : 0;

    // spherical coordinates; TODO ignoring staggering
    xr0 = (D>=1) ? coord.mid().x(iglob) : 0;
    yr0 = (D>=2) ? coord.mid().y(jglob) : 0;
    zr0 = (D>=3) ? coord.mid().z(kglob) : 0;

    rad  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0);
    sint = D == 2 ? abs(xr0)/rad : sqrt(xr0*xr0 + yr0*yr0)/rad; // sin\theta

    float_m eta = rad/Rbc; // dimensionless dipole coordinate radius

    //if( eta < 1.4*sint*sint) { // TODO smooth or not?
    if( eta < 1.0*sint*sint) {
      exi = yee.ex(i,j,k);
      eyi = yee.ey(i,j,k);
      ezi = yee.ez(i,j,k);

      bxi = yee.bx(i,j,k);
      byi = yee.by(i,j,k);
      bzi = yee.bz(i,j,k);
      bn  = std::sqrt( bxi*bxi + byi*byi + bzi*bzi ) + EPS;

      // E_\parallel
      epar = (exi*bxi + eyi*byi + ezi*bzi)/bn;

      // take out eparallel component from electric field
      exnew = exi - epar*bxi/bn;
      eynew = eyi - epar*byi/bn;
      eznew = ezi - epar*bzi/bn;

      // smoothing function
      //s = shape( rad, 1.05*Rbc*sint*sint, 50.0); // location of the open-closed field line bc
      s = 1.0; // no smoothing along the polar cap rims

      // blend solution in with a smoothing function
      yee.ex(i,j,k) = s*exnew + (1.0f - s)*exi;
      yee.ey(i,j,k) = s*eynew + (1.0f - s)*eyi;
      yee.ez(i,j,k) = s*eznew + (1.0f - s)*ezi;
    }
  }


}


// helper script to iterate yee container
inline void iterate_yee(
    fields::YeeLattice& yee,
    int imin,
    int imax,
    int jmin,
    int jmax,
    int kmin,
    int kmax,
    const float_m val,
    int mode)
{

  if( mode == 0 ) {

    for(int k=kmin; k<kmax; k++) 
    for(int j=jmin; j<jmax; j++) 
    for(int i=imin; i<imax; i++) {
      yee.ex(i,j,k) = val;
      yee.ey(i,j,k) = val;
      yee.ez(i,j,k) = val;
    }

  } else if( mode == 1 ) {

    for(int k=kmin; k<kmax; k++) 
    for(int j=jmin; j<jmax; j++) 
    for(int i=imin; i<imax; i++) {
      yee.bx(i,j,k) = val;
      yee.by(i,j,k) = val;
      yee.bz(i,j,k) = val;
    }
  }

}

template<>
void fields::Conductor<2>::null_edges(
    fields::Tile<2>& tile,
    int mode
    )
{
  const float_m eval = 0.0f; // edge value

  auto mins = tile.mins;
  auto maxs = tile.maxs;
  
  auto& yee = tile.get_yee();

  const int H = 2;
    
  // 2D
  const int kmin = 0;
  const int kmax = 1;
  int imin, imax, jmin, jmax;

  //-------------------------------------------------- 
  // null sides to prevent periodic bc conditions
  bool left  = false;
  bool right = false;
  bool top   = false;
  bool bot   = false;
  
  if( mins[1] < 1 )    bot   = true; 
  if( mins[0] < 1 )    left  = true; 
  if( maxs[1] > Ny-1 ) top   = true; 
  if( maxs[0] > Nx-1 ) right = true; 


  //--------------------------------------------------
  if(bot) {
    //std::cout << "bottom " << mins[0] << " " << mins[1] << " " << maxs[0] << " " << maxs[1] << std::endl;

    jmin = -H;
    jmax =  0;

    imin = -H;
    imax = static_cast<int>(tile.mesh_lengths[0]) + H;

    iterate_yee(yee, imin, imax, jmin, jmax, kmin, kmax, eval, mode);
  }

  //--------------------------------------------------
  if(left) {
    //std::cout << "left " << mins[0] << " " << mins[1] << " " << maxs[0] << " " << maxs[1] << std::endl;

    jmin = -H;
    jmax = static_cast<int>(tile.mesh_lengths[1]) + H;

    imin = -H;
    imax =  0;

    iterate_yee(yee, imin, imax, jmin, jmax, kmin, kmax, eval, mode);
  }

  //--------------------------------------------------
  if(top) {
    //std::cout << "top " << mins[0] << " " << mins[1] << " " << maxs[0] << " " << maxs[1] << std::endl;

    jmin = static_cast<int>(tile.mesh_lengths[1]);
    jmax = static_cast<int>(tile.mesh_lengths[1]) + H;

    imin = -H;
    imax = static_cast<int>(tile.mesh_lengths[0]) + H;

    iterate_yee(yee, imin, imax, jmin, jmax, kmin, kmax, eval, mode);
  }

  //--------------------------------------------------
  if(right) {
    //std::cout << "right " << mins[0] << " " << mins[1] << " " << maxs[0] << " " << maxs[1] << std::endl;

    jmin = -H;
    jmax = static_cast<int>(tile.mesh_lengths[1]) + H;

    imin = static_cast<int>(tile.mesh_lengths[0]);
    imax = static_cast<int>(tile.mesh_lengths[0]) + H;

    iterate_yee(yee, imin, imax, jmin, jmax, kmin, kmax, eval, mode);
  }

  
}


//--------------------------------------------------
// explicit template instantiation

template class fields::Conductor<2>; // 2D
template class fields::Conductor<3>; // 3D

