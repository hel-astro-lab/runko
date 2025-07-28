#include <cmath> 
#include <cassert>
#include <string>
#include <map>

#include "core/pic/boundaries/star_surface_injector.h"
#include "core/emf/boundaries/conductor.h"
#include "tools/vector.h"
#include "tools/signum.h"
#include "tools/staggered_grid.h"


using std::min;
using std::max;
using std::abs;
using std::sqrt;
using std::sin;
using std::cos;

using toolbox::sign;
using toolbox::norm;
using toolbox::norm1d;
using toolbox::norm2d;
using toolbox::cross;
using toolbox::dot;
using toolbox::Vec3;

using toolbox::StaggeredSphericalCoordinates;


// simple pseudo-random floats with C library rand() (outputting int's).
// note that we do not call srand( seed ) so it is set to seed(1). 
//
// NOTE: we now use the full Mersenne Twister from std
//
//inline float rand_uni(float a, float b) {
//  return ((b - a) * ((float)rand() / RAND_MAX)) + a;
//}


template<size_t D>
void pic::Star<D>::solve(
    pic::Tile<D>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;
    
  //bool left  = false;
  //bool right = false;
  bool top   = false;
  bool bot   = false;

  const int H = 2; // halo size

  if( D == 1 ) { // 1D boudaries
    if( mins[0] < 1.1*radius + cenx ) bot   = true; 
    if( maxs[0] > Nx-1 )              top   = true; 
  } else if( D == 2 ) { // 2D BCs
    if( mins[1] < 1.1*radius + ceny ) bot   = true; 
    if( maxs[1] > Ny-1 )              top   = true; 
  } else if( D == 3 ){ // 3D BCs
    if( mins[2] < 1.1*radius + cenz ) bot   = true; 
    if( maxs[2] > Nz-1 )              top   = true; 
  }


  //std::cout << " insider check: " 
  //          << " mins:" << mins[0] 
  //          << " maxs:" << maxs[0] 
  //          << " r1.1:" << 1.1*radius 
  //          << " cenx:" << cenx 
  //          << " val:" << 1.1*radius + cenx 
  //          << " bot " << bot 
  //          << " top " << top 
  //          << "\n";

  //--------------------------------------------------
  // operate only on roughly correct tiles 
  if(!(top || bot)) return;


  std::map<std::string, ConPtr> cons;
  for(auto&& con : tile.containers) cons.emplace(con.type, &con );

  // get charge (assume q_- = q_+)
  const float q = cons["e-"]->q;
  //float wep = 1.0f;
  //float wph = 1.0f;

  //--------------------------------------------------
  // main rouutine starts here now that we have the right tile
    
  Vec3<float> Om; 

  float Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check

  if(D == 1) Om.set(Omega,              0.0,   0.0); // Omega unit vector along x-axis
  if(D == 2) Om.set(0.0,                Omega, 0.0); // Omega unit vector along y-axis
  //if(D == 3) Om.set( sin(chi_om)*Omega, 0.0,   cos(chi_om)*Omega ); // general Omega unit vector 
  if(D == 3) Om.set( sin(chi_om)*cos(phase_om)*Omega, sin(chi_om)*sin(phase_om)*Omega, cos(chi_om)*Omega ); 

  //--------------------------------------------------
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  auto& gs = tile.get_grids();
  //const auto c = tile.cfl;

  //--------------------------------------------------
  // inject to each cell

  int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;


  //UniIter::iterate3D(
  //[=] DEVCALLABLE (int i, int j, int k, 
  //                Grids &gs)
  //{
  //}, 
  //  nx_tile, ny_tile, nz_tile,
  //  gs);

  for(int k=0; k<nz_tile; k++) 
  for(int j=0; j<ny_tile; j++) 
  for(int i=0; i<nx_tile; i++) {

    // global grid coordinates
    float iglob = (D>=1) ? i + mins[0] : 0.0;
    float jglob = (D>=2) ? j + mins[1] : 0.0;
    float kglob = (D>=3) ? k + mins[2] : 0.0;

    // spherical coordinates 
    auto rvec = coord.rh().vec(iglob, jglob, kglob, D); //cartesian radius vector in star's coords

    //--------------------------------------------------
    // check if we are inside star

    // inject top of star
    //const int height_atms = 1; // height of the atmosphere in cells
    //bool inside_star  = norm(rvec) < 1.0*radius;
    //bool inside_atmos = norm(rvec) < 1.0*radius + height_atms;

    //--------------------------------------------------
    // flat surface; slab on the top
    //bool inside_star  = (D==2) ? abs(rvec(1)) <= 1.0*radius + 0.0         : abs(rvec(2)) <= 1.0*radius + 0.0;
    //bool inside_atmos = (D==2) ? abs(rvec(1)) <= 1.0*radius + height_atms : abs(rvec(2)) <= 1.0*radius + height_atms;

    // flat surface; slab higher up
    //bool inside_star  = (D==2) ? abs(rvec(1)) <= 1.0*radius + 1.0*height_atms : abs(rvec(2)) <= 1.0*radius + 1.0*height_atms;
    //bool inside_atmos = (D==2) ? abs(rvec(1)) <= 1.0*radius + 2.0*height_atms : abs(rvec(2)) <= 1.0*radius + 2.0*height_atms;

    bool inside_star  = false;
    bool inside_atmos = false;

    // flat surface; slab higher up
    // NOTE: we have +1 height so that the update_e does not operate in the region we inject particles into
    if(D == 1){
      inside_star  = abs(rvec(0)) <= 1.0*radius + 1.0 + 1.0*height_atms;
      inside_atmos = abs(rvec(0)) <= 1.0*radius + 1.0 + 2.0*height_atms;
    } else if (D == 2){
      inside_star  = abs(rvec(1)) <= 1.0*radius + 1.0 + 1.0*height_atms;
      inside_atmos = abs(rvec(1)) <= 1.0*radius + 1.0 + 2.0*height_atms;
    } else if (D == 3){
      inside_star  = abs(rvec(2)) <= 1.0*radius + 1.0 + 1.0*height_atms;
      inside_atmos = abs(rvec(2)) <= 1.0*radius + 1.0 + 2.0*height_atms;
    }


    //--------------------------------------------------
    // inject below surface
    //bool inside_star  = norm(rvec) <= 1.0*radius - 3.0;
    //bool inside_atmos = norm(rvec) <= 1.0*radius + 0.0;


    //--------------------------------------------------
    // inject exactly the given cell thickness

    //int th = 0;
    //if(!inside_star && inside_atmos) {

    //  // axis of depth
    //  int dir = 0;
    //  if(      abs(rvec(0)) > max(abs(rvec(1)), abs(rvec(2))) ) dir = 0; // x dir
    //  else if( abs(rvec(1)) > max(abs(rvec(0)), abs(rvec(2))) ) dir = 1; // y dir
    //  else if( abs(rvec(2)) > max(abs(rvec(0)), abs(rvec(1))) ) dir = 2; // z dir
    //                                                                
    //  // drill below and check how many cells we can go deeper to be inside star
    //  auto rtmp = rvec;
    //  for(int h=0; h<height_atms; h++) {
    //    rtmp(dir) += -1.0f*sign(rtmp(dir)); // go one below to the direction of negative r
    //    if( norm(rtmp) < radius) th += 1;
    //  }

    //  // check if thickness is the wanted value and then inject
    //  inside_atmos = false;
    //  if( th-1 < height_atms ) inside_atmos = true;
    //}


    //--------------------------------------------------
    const float offs = 2.0f*delta_pc; // expand polar cap a bit
    //bool inside_pcap = (D==2) ? norm1d(rvec) < radius_pc + offs : norm2d(rvec) < radius_pc + offs;

    bool inside_pcap = false;
    if(D == 1){
      inside_pcap = true; // always inside pcap in 1D
    } else if (D == 2){
      inside_pcap = norm1d(rvec) < radius_pc + offs;
    } else if (D == 3){
      inside_pcap = norm2d(rvec) < radius_pc + offs;
    }
  

    //std::cout << "ijk " << i << " " << j << " " << k 
    //          << " conds:" << inside_atmos 
    //          << " " << inside_star 
    //          << " " << inside_pcap 
    //          << " full:" << (inside_atmos && !inside_star && inside_pcap)
    //          << "\n";

    //--------------------------------------------------
    // we are inside a thin layer above the star
    if( inside_atmos && !inside_star && inside_pcap ) {

      // debug to show the injection region
      //gs.ex(i,j,k) = 10.0;
      //gs.ey(i,j,k) = 10.0;
      //gs.ez(i,j,k) = 10.0;

      // read current from this many cells ahead of the atmopshere
      // this leap is needed so that electric currents have effect on the E field and we measure the right
      // electric field. 
      // NOTE: the value here should be about 2x delta ~ 2 cells
      int isk = D==1 ? 3 : 0;
      int jsk = D==2 ? 3 : 0;
      int ksk = D==3 ? 3 : 0;


      // get epar (TODO not on the right staggering)
      auto ex = gs.ex(i+isk,j+jsk,k+ksk);
      auto ey = gs.ey(i+isk,j+jsk,k+ksk);
      auto ez = gs.ez(i+isk,j+jsk,k+ksk);

      auto bx = gs.bx(i+isk,j+jsk,k+ksk);
      auto by = gs.by(i+isk,j+jsk,k+ksk);
      auto bz = gs.bz(i+isk,j+jsk,k+ksk);

      auto b    = sqrt( bx*bx + by*by + bz*bz ) + EPS;
      auto epar = ( ex*bx + ey*by + ez*bz )/b;

      // vectors for calculation of pcap rotation velocity
      Vec3 E(ex, ey, ez);
      Vec3 B(bx, by, bz);

      //--------------------------------------------------
      // ver 1: 
      // E is normalized with e n_GJ R_pc from which we can solve n_inj
        
      float ninj = ninj_pairs*abs(epar/q)/radius_pc;
      //std::cout << " ninj " << ninj << " epar" << epar << " epar/q" << epar/q << "\n";


      //--------------------------------------------------
      // ver2; current dependent inj; NOTE does not work because current arrays are emptied
      //float jx = gs.jx(i,j,k);
      //float jy = gs.jy(i,j,k);
      //float jz = gs.jz(i,j,k);
      //float j = sqrt(jx*jx + jy*jy + jz*jz);
      
      // current is given as j = e*n_pcc*c so that n_ppc = j*c/e
      // We supply M_atms = 10x atmospheric particles required to screen the local current
      // n_atms = M_atms c (j/e)
      //ninj = 10.0*0.45*abs(j/q); 

      //std::cout << " ninj " << ninj << " j" << j << " j/q" << j/q << "\n";
      ninj = max( (float)ninj_min_pairs, ninj);

      //--------------------------------------------------
      // add ninj pairs with MC injection; results on average in ninj injections
      float ncop = 0.0f; // number of pairs added
      float z1 = rand();

      auto n_to_be_inj = static_cast<size_t>(max(1.0f, std::ceil(ninj)));

      // pre-created arrays for particles
      auto x_to_be_inj  = std::vector<float>(n_to_be_inj);
      auto y_to_be_inj  = std::vector<float>(n_to_be_inj);
      auto z_to_be_inj  = std::vector<float>(n_to_be_inj);
      auto ux_to_be_inj = std::vector<float>(n_to_be_inj);
      auto uy_to_be_inj = std::vector<float>(n_to_be_inj);
      auto uz_to_be_inj = std::vector<float>(n_to_be_inj);

      //std::cout << " inj ncop:" << ninj << "\n";

      //--------------------------------------------------
      //add particles
      while( ninj > z1 + ncop ) {

        // inject location is set randomly inside the cell
        float dx = (D >= 1) ? rand() : 0.0f; 
        float dy = (D >= 2) ? rand() : 0.0f; 
        float dz = (D >= 3) ? rand() : 0.0f; 

        //--------------------------------------------------
        // sample velocity from thermal distribution
        // using Box-Muller method to draw thermal velocities; valid for v <~ 0.2c
        float rr1 = rand(); //zeta4[ncop];
        float vr = sqrt( -2.0f*log(rr1))*temp_pairs;
        
        // 1D distribution along B-field
        // NOTE: same "random" velocity taken for both particles as an approx
        //auto ux1 = vr*bx/b;
        //auto uy1 = vr*by/b;
        //auto uz1 = vr*bz/b;

        // Using now a random 3D distribution instead:
        float zeta = 2.0f*PI*rand();
        float mu = -1.0f + 2.0f*rand();
        float sin_theta = sqrt(1.0f-pow(mu,2));

        auto ux1 = vr*sin_theta*cos(zeta);
        auto uy1 = vr*sin_theta*sin(zeta);
        auto uz1 = vr*mu;

        //--------------------------------------------------
        // pcap rotation vector
        //auto r = rvec;  // copy
        //r(0) += dx;
        //r(1) += dy;
        //r(2) += dz;
        //auto vrot = cross(Om, r); 

        // ExB version
        //auto B2E2 = 1.0f/( dot(B, B) ); //+ dot(Ev, Ev) );
        //auto vrot2 = B2E2*cross(E, B);
        // NOTE matches to about ~30% the surface rotation velocity
          
        //std::cout << "xr0  :   " << xr0   << "\n";
        //std::cout << "u    : " << ux1 << " " << uy1 << " " << uz1  << "\n";
        //std::cout << "vrot1: " << vrot  << "\n";
        //std::cout << "vrot2: " << vrot2 << "\n";
        //std::cout << "  rat: " << vrot(0)/vrot2(0) << " " << vrot(1)/vrot2(1) << "\n\n";

        // Newtonian boost to the disk frame; should be done with Lorentz boost to be more correct
        //ux1 += vrot2(0);
        //uy1 += vrot2(1);
        //uz1 += vrot2(2);

        x_to_be_inj[ncop]  = iglob + dx;
        y_to_be_inj[ncop]  = jglob + dy;
        z_to_be_inj[ncop]  = kglob + dz;
        ux_to_be_inj[ncop] = ux1;
        uy_to_be_inj[ncop] = uy1;
        uz_to_be_inj[ncop] = uz1;

        //cons["e-"]->add_particle( 
        //    {{iglob + dx, jglob + dy, kglob + dz}}, 
        //    {{ux1, uy1, uz1}}, wep); 

        //cons["e+"]->add_particle( 
        //    {{iglob + dx, jglob + dy, kglob + dz}}, 
        //    {{ux1, uy1, uz1}}, wep); 

        ncop += 1;
      }


      //--------------------------------------------------
      // add the pre-created particles
      for(int n=0; n<ncop; n++){
        //std::cout <<" x y z " << 
        //  x_to_be_inj[n] << " " <<
        //  y_to_be_inj[n] << " " <<
        //  z_to_be_inj[n] << " " 
        //  << " ux uy uz " <<
        //  ux_to_be_inj[n] << " " <<
        //  uy_to_be_inj[n] << " " <<
        //  uz_to_be_inj[n] << "\n";
        cons["e-"]->add_particle( {{  x_to_be_inj[n],  y_to_be_inj[n],  z_to_be_inj[n] }}, 
                                  {{ ux_to_be_inj[n], uy_to_be_inj[n], uz_to_be_inj[n] }}, wep); 
      }

      for(int n=0; n<ncop; n++){
        cons["e+"]->add_particle( {{  x_to_be_inj[n],  y_to_be_inj[n],  z_to_be_inj[n] }}, 
                                  {{ ux_to_be_inj[n], uy_to_be_inj[n], uz_to_be_inj[n] }}, wep); 
      }


      //--------------------------------------------------
      // photon injection

      ninj = ninj_phots; //*abs(epar/q)/radius_pc;
      ninj = max( (float)ninj_min_phots, ninj);

      // maximum amount to be injected
      n_to_be_inj = static_cast<size_t>(max(1.0f, std::ceil(ninj)));

      // pre-created arrays for particles; not used
      x_to_be_inj.resize( n_to_be_inj);
      y_to_be_inj.resize( n_to_be_inj);
      z_to_be_inj.resize( n_to_be_inj);
      ux_to_be_inj.resize(n_to_be_inj);
      uy_to_be_inj.resize(n_to_be_inj);
      uz_to_be_inj.resize(n_to_be_inj);


      ncop = 0.0f; // number of phots added
      z1 = rand();

      //--------------------------------------------------
      while( ninj > z1 + ncop ) {

        //--------------------------------------------------
        // draw random isotropic 3d vector
        //float vz = 2.0f*rand() - 1.0f; // TODO maybe only upper half of the sphere is more realistic?
        //float xi0 = 2.0f*PI*rand();
        //float vx = sqrt(1.0f-vz*vz)*cos(xi0);
        //float vy = sqrt(1.0f-vz*vz)*sin(xi0);

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

        //auto vl = sqrt( vx*vx + vy*vy + vz*vz );
        //std::cout << " ph : " << vx << " " << vy << " " << vz << " len: " << vl << "\n";

        //cons["ph"]->add_particle( {{iglob, jglob, kglob}}, {{ux, uy, uz}}, wph );

        x_to_be_inj[ncop]  = iglob;
        y_to_be_inj[ncop]  = jglob;
        z_to_be_inj[ncop]  = kglob;
        ux_to_be_inj[ncop] = ux;
        uy_to_be_inj[ncop] = uy;
        uz_to_be_inj[ncop] = uz;

        ncop += 1;
      }

      //--------------------------------------------------
      // add the pre-created particles
      for(int n=0; n<ncop; n++){
        cons["ph"]->add_particle( {{  x_to_be_inj[n],  y_to_be_inj[n],  z_to_be_inj[n] }}, 
                                  {{ ux_to_be_inj[n], uy_to_be_inj[n], uz_to_be_inj[n] }}, wph); 
      }



    }
  }

  //--------------------------------------------------
  // remove outflowing particles

  //float tile_len = (D == 2 ) ? tile.mesh_lengths[1] : tile.mesh_lengths[2];
  //float rbox = (D == 2) ? 0.5*Nx - 0.5*tile_len : 0.5*Nx - H - 1; // maximum cylindrical radius
  //float tile_height = (D==2) ? tile.mesh_lengths[1] : tile.mesh_lengths[2]; // height of the tile

  float rbox      = 0.0;
  if(D == 1) rbox = 0.0; // N/A
  if(D == 2) rbox = 0.5*Nx - 0.5*tile.mesh_lengths[1]; // maximum cylindrical radius
  if(D == 3) rbox = 0.5*Nx - H - 1;                    // maximum cylindrical radius


  // call pre-iteration functions to update internal arrays 
  //for(auto&& con : tile.containers) {
  //    con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
  //}


  for(auto&& con : tile.containers) {
    for(size_t n=0; n<con.size(); n++) {

      float iglob = (D>=1) ? con.loc(0,n) : 0.0;
      float jglob = (D>=2) ? con.loc(1,n) : 0.0;
      float kglob = (D>=3) ? con.loc(2,n) : 0.0;

      float xr0 = (D>=1) ? coord.mid().x(iglob) : 0.0;
      float yr0 = (D>=2) ? coord.mid().y(jglob) : 0.0;
      float zr0 = (D>=3) ? coord.mid().z(kglob) : 0.0;

      // remove particles based on following regimes
      //bool inside_star    = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius;
      //bool below_star     = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius-2; // NOTE hard-coded thickness of 2

      // flat surface
      //bool inside_star    = (D == 2) ? abs(yr0) < 1.0*radius : abs(zr0) < 1.0*radius; 

      //bool below_star     = (D == 2) ? abs(yr0) < 1.0*radius : abs(zr0) < 1.0*radius - 1; 
      //bool inside_bot     = (D == 2) ? jglob < H    : kglob < H; // y or z direction flip 
      //bool inside_top     = (D == 2) ? jglob > Ny - 0.75*tile_height : kglob > Nz - 0.75*tile_height; // y or z direction flip 
      //bool outside_pcap   = (D == 2) ? abs(xr0) > radius_pc : sqrt( xr0*xr0 + yr0*yr0 ) > radius_pc; // same here
      //bool inside_cyl_bcs = (D == 2) ? abs(xr0) > rbox : sqrt(xr0*xr0 + yr0*yr0) > rbox; // box sides covering a cylindrical region


      bool below_star = false;
      if(D == 1) below_star = abs(xr0) < 1.0*radius;
      if(D == 2) below_star = abs(yr0) < 1.0*radius;
      if(D == 3) below_star = abs(zr0) < 1.0*radius;

      bool inside_bot = false;
      if(D == 1) inside_bot = iglob < H; // x direction 
      if(D == 2) inside_bot = jglob < H; // y direction 
      if(D == 3) inside_bot = kglob < H; // z direction 

      bool inside_top = false;
      if(D == 1) inside_top = iglob > Nx - 0.75*tile.mesh_lengths[0];
      if(D == 2) inside_top = jglob > Ny - 0.75*tile.mesh_lengths[1];
      if(D == 3) inside_top = kglob > Nz - 0.75*tile.mesh_lengths[2]; 


      // additional, more aggressive removal of outflowing particles in 1D
      if(D == 1) {
	      if( iglob > Nx - 20.0*tile.mesh_lengths[0] ) { 
	        // negative (inwards-going) particle
	        if( con.vel(0, n) < 0.0f ) con.info(n) = -1;
	      }
      }

      bool inside_cyl_bcs = false;
      if(D == 1) inside_cyl_bcs = false; // no boundaries in 1D
      if(D == 2) inside_cyl_bcs = abs(xr0) > rbox; // box sides covering a cylindrical region
      if(D == 3) inside_cyl_bcs = sqrt(xr0*xr0 + yr0*yr0) > rbox; // box sides covering a cylindrical region

      if( inside_bot                  ||
          //inside_star && outside_pcap ||
          below_star                  ||
          inside_cyl_bcs              || 
          (top && inside_top)              // particles at the very top
          ) {
        //con.to_other_tiles.push_back( {1,1,1,n} );
        con.info(n) = -1;
      }
    }
  }

  // remove annihilated prtcls; this transfer storage 
  // is used as a tmp container for storing the indices
  for(auto&& con : tile.containers) {
    con.delete_transferred_particles(); 
  }


  return;
}


template class pic::Star<1>; // 1D
template class pic::Star<2>; // 2D
template class pic::Star<3>; // 3D
