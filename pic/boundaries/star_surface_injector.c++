#include "star_surface_injector.h"
#include "../../em-fields/boundaries/conductor.h"

#include "../../tools/vector.h"

#include <cmath> 
#include <cassert>
#include <string>
#include <map>

using std::min;
using std::max;
using std::abs;
using std::sqrt;

using toolbox::cross;
using toolbox::dot;
using toolbox::Vec3;


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

  if( D == 2 ) { // 2D boudaries
    if( mins[1] < 1 )    bot   = true; 
    if( maxs[1] > Ny-1 ) top   = true; 
  } else if( D == 3 ){
    if( mins[2] < 1 )    bot   = true; 
    if( maxs[2] > Nz-1 ) top   = true; 
  }

  //--------------------------------------------------
  // operate only on roughly correct tiles 
  if(!(top || bot)) return;


  std::map<std::string, ConPtr> cons;
  for(auto&& con : tile.containers) cons.emplace(con.type, &con );

  // get charge (assume q_- = q_+)
  const float_m q = cons["e-"]->q;
  float_p wep = 1.0f;

  //--------------------------------------------------
  // main rouutine starts here now that we have the right tile
    
  float_m iglob, jglob, kglob;
  float_m xr0,yr0,zr0;
  float_m ex, ey, ez, bx, by, bz;
  float_m epar, ninj;

  Vec3<float_m> Bv, Ev; // tmp variables
  Vec3<float_m> r, Om; // tmp variables

  float_m Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check

  if(D == 2) Om.set(0.0,                Omega, 0.0); // Omega unit vector along y-axis
  if(D == 3) Om.set( sin(chi_om)*Omega, 0.0,   cos(chi_om)*Omega ); // general Omega unit vector 

  //--------------------------------------------------
  // four-velocity components of the particles to be injected
  float_p ux1=0.0, uy1=0.0, uz1=0.0;

  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  auto& yee = tile.get_yee();
  const float_p c = tile.cfl;


  //--------------------------------------------------
  // inject to each cell

  int nx_tile = (D>=1) ? tile.mesh_lengths[0] : 1;
  int ny_tile = (D>=2) ? tile.mesh_lengths[1] : 1;
  int nz_tile = (D>=3) ? tile.mesh_lengths[2] : 1;

  for(int k=0; k<nz_tile; k++) 
  for(int j=0; j<ny_tile; j++) 
  for(int i=0; i<nx_tile; i++) {

    // global grid coordinates
    iglob = (D>=1) ? i + mins[0] : 0.0;
    jglob = (D>=2) ? j + mins[1] : 0.0;
    kglob = (D>=3) ? k + mins[2] : 0.0;

    // spherical coordinates 
    xr0 = (D>=1) ? coord.rh().x(iglob) : 0.0;
    yr0 = (D>=2) ? coord.rh().y(jglob) : 0.0;
    zr0 = (D>=3) ? coord.rh().z(kglob) : 0.0;

    // check if we are inside star
    bool inside_star  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius - 3.0;
    bool inside_atmos = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius + 0.0;
    bool inside_pcap = (D == 2) ? abs(xr0) < radius_pc : sqrt( xr0*xr0 + yr0*yr0 ) < radius_pc; // same here

    // inside a thin layer above the star
    if( inside_atmos && !inside_star && inside_pcap ) {

      // get epar (TODO not on the right staggering)
      ex = yee.ex(i,j,k);
      ey = yee.ey(i,j,k);
      ez = yee.ez(i,j,k);

      bx = yee.bx(i,j,k);
      by = yee.by(i,j,k);
      bz = yee.bz(i,j,k);

      float_m b = sqrt( bx*bx + by*by + bz*bz );
      epar      = ( ex*bx + ey*by + ez*bz )/b;


      // vectors for calculation of pcap rotation velocity
      Ev.set(ex, ey, ez);
      Bv.set(bx, by, bz);

      //--------------------------------------------------
      // ver 1: 
      // E is normalized with e n_GJ R_pc from which we can solve n_inj
        
      //ninj = 0.1*std::abs(epar/q/radius_pc);
      ninj = ninj_pairs*std::abs(epar/q)/radius_pc;
      //std::cout << " ninj " << ninj << " epar" << epar << " epar/q" << epar/q << "\n";

      //--------------------------------------------------
      // ver2; current dependent inj; NOTE does not work because current arrays are emptied
      //float_m jx = yee.jx(i,j,k);
      //float_m jy = yee.jy(i,j,k);
      //float_m jz = yee.jz(i,j,k);
      //float_m j = sqrt(jx*jx + jy*jy + jz*jz);
      
      // current is given as j = e*n_pcc*c so that n_ppc = j*c/e
      // We supply M_atms = 10x atmospheric particles required to screen the local current
      // n_atms = M_atms c (j/e)
      //ninj = 10.0*0.45*abs(j/q); 

      //std::cout << " ninj " << ninj << " j" << j << " j/q" << j/q << "\n";
        
      ninj = std::max( (float_m)ninj_min_pairs, ninj);

      //--------------------------------------------------
      // add ninj pairs with MC injection; results on average in ninj injections
      double ncop = 0.0; // number of pairs added
      double z1 = rand();

      while( ninj > z1 + ncop ) {

        float dx = rand(); // inject location is set randomly inside the cell
        float dy = rand(); // inject location is set randomly inside the cell

        //--------------------------------------------------
        // sample velocity from thermal distribution
        // using Box-Muller method to draw thermal velocities; valid for v <~ 0.2c
        double rr1 = rand();
        double vr = sqrt( -2.0*log(rr1))*temp_pairs;
        
        // 1D distribution along B-field
        ux1 = vr*bx/b;
        uy1 = vr*by/b;
        uz1 = vr*bz/b;

        //--------------------------------------------------
        // pcap rotation vector
        r.set(xr0 + dx, yr0 + dy, zr0);
        auto vrot = cross(Om, r); 

        // ExB version
        auto B2E2 = 1.0f/( dot(Bv, Bv) ); //+ dot(Ev, Ev) );
        auto vrot2 = B2E2*cross(Ev, Bv);
        // NOTE amtches to about ~30% the surface rotation velocity

        //std::cout << "vrot1: " << vrot  << "\n";
        //std::cout << "vrot2: " << vrot2 << "\n";
        //std::cout << "  rat: " << vrot(0)/vrot2(0) << " " << vrot(1)/vrot2(1) << "\n\n";

        // Newtonian boost to the disk frame; should be done with Lorentz boost to be more correct
        ux1 += vrot2(0);
        uy1 += vrot2(1);
        uz1 += vrot2(2);

        // TODO same "random" velocity taken for both particles

        cons["e-"]->add_particle( 
            {{iglob + dx, jglob + dy, kglob}}, 
            {{ux1, uy1, uz1}}, wep); 

        cons["e+"]->add_particle( 
            {{iglob + dx, jglob + dy, kglob}}, 
            {{ux1, uy1, uz1}}, wep); 

        ncop += 1.0;
      }

      //--------------------------------------------------
      // TODO add photon injection





    }
  }

  //--------------------------------------------------
  // remove outflowing 
    
  // call pre-iteration functions to update internal arrays 
  for(auto&& con : tile.containers) {
      con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
  }

  for(auto&& container : tile.containers) {
    for(size_t n=0; n<container.size(); n++) {

      iglob = (D>=1) ? container.loc(0,n) : 0.0;
      jglob = (D>=2) ? container.loc(1,n) : 0.0;
      kglob = (D>=3) ? container.loc(2,n) : 0.0;

      xr0 = (D>=1) ? coord.rh().x(iglob) : 0.0;
      yr0 = (D>=2) ? coord.rh().y(jglob) : 0.0;
      zr0 = (D>=3) ? coord.rh().z(kglob) : 0.0;

      // remove particles based on following regimes
      bool inside_star  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius;
      bool inside_bot   = (D == 2) ? jglob < 3    : kglob < 3; // y or z direction flip 
      bool inside_top   = (D == 2) ? jglob > Ny-1 : kglob > Nz-1; // y or z direction flip 
      bool outside_pcap = (D == 2) ? abs(xr0) > radius_pc : sqrt( xr0*xr0 + yr0*yr0 ) > radius_pc; // same here

      if( inside_bot ||
          inside_star && outside_pcap ) {
        container.to_other_tiles.push_back( {1,1,1,n} );
      }

      // remove particles from the very top
      if( top && inside_top ) {
        container.to_other_tiles.push_back( {1,1,1,n} );
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


template class pic::Star<2>; // 2D
template class pic::Star<3>; // 3D
