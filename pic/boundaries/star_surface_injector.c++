#include "star_surface_injector.h"
#include "../../em-fields/boundaries/conductor.h"

#include "../../tools/vector.h"

#include <cmath> 
#include <cassert>
#include <string>
#include <map>

using std::min;
using std::max;

// simple pseudo-random floats with C library rand() (outputting int's).
// note that we do not call srand( seed ) so it is set to seed(1). 
inline float rand_uni(float a, float b) {
  return ((b - a) * ((float)rand() / RAND_MAX)) + a;
}


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

  if( mins[1] < 1 )    bot   = true; 
  //if( mins[0] < 1 )    left  = true; 
  if( maxs[1] > Ny-1 ) top   = true; 
  //if( maxs[0] > Nx-1 ) right = true; 


  //--------------------------------------------------
  // operate only on roughly correct tiles 
  //float cenr = std::sqrt( cenx*cenx + ceny*ceny + cenz*cenz );
  //if(!(cenr + 5 <= maxs[1])) return;
  //if( maxs[1] - ceny > radius*1.02 ) return;
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

  // four-velocity components of the particles to be injected
  float_p ux1=0.0, uy1=0.0, uz1=0.0;

  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);

  auto& yee = tile.get_yee();

  const double c = tile.cfl;
  //const double q = container.q; // TODO needed

  // operate only on the bottom-y row tiles
  // TODO need to be generalized to arbitrary position

  // TODO add the container reference from QED codes

  //--------------------------------------------------
  // inject 

  const int k = 0; 
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    //-----------
    // global grid coordinates
    iglob = static_cast<float_m>(i) + mins[0];
    jglob = static_cast<float_m>(j) + mins[1];
    kglob = 0; //static_cast<float_m>(k) + mins[2];

    // spherical coordinates 
    xr0 = coord.rh().x(iglob);
    yr0 = coord.rh().y(jglob);
    zr0 = 0.0; 

    // check if we are inside star
    bool inside_star  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius + 0.0;
    bool inside_atmos = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius + 2.0;

    // approximate as flat surface
    //bool inside_star  = std::sqrt(yr0*yr0) <= 1.0*radius - 0.0;
    //bool inside_atmos = std::sqrt(yr0*yr0) <= 1.0*radius + 2.0;

    bool inside_pcap = std::abs(xr0) < 1.0*radius_pc;

    //if(xr0 > 0.0) continue;

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

      // injection rate
      // E is normalized with e n_GJ R_pc from which we can solve n_inj
        
      //ninj = 0.1*std::abs(epar/q/radius_pc);
      ninj = 0.05*std::abs(epar/q)/radius_pc;
      //std::cout << " ninj " << ninj << " epar" << epar << " epar/q" << epar/q << "\n";


      //ninj = std::max(0.002f, ninj);

      //--------------------------------------------------
      // ver2; current dependent inj
      //float_m jx = yee.jx(i,j,k);
      //float_m jy = yee.jy(i,j,k);
      //float_m jz = yee.jz(i,j,k);
      //float_m j = sqrt(jx*jx + jy*jy + jz*jz);
      
      // current is given as j = e*n_pcc*c so that n_ppc = j*c/e
      // We supply M_atms = 10x atmospheric particles required to screen the local current
      // n_atms = M_atms c (j/e)
      //ninj = 10.0*0.45*abs(j/q); 

      //std::cout << " ninj " << ninj << " j" << j << " j/q" << j/q << "\n";
      
      // TODO no need to smooth since epar is set zero outside polarcap
      //ninj *= shape( abs(xr0), radius_pc, delta_pc); // damp injection smoothly to zero outside polar cap
        
      ninj = std::max(0.01f, ninj);

      //--------------------------------------------------
      // add ninj pairs with MC injection; results on average in ninj injections
      double ncop = 0.0; // number of pairs added
      double z1 = rand_uni(0.0, 1.0);

      while( ninj > z1 + ncop ) {

        float dx = rand_uni(0.0, 1.0); // inject location is set randomly inside the cell

        //--------------------------------------------------
        // sample velocity from thermal distribution
        // using Box-Muller method to draw thermal velocities; valid for v <~ 0.2c
        double vth = 0.2;
        double rr1 = rand_uni(0.0, 1.0);
        double vr = sqrt( -2.0*log(rr1))*vth;
        
        // 1D distribution along B-field
        ux1 = vr*bx/b;
        uy1 = vr*by/b;
        uz1 = vr*bz/b;

        // TODO same "random" velocity taken for both particles

        cons["e-"]->add_particle( 
            {{iglob + dx, jglob , kglob}}, 
            {{ux1, uy1, uz1}}, wep); 

        cons["e+"]->add_particle( 
            {{iglob + dx, jglob, kglob}}, 
            {{ux1, uy1, uz1}}, wep); 

        ncop += 1.0;
      }

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

      iglob  = container.loc(0,n) + mins[0];
      jglob  = container.loc(1,n) + mins[1];
      kglob  = container.loc(2,n); // + mins[2];

      xr0 = coord.rh().x(iglob);
      yr0 = coord.rh().y(jglob);
      zr0 = 0.0; 

      // remove particles exiting from bottom
      bool inside_star = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius;
      if( //inside_star || 
          jglob < 3   || 
          inside_star && (std::abs(xr0) > radius_pc) ) {
        container.to_other_tiles.push_back( {1,1,1,n} );
      }

      // remove particles from the very top
      if( top && jglob > maxs[1]-1 ) {
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
