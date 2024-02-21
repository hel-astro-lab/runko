#include "star_surface_injector.h"
#include "../../em-fields/boundaries/conductor.h"

#include "../../tools/vector.h"

#include <cmath> 
#include <cassert>
#include <string>
#include <map>

using std::min;
using std::max;


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
    xr0 = coord.x(iglob, 0.0);
    yr0 = coord.y(jglob, 0.0);
    zr0 = 0.0; //coord.z(kglob, 0.5);

    // check if we are inside star
    bool inside_star  = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius - 0.0;
    bool inside_atmos = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius + 2.0;

    // approximate as flat surface
    //bool inside_star  = std::sqrt(yr0*yr0) <= 1.0*radius - 0.0;
    //bool inside_atmos = std::sqrt(yr0*yr0) <= 1.0*radius + 2.0;

    // inside a thin layer above the star
    if( inside_atmos && !inside_star ) {

      // get epar
      ex = yee.ex(i,j,k);
      ey = yee.ey(i,j,k);
      ez = yee.ez(i,j,k);

      bx = yee.bx(i,j,k);
      by = yee.by(i,j,k);
      bz = yee.bz(i,j,k);

      epar = ( ex*bx + ey*by + ez*bz )/( bx*bx + by*by + bz*bz );

      // injection rate
      ninj = 0.2*std::abs(epar/q);
      //std::cout << " ninj " << ninj << " epar" << epar << " epar/q" << epar/q << "\n";

      //ninj = 10; // FIXME

      // TODO no need to smooth since epar is set zero outside polarcap
      //ninj *= shape( abs(xr0), radius_pc, delta_pc); // damp injection smoothly to zero outside polar cap
                                                       
      //--------------------------------------------------
      // add ninj pairs with MC injection; results on average in ninj injections
      double ncop = 0.0; // number of pairs added
      double z1 = rand();
      while( ninj > z1 + ncop ) {

        cons["e-"]->add_particle( 
            {{iglob, jglob, kglob}}, 
            {{ux1, uy1, uz1}}, wep); 

        cons["e+"]->add_particle( 
            {{iglob, jglob, kglob}}, 
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

      xr0 = coord.x(iglob, 0.0);
      yr0 = coord.y(jglob, 0.0);
      zr0 = 0.0; //coord.z(kglob, 0.5);

      // remove particles exiting from bottom
      bool inside_star = std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) <= 1.0*radius;
      if( //inside_star || 
          jglob < 3   || 
          inside_star && (abs(xr0) > radius_pc) ) {
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
