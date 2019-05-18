#include "piston.h"
#include "../../tools/signum.h"

#include <cmath> 
#include <cassert>

template<size_t D>
void pic::Piston<D>::zigzag(
    pic::Tile<D>& tile,
    double x0, 
    double y0, 
    double z0, 
    double x1, 
    double y1, 
    double z1, 
    double q)
{
  auto& yee = tile.get_yee();


}

template<size_t D>
void pic::Piston<D>::solve(
    pic::Tile<D>& tile)
{

  // get reference to the Yee grid 
  //auto& yee = tile.get_yee();

  for(auto&& container : tile.containers) {
    int nparts = container.size();

    // initialize pointers to particle arrays
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );

    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    double c = tile.cfl;
    //double qm = sign(container.q);

    double x0, y0, z0, gamma, tfrac, walloc0;
    double xcolis, ycolis, zcolis;
    for(int n=n1; n<n2; n++) {

      // left side of the wall boundary
      if(loc[0][n] < walloc) {

        gamma = sqrt(1.0
            + vel[0][n]*vel[0][n]
            + vel[1][n]*vel[1][n]
            + vel[2][n]*vel[2][n]);

        x0 = loc[0][n] - vel[0][n]/gamma*c;
        y0 = loc[1][n];
        z0 = loc[2][n];

        // unwind wall location
        walloc0 = walloc - betawall*c;

        // compute crossing point
        tfrac = abs((x0-walloc0)/(betawall*c - vel[0][n]/gamma*c));
        xcolis = x0 + vel[0][n]/gamma*c*tfrac;
        ycolis = y0;
        zcolis = z0;

        // deposit current up to intersection point
        zigzag(tile, xcolis, ycolis, zcolis, x0, y0, z0, container.q);

        // reset particle momentum, getting kick from the wall
        vel[0][n] = gammawall*gammawall*gamma
          *(2.*betawall - vel[0][n]/gamma*(1.0+betawall*betawall));

        // move particle from the location of the wall with new velocity
        loc[0][n] = xcolis + vel[0][n]/gamma*c * tfrac;
        loc[1][n] = ycolis;
        loc[2][n] = zcolis;

        // clean up the part of trajectory behind the wall
        // that will be added by the deposition routine that unwinds
        // the particle location by a full time step.
        zigzag(tile, xcolis, ycolis, zcolis,
            loc[0][n] - vel[0][n]/gamma*c,
            loc[1][n] - vel[1][n]/gamma*c,
            loc[2][n] - vel[2][n]/gamma*c,
            -container.q);

      }
    }

  } // end of loop over species

  return;
}


//--------------------------------------------------
// explicit template instantiation

//template class pic::Piston<1>; // 1D3V
template class pic::Piston<2>; // 2D3V
//template class pic::Piston<3>; // 3D3V

