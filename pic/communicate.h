#pragma once

#include <algorithm>

#include "cell.h"

using std::min;
using std::max;


namespace pic {

class Communicator {
  public:


  // expose particle memory of external neighbors
  pic::ParticleBlock& get_external_data(
    int i, int j,
    pic::PicCell& cell,
    corgi::Node& grid)
  { 
    auto ind = cell.neighs(i, j); 
    uint64_t cid = grid.cellId( std::get<0>(ind), std::get<1>(ind) );
    pic::PicCell& external_cell = dynamic_cast<pic::PicCell&>( grid.getCell(cid) );

    return external_cell.container;
  }


  void transfer( pic::PicCell& cell, corgi::Node& grid)
  {

    // initialize neighbors
    pic::ParticleBlock& left  = get_external_data(-1,0, cell, grid);
    pic::ParticleBlock& right = get_external_data(+1,0, cell, grid);


    // initialize pointers to particle arrays
    int nparts = cell.container.size();

    // block limits
    double xmin = 0.0;
    double xmax = 1.0*cell.container.Nx;

    double ymin = 0.0;
    double ymax = 1.0*cell.container.Ny;

    double zmin = 0.0;
    double zmax = 1.0*cell.container.Nz;


    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double x0, y0, z0;
    double dx, dy, dz;

    // loop and check particles
    int n1 = 0;
    int n2 = nparts;

    for(int n=n1; n<n2; n++) {

      x0 = loc[0][n];
      y0 = loc[1][n];
      z0 = loc[2][n];

      // left wrap
      if(x0 < xmin) {
        dx = abs(xmin - x0);
        loc[0][n] = xmax - dx;
      }

      //right wrap
      if(x0 > xmin) {
        dx = abs(x0 - xmax);
        loc[0][n] = xmin + dx;
      }


    }

  }



};


} // end of namespace pic
