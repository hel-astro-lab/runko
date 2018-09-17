#include "linear.h"

#include <cmath> 
#include <cassert>

template<size_t D, size_t V>
void pic::LinearInterpolator<D,V>::solve(
    pic::Tile<D>& tile)
{

  // get reference to the Yee grid 
  auto& yee = tile.getYee();


  for (size_t ispc=0; ispc<tile.Nspecies(); ispc++) {
    ParticleBlock& container = tile.get_container(ispc);

    int nparts = container.size();
    container.resizeEM(nparts, 3); // make EM containers ready for insertion


    // initialize pointers to particle arrays
    Realf* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    // 1-d arrays
    //Realf* ex = &( (*tile.container.Epart)[0*nparts] );
    //Realf* ey = &( (*tile.container.Epart)[1*nparts] );
    //Realf* ez = &( (*tile.container.Epart)[2*nparts] );

    // multiD array version
    //Realf *efield[3], *bfield[3];
    //for( int i=0; i<3; i++) {
    //  efield[i] = &( tile.container.Epart[i][0] );
    //  bfield[i] = &( tile.container.Bpart[i][0] );
    //}
      
    Realf *ex, *ey, *ez, *bx, *by, *bz;
    ex = &( container.Epart[0][0] );
    ey = &( container.Epart[1][0] );
    ez = &( container.Epart[2][0] );

    bx = &( container.Bpart[0][0] );
    by = &( container.Bpart[1][0] );
    bz = &( container.Bpart[2][0] );


    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    //Realf c = tile.cfl;
    //Realf cinv = 1.0/c;

    int i=0, j=0, k=0;
    Realf dx=0.0, dy=0.0, dz=0.0;
    Realf f,g;

    int iz = 0; // flip switch for making array queries 2D

    auto mins = tile.mins;
    //auto maxs = tile.maxs;

    // TODO: think SIMD (not possible due to ijk writing to yee)
    for(int n=n1; n<n2; n++) {

      // particle location in the grid
        
      if (D >= 1) {
	      i  = floor( loc[0][n]-mins[0] );
	      dx = (loc[0][n]-mins[0]) - i;
      }

      if (D >= 2) {
	      j  = floor( loc[1][n]-mins[1] );
	      dy = (loc[1][n]-mins[1]) - j;
      }

      if (D >= 3) {
	      k  = floor( loc[2][n]-mins[2] );
	      dz = (loc[2][n]-mins[2]) - k;
      }

      /*
      std::cout << '\n';
      std::cout << "x: " << loc[0][n] << " y: " << loc[1][n] << " z:" << loc[2][n] << '\n';
      std::cout << " ijk " << i << "," << j << "," << k << '\n';
      std::cout << " ds " << dx << "," << dy << "," << dz << '\n';
      */
        
	    //l = i; // + iy*(j-1) + iz*(k-1);

      if (D >= 1) assert(i >= 0 && i < static_cast<int>(tile.mesh_lengths[0]) );
      if (D >= 2) assert(j >= 0 && j < static_cast<int>(tile.mesh_lengths[1]) );
      if (D >= 3) assert(k >= 0 && k < static_cast<int>(tile.mesh_lengths[2]) );


      // TODO: these can be optimized further when we know D
      f = yee.ex(i,j,k) + yee.ex(i-1,j,k) +    dx*(yee.ex(i+1,j  ,k   ) - yee.ex(i-1,j  ,k  ));
      f+=                                      dy*(yee.ex(i  ,j+1,k   ) + yee.ex(i-1,j+1,k  )+
                                               dx*(yee.ex(i+1,j+1,k   ) - yee.ex(i-1,j+1,k  ))-f);
      g = yee.ex(i,j,k+1)+yee.ex(i-1,j,k+iz)+  dx*(yee.ex(i+1,j  ,k+iz) - yee.ex(i-1,j  ,k+iz));
      g+=                                      dy*(yee.ex(i  ,j+1,k+iz) + yee.ex(i-1,j+1,k+iz)
                                           +   dx*(yee.ex(i+1,j+1,k+iz) - yee.ex(i-1,j+1,k+iz))-g);
      ex[n] = (f+dz*(g-f))*0.5;

      f = yee.ey(i,j,k)+yee.ey(i,j-1,k)+       dy*(yee.ey(i  ,j+1,k   ) - yee.ey(i  ,j-1,k  ));
      f+=                                      dz*(yee.ey(i  ,j  ,k+iz) + yee.ey(i  ,j-1,k+iz)+      
                                               dy*(yee.ey(i  ,j+1,k+iz) - yee.ey(i  ,j-1,k+iz))-f);
      g = yee.ey(i+1,j,k)+yee.ey(i+1,j-1,k)+   dy*(yee.ey(i+1,j+1,k   ) - yee.ey(i+1,j-1,k   ));
      g+=                                      dz*(yee.ey(i+1,j  ,k+iz) + yee.ey(i+1,j-1,k+iz)+ 
                                               dy*(yee.ey(i+1,j+1,k+iz) - yee.ey(i+1,j-1,k+iz))-g);
      ey[n]=(f+dx*(g-f))*0.5;

      f = yee.ez(i,j,k)+yee.ez(i,j,k-iz)+      dz*(yee.ez(i  ,j  ,k+iz) - yee.ez(i  ,j  ,k-iz));
      f+=                                      dx*(yee.ez(i+1,j  ,k   ) + yee.ez(i+1,j  ,k-iz)+
                                               dz*(yee.ez(i+1,j  ,k+iz) - yee.ez(i+1,j  ,k-iz))-f);
      g = yee.ez(i,j+1,k)+ yee.ez(i,j+1,k-iz)+ dz*(yee.ez(i  ,j+1,k+iz) - yee.ez(i  ,j+1,k-iz));
      g+=                                      dx*(yee.ez(i+1,j+1,k   ) + yee.ez(i+1,j+1,k-iz)+
                                               dz*(yee.ez(i+1,j+1,k+iz) - yee.ez(i+1,j+1,k-iz))-g);
      ez[n]=(f+dy*(g-f))*0.5;

      f = yee.bx(i,j-1,k)  +yee.bx(i,j-1,k-iz )   +dz*(yee.bx(i,j-1,k+iz)   - yee.bx(i,j-1,k-iz));
      f = yee.bx(i,j,k)    +yee.bx(i,j,k-iz)      +dz*(yee.bx(i,j,k+iz)     - yee.bx(i,j,k-iz))+f+dy 
        *(yee.bx(i,j+1,k)  +yee.bx(i,j+1,k-iz)    +dz*(yee.bx(i,j+1,k+iz)   - yee.bx(i,j+1,k-iz))-f);
      g = yee.bx(i+1,j-1,k)+yee.bx(i+1,j-1,k-iz)  +dz*(yee.bx(i+1,j-1,k+iz) - yee.bx(i+1,j-1,k-iz));
      g = yee.bx(i+1,j,k)  +yee.bx(i+1,j,k-iz)    +dz*(yee.bx(i+1,j,k+iz)   - yee.bx(i+1,j,k-iz))
                                           +g     +dy*(yee.bx(i+1,j+1,k)    + yee.bx(i+1,j+1,k-iz)
                                                  +dz*(yee.bx(i+1,j+1,k+iz) - yee.bx(i+1,j+1,k-iz))-g);
      bx[n]=(f+dx*(g-f))*(.25);

      f = yee.by(i,j,k-iz)+yee.by(i-1,j,k-iz)     +dx*(yee.by(i+1,j,k-iz)   - yee.by(i-1,j,k-iz));
      f = yee.by(i,j,k)+yee.by(i-1,j,k)           +dx*(yee.by(i+1,j,k)      - yee.by(i-1,j,k))+f+dz 
        *(yee.by(i,j,k+iz)+yee.by(i-1,j,k+iz)     +dx*(yee.by(i+1,j,k+iz)   - yee.by(i-1,j,k+iz))-f);
      g = yee.by(i,j+1,k-iz)+yee.by(i-1,j+1,k-iz) +dx*(yee.by(i+1,j+1,k-iz) - yee.by(i-1,j+1,k-iz));
      g = yee.by(i,j+1,k)+yee.by(i-1,j+1,k)       +dx*(yee.by(i+1,j+1,k)    - yee.by(i-1,j+1,k))
                                               +g +dz*(yee.by(i,j+1,k+iz)   + yee.by(i-1,j+1,k+iz)
                                                  +dx*(yee.by(i+1,j+1,k+iz) - yee.by(i-1,j+1,k+iz))-g);
      by[n]=(f+dy*(g-f))*(.25);

      f = yee.bz(i-1,j,k)+yee.bz(i-1,j-1,k )      +dy*(yee.bz(i-1,j+1,k)    - yee.bz(i-1,j-1,k));
      f = yee.bz(i,j,k)+yee.bz(i,j-1,k)           +dy*(yee.bz(i,j+1,k)      - yee.bz(i,j-1,k))+f+dx 
        * (yee.bz(i+1,j,k)+yee.bz(i+1,j-1,k)      +dy*(yee.bz(i+1,j+1,k)    - yee.bz(i+1,j-1,k))-f);
      g = yee.bz(i-1,j, k+iz)+yee.bz(i-1,j-1,k+iz)+dy*(yee.bz(i-1,j+1,k+iz) - yee.bz(i-1,j-1,k+iz));
      g = yee.bz(i,j,k+iz)+yee.bz(i,j-1,k+iz )    +dy*(yee.bz(i,j+1,k+iz)   - yee.bz(i,j-1,k+iz))
                                               +g +dx*(yee.bz(i+1,j,k+iz)   + yee.bz(i+1,j-1,k+iz)
                                                  +dy*(yee.bz(i+1,j+1,k+iz) - yee.bz(i+1,j-1,k+iz))-g);
      bz[n]=(f+dz*(g-f))*(.25);
    }

  } // end of loop over species

  return;
}


//--------------------------------------------------
// explicit template instantiation

//template class pic::LinearInterpolator<1,3>; // 1D3V
template class pic::LinearInterpolator<2,3>; // 2D3V
template class pic::LinearInterpolator<3,3>; // 3D3V

