#include "linear.h"

#include <cmath> 
#include <cassert>

template<size_t D, size_t V>
void pic::LinearInterpolator<D,V>::solve(
    pic::Tile<D>& tile)
{

  // get reference to the Yee grid 
  auto& yee = tile.get_yee();

  for(auto&& container : tile.containers) {

    int nparts = container.size();

    // initialize pointers to particle arrays
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    /// resize internal arrays
    container.Epart.resize(3*nparts);
    container.Bpart.resize(3*nparts);
      
    double *ex, *ey, *ez, *bx, *by, *bz;
    ex = &( container.Epart[0*nparts] );
    ey = &( container.Epart[1*nparts] );
    ez = &( container.Epart[2*nparts] );

    bx = &( container.Bpart[0*nparts] );
    by = &( container.Bpart[1*nparts] );
    bz = &( container.Bpart[2*nparts] );


    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    int i=0, j=0, k=0;
    double dx=0.0, dy=0.0, dz=0.0;
    double f,g;

    int iz = 1;
    if (D<=2) iz = 0; // flip switch for making array queries 2D

    auto mins = tile.mins;
    //auto maxs = tile.maxs;

    // TODO: think SIMD (not possible due to ijk writing to yee)
    for(int n=n1; n<n2; n++) {

      // particle location in the grid
	    if (D >= 1) i  = static_cast<int>(floor( loc[0][n] - mins[0] ));
	    if (D >= 2) j  = static_cast<int>(floor( loc[1][n] - mins[1] ));
	    if (D >= 3) k  = static_cast<int>(floor( loc[2][n] - mins[2] ));

	    if (D >= 1) dx = (loc[0][n]-mins[0]) - i;
	    if (D >= 2) dy = (loc[1][n]-mins[1]) - j;
	    if (D >= 3) dz = (loc[2][n]-mins[2]) - k;


      // check section; TODO; remove
      bool debug_flag = false;
      if(D >= 1) { if(! (i >= 0 && i <= static_cast<int>(tile.mesh_lengths[0]) )) debug_flag = true;}
      if(D >= 2) { if(! (j >= 0 && j <= static_cast<int>(tile.mesh_lengths[1]) )) debug_flag = true;}
      if(D >= 3) { if(! (k >= 0 && k <= static_cast<int>(tile.mesh_lengths[2]) )) debug_flag = true;}

      if(debug_flag) {
        std::cout << "--------------------------------------------------\n";
        std::cout << "n=" << n;
        std::cout << " i: " << i;
        std::cout << " j: " << j;
        std::cout << " k: " << k;
        std::cout << "\n";

        std::cout << " mins0: " << mins[0];
        std::cout << " mins1: " << mins[1];
        std::cout << " mins2: " << mins[2];

        std::cout << " x: " << loc[0][n];
        std::cout << " y: " << loc[1][n];
        std::cout << " z: " << loc[2][n];
        std::cout << "\n";

        std::cout << " dx: " << dx;
        std::cout << " dy: " << dy;
        std::cout << " dz: " << dz;
        std::cout << "\n";

        std::cout << std::flush;
        // always fail
        assert(false);
      }


      // TODO: these can be optimized further when we know D
      // TODO: can also be optimized by using 1D indexing (since we can pre-compute index with
      // l = i + iy*(j-1)+iz*(k-1)
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
template class pic::LinearInterpolator<3,3>; // 3D3V //TODO; validate

