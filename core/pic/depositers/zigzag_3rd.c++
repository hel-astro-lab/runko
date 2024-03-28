#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/pic/depositers/zigzag_3rd.h"
#include "core/pic/shapes.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;



template<size_t D, size_t V>
void pic::ZigZag_3rd<D,V>::solve( pic::Tile<D>& tile )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& gs = tile.get_grids();
  const auto mins = tile.mins;

  //clear arrays before new update
  gs.jx.clear();
  gs.jy.clear();
  gs.jz.clear();


  for(auto&& con : tile.containers) {

    const double c = tile.cfl;    // speed of light
    const double q = con.q; // charge
                            //
    // skip particle species if zero charge
    if (q == 0.0) continue;

    //UniIter::iterate([=] DEVCALLABLE (
    //            size_t n, 
    //            emf::Grids &gs,
    //            pic::ParticleContainer<D>& con
    //            ){
    for(size_t n=0; n<con.size(); n++) {

      //--------------------------------------------------
      double u = con.vel(0,n);
      double v = con.vel(1,n);
      double w = con.vel(2,n);
      double invgam = 1.0/sqrt(1.0 + u*u + v*v + w*w);

      //--------------------------------------------------
      // new (normalized) location, x_{n+1}
        
      //double x1, y1, z1;
      double x2 = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      double y2 = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      double z2 = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      // previos location, x_n
      double x1 = x2 - u*invgam*c;
      double y1 = y2 - v*invgam*c;
      double z1 = z2 - w*invgam*c; 

      //--------------------------------------------------
      int i1  = D >= 1 ? floor(x1) : 0;
      int i2  = D >= 1 ? floor(x2) : 0;
      int j1  = D >= 2 ? floor(y1) : 0;
      int j2  = D >= 2 ? floor(y2) : 0;
      int k1  = D >= 3 ? floor(z1) : 0;
      int k2  = D >= 3 ? floor(z2) : 0;

      // 1st order relay point; +1 is equal to +\Delta x
      //double xr = min( double(min(i1,i2)+1), max( double(max(i1,i2)), double(0.5*(x1+x2)) ) );
      //double yr = min( double(min(j1,j2)+1), max( double(max(j1,j2)), double(0.5*(y1+y2)) ) );
      //double zr = min( double(min(k1,k2)+1), max( double(max(k1,k2)), double(0.5*(z1+z2)) ) );
        
      // 2nd order relay point; +1 is equal to +\Delta x
      double xr = min( double(min(i1,i2)+1), max( double(i1+i2)*0.5, double(0.5*(x1+x2)) ) );
      double yr = min( double(min(j1,j2)+1), max( double(j1+j2)*0.5, double(0.5*(y1+y2)) ) );
      double zr = min( double(min(k1,k2)+1), max( double(k1+k2)*0.5, double(0.5*(z1+z2)) ) );

      //--------------------------------------------------

      // \Delta x from grid points
      double dx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.5; 
      double dx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.5; 
                                                    
      double dy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.5; 
      double dy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.5; 
                                                    
      double dz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.5; 
      double dz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.5; 

      //--------------------------------------------------
      // Lorentz contract lenghts
      // https://physics.stackexchange.com/questions/56078/lorentz-boost-matrix-for-an-arbitrary-direction-in-terms-of-rapidity/588284
      //double gam = sqrt(1.0 + u*u + v*v + w*w); // \gamma
      //double betax2 = pow(u/gam, 2);            // v_x^2
      //double betay2 = pow(v/gam, 2);            // v_y^2
      //double betaz2 = pow(w/gam, 2);            // v_z^2
      //double beta2  = betax2 + betay2 + betaz2; // |v|^2

      //dx1 = dx1/(1.0 + (gam-1.0)*betax2/(beta2 + EPS));
      //dx2 = dx2/(1.0 + (gam-1.0)*betax2/(beta2 + EPS));
      //          
      //dy1 = dy1/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dy2 = dy2/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //          
      //dz1 = dz1/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));
      //dz2 = dz2/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));


      //--------------------------------------------------
      // particle weights 
      double Wxx1[5] = {1.0}, 
             Wxx2[5] = {1.0}, 
             Wyy1[5] = {1.0}, 
             Wyy2[5] = {1.0}, 
             Wzz1[5] = {1.0}, 
             Wzz2[5] = {1.0};

      if(D >= 1) W3rd(dx1, Wxx1);
      if(D >= 1) W3rd(dx2, Wxx2);

      if(D >= 2) W3rd(dy1, Wyy1);
      if(D >= 2) W3rd(dy2, Wyy2);

      if(D >= 3) W3rd(dz1, Wzz1);
      if(D >= 3) W3rd(dz2, Wzz2);

      //--------------------------------------------------
      // staggered grid along the motion
      double Wx1[3] = {1.0}, 
             Wx2[3] = {1.0},
             Wy1[3] = {1.0},
             Wy2[3] = {1.0},
             Wz1[3] = {1.0},
             Wz2[3] = {1.0};

      // Shifting \Delta x by -0.5 to accommodate staggering
      if(D >= 1) W2nd(dx1-0.5, Wx1);
      if(D >= 1) W2nd(dx2-0.5, Wx2);

      if(D >= 2) W2nd(dy1-0.5, Wy1);
      if(D >= 2) W2nd(dy2-0.5, Wy2);

      if(D >= 3) W2nd(dz1-0.5, Wz1);
      if(D >= 3) W2nd(dz2-0.5, Wz2);


      //--------------------------------------------------
      // q v = q (x_{i+1} - x_i)/dt
      //
      // NOTE: +q since - sign is already included in the Ampere's equation
      // NOTE: extra c to introduce time step immediately; therefore we store on grid J -> J\Delta t
      // NOTE: More generally we should have: q = weight*qe;
      double qvx1 = D >= 1 ? +q*(xr - x1) : +q*c;
      double qvy1 = D >= 2 ? +q*(yr - y1) : +q*c;
      double qvz1 = D >= 3 ? +q*(zr - z1) : +q*c;

      double qvx2 = D >= 1 ? +q*(x2 - xr) : +q*c;
      double qvy2 = D >= 2 ? +q*(y2 - yr) : +q*c;
      double qvz2 = D >= 3 ? +q*(z2 - zr) : +q*c;

      //--------------------------------------------------
      // dimension dependent loop boundaries
      const int xlim = D >= 1 ? 1 : 0;
      const int ylim = D >= 2 ? 1 : 0;
      const int zlim = D >= 3 ? 1 : 0;
        

      // NOTE: incrementing loop counter with ++i to enforce at least 1 iteration always
      // NOTE: limits are from -zlim to zlim+1 to get (-1,0,1,2) 3rd order access pattern
        
      //jx
      for(int zi=-zlim; zi <=zlim+1; ++zi)
      for(int yi=-ylim; yi <=ylim+1; ++yi){

        //std::cout << "jx: injecting into" <<
        //"(" << i1-1 <<","<< j1+yi <<","<< k1+zi <<") " <<
        //"(" << i1   <<","<< j1+yi <<","<< k1+zi <<") " <<
        //"(" << i2-1 <<","<< j2+yi <<","<< k2+zi <<") " <<
        //"(" << i2   <<","<< j2+yi <<","<< k2+zi <<") " << "\n";

        if(D >= 2) atomic_add( gs.jx(i1-1, j1+yi, k1+zi), qvx1* Wx1[0]* Wyy1[yi+1]*Wzz1[zi+1] );
        if(D >= 1) atomic_add( gs.jx(i1  , j1+yi, k1+zi), qvx1* Wx1[1]* Wyy1[yi+1]*Wzz1[zi+1] );
        if(D >= 2) atomic_add( gs.jx(i1+1, j1+yi, k1+zi), qvx1* Wx1[2]* Wyy1[yi+1]*Wzz1[zi+1] );


        if(D >= 2) atomic_add( gs.jx(i2-1, j2+yi, k2+zi), qvx2* Wx2[0]* Wyy2[yi+1]*Wzz2[zi+1] );
        if(D >= 1) atomic_add( gs.jx(i2,   j2+yi, k2+zi), qvx2* Wx2[1]* Wyy2[yi+1]*Wzz2[zi+1] );
        if(D >= 2) atomic_add( gs.jx(i2+1, j2+yi, k2+zi), qvx2* Wx2[2]* Wyy2[yi+1]*Wzz2[zi+1] );
      }


      //jy
      for(int zi=-zlim; zi <=zlim+1; ++zi)
      for(int xi=-xlim; xi <=xlim+1; ++xi){
        //std::cout << "jy: injecting into" <<
        //"(" << i1+xi <<","<< j1-1 <<","<< k1+zi <<") " <<
        //"(" << i1+xi <<","<< j1   <<","<< k1+zi <<") " <<
        //"(" << i2+xi <<","<< j2-1 <<","<< k2+zi <<") " <<
        //"(" << i2+xi <<","<< j2   <<","<< k2+zi <<") " << "\n";

        if(D >= 2) atomic_add( gs.jy(i1+xi, j1-1, k1+zi), qvy1* Wy1[0] *Wxx1[xi+1]*Wzz1[zi+1] );
        if(D >= 1) atomic_add( gs.jy(i1+xi, j1,   k1+zi), qvy1* Wy1[1] *Wxx1[xi+1]*Wzz1[zi+1] );
        if(D >= 2) atomic_add( gs.jy(i1+xi, j1+1, k1+zi), qvy1* Wy1[2] *Wxx1[xi+1]*Wzz1[zi+1] );


        if(D >= 2) atomic_add( gs.jy(i2+xi, j2-1, k2+zi), qvy2* Wy2[0] *Wxx2[xi+1]*Wzz2[zi+1] );
        if(D >= 1) atomic_add( gs.jy(i2+xi, j2,   k2+zi), qvy2* Wy2[1] *Wxx2[xi+1]*Wzz2[zi+1] );
        if(D >= 2) atomic_add( gs.jy(i2+xi, j2+1, k2+zi), qvy2* Wy2[2] *Wxx2[xi+1]*Wzz2[zi+1] );
      }
                                                                                       
      //jz                                                                             
      for(int yi=-ylim; yi <=ylim+1; ++yi)                                                     
      for(int xi=-xlim; xi <=xlim+1; ++xi){                                                    
        //std::cout << "jz: injecting into" <<                                         
        //"(" << i1+xi <<","<< j1+yi <<","<< k1-1 <<") " <<                            
        //"(" << i1+xi <<","<< j1+yi <<","<< k1   <<") " <<                            
        //"(" << i2+xi <<","<< j2+yi <<","<< k2-1 <<") " <<                            
        //"(" << i2+xi <<","<< j2+yi <<","<< k2   <<") " << "\n";                      
                                                                                       
        if(D >= 3) atomic_add( gs.jz(i1+xi, j1+yi, k1-1), qvz1* Wz1[0] *Wxx1[xi+1]*Wyy1[yi+1] );
        if(D >= 1) atomic_add( gs.jz(i1+xi, j1+yi, k1  ), qvz1* Wz1[1] *Wxx1[xi+1]*Wyy1[yi+1] );
        if(D >= 3) atomic_add( gs.jz(i1+xi, j1+yi, k1+1), qvz1* Wz1[2] *Wxx1[xi+1]*Wyy1[yi+1] );

        if(D >= 3) atomic_add( gs.jz(i2+xi, j2+yi, k2-1), qvz2* Wz2[0] *Wxx2[xi+1]*Wyy2[yi+1] );
        if(D >= 1) atomic_add( gs.jz(i2+xi, j2+yi, k2  ), qvz2* Wz2[1] *Wxx2[xi+1]*Wyy2[yi+1] );
        if(D >= 3) atomic_add( gs.jz(i2+xi, j2+yi, k2+1), qvz2* Wz2[2] *Wxx2[xi+1]*Wyy2[yi+1] );
      }

    }
    //}, con.size(), gs, con);

  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag_3rd<1,3>; // 1D3V
template class pic::ZigZag_3rd<2,3>; // 2D3V
template class pic::ZigZag_3rd<3,3>; // 3D3V

