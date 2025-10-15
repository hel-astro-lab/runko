#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/pic/depositers/zigzag_2nd.h"
#include "core/pic/shapes.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;



template<size_t D, size_t V>
void pic::ZigZag_2nd<D,V>::solve( pic::Tile<D>& tile )
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
      double x2 = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      double y2 = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      double z2 = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      // previos location, x_n
      double x1 = x2 - u*invgam*c;
      double y1 = y2 - v*invgam*c;
      double z1 = z2 - w*invgam*c; 

      //--------------------------------------------------
      // primary grid; -1/2 to +1/2
      int i1p = D >= 1 ? round(x1) : 0;
      int i2p = D >= 1 ? round(x2) : 0;
      int j1p = D >= 2 ? round(y1) : 0;
      int j2p = D >= 2 ? round(y2) : 0;
      int k1p = D >= 3 ? round(z1) : 0;
      int k2p = D >= 3 ? round(z2) : 0;

      // dual grid 0 to 1
      int i1d = D >= 1 ? floor(x1) : 0;
      int i2d = D >= 1 ? floor(x2) : 0;
      int j1d = D >= 2 ? floor(y1) : 0;
      int j2d = D >= 2 ? floor(y2) : 0;
      int k1d = D >= 3 ? floor(z1) : 0;
      int k2d = D >= 3 ? floor(z2) : 0;


      // 1st order relay point; +1 is equal to +\Delta x
      //double xr = min( double(min(i1,i2)+1), max( double(max(i1,i2)), double(0.5*(x1+x2)) ) );
      //double yr = min( double(min(j1,j2)+1), max( double(max(j1,j2)), double(0.5*(y1+y2)) ) );
      //double zr = min( double(min(k1,k2)+1), max( double(max(k1,k2)), double(0.5*(z1+z2)) ) );
        
      // 2nd order primary grid relay point
      double xr = min( double(min(i1p,i2p)+1), max( double(i1p+i2p)*0.5, double(0.5*(x1+x2)) ) );
      double yr = min( double(min(j1p,j2p)+1), max( double(j1p+j2p)*0.5, double(0.5*(y1+y2)) ) );
      double zr = min( double(min(k1p,k2p)+1), max( double(k1p+k2p)*0.5, double(0.5*(z1+z2)) ) );


      // 1st order staggered grid relay point
      double xrd = min( double(min(i1d,i2d)+1), max( double(max(i1d,i2d)), double(0.5*(x1+x2)) ) );
      double yrd = min( double(min(j1d,j2d)+1), max( double(max(j1d,j2d)), double(0.5*(y1+y2)) ) );
      double zrd = min( double(min(k1d,k2d)+1), max( double(max(k1d,k2d)), double(0.5*(z1+z2)) ) );


      //--------------------------------------------------
      // \Delta x from grid points
      double dx1 = D >= 1 ? 0.5*(x1 + xr) - i1p : 0.5; 
      double dx2 = D >= 1 ? 0.5*(x2 + xr) - i2p : 0.5; 

      double dy1 = D >= 2 ? 0.5*(y1 + yr) - j1p : 0.5; 
      double dy2 = D >= 2 ? 0.5*(y2 + yr) - j2p : 0.5; 

      double dz1 = D >= 3 ? 0.5*(z1 + zr) - k1p : 0.5; 
      double dz2 = D >= 3 ? 0.5*(z2 + zr) - k2p : 0.5; 

      //--------------------------------------------------
      // staggered grid
      // note: relay point already has staggering so no -0.5 factor here
      double dx1d = D >= 1 ? 0.5*(x1 + xrd) - i1d : 0.0; 
      double dx2d = D >= 1 ? 0.5*(x2 + xrd) - i2d : 0.0; 
      double dy1d = D >= 2 ? 0.5*(y1 + yrd) - j1d : 0.0; 
      double dy2d = D >= 2 ? 0.5*(y2 + yrd) - j2d : 0.0; 
      double dz1d = D >= 3 ? 0.5*(z1 + zrd) - k1d : 0.0; 
      double dz2d = D >= 3 ? 0.5*(z2 + zrd) - k2d : 0.0; 

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
      //dy1 = dy1/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dy2 = dy2/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dz1 = dz1/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));
      //dz2 = dz2/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));

      //dx1d = dx1d/(1.0 + (gam-1.0)*betax2/(beta2 + EPS));
      //dx2d = dx2d/(1.0 + (gam-1.0)*betax2/(beta2 + EPS));
      //dy1d = dy1d/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dy2d = dy2d/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dz1d = dz1d/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));
      //dz2d = dz2d/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));

      //--------------------------------------------------
      // particle weights on primary grid
      double Wxx1[3] = {1.0}, Wxx2[3] = {1.0}, Wyy1[3] = {1.0}, Wyy2[3] = {1.0}, Wzz1[3] = {1.0}, Wzz2[3] = {1.0};

      if(D >= 1) W2nd(dx1, Wxx1);
      if(D >= 1) W2nd(dx2, Wxx2);

      if(D >= 2) W2nd(dy1, Wyy1);
      if(D >= 2) W2nd(dy2, Wyy2);

      if(D >= 3) W2nd(dz1, Wzz1);
      if(D >= 3) W2nd(dz2, Wzz2);

      //--------------------------------------------------
      // staggered grid along the motion
      double Wx1[3] = {1.0}, Wx2[3] = {1.0}, Wy1[3] = {1.0}, Wy2[3] = {1.0}, Wz1[3] = {1.0}, Wz2[3] = {1.0};


      // Shifting \Delta x by -0.5 to accommodate staggering
      // DONE working ver
      if(D >= 1) W1st(dx1d, Wx1);
      if(D >= 1) W1st(dx2d, Wx2);
      if(D >= 2) W1st(dy1d, Wy1);
      if(D >= 2) W1st(dy2d, Wy2);
      if(D >= 3) W1st(dz1d, Wz1);
      if(D >= 3) W1st(dz2d, Wz2);


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


      //Wx1_ip1 = 0.5f*(x + xr) - i; 
      //jx(i-1, j-1, k) = qvx * (0.5 - Wx1_ip1)*W2_jm1
      //jx(i  , j-1, k) = qvx * (0.5 + Wx1_ip1)*W2_jm1
      //jx(i-1, j  , k) = qvx * (0.5 - Wx1_ip1)*W2_j  
      //jx(i  , j  , k) = qvx * (0.5 + Wx1_ip1)*W2_j  
      //jx(i-1, j+1, k) = qvx * (0.5 - Wx1_ip1)*W2_jp1
      //jx(i  , j+1, k) = qvx * (0.5 + Wx1_ip1)*W2_jp1


#ifdef DEBUG

      bool debug = false;
      for(int iii=0; iii<3; iii++){
          if(Wxx1[iii] < 0.0)  debug = true;
          if(Wxx2[iii] < 0.0)  debug = true;

          if(Wyy1[iii] < 0.0)  debug = true;
          if(Wyy2[iii] < 0.0)  debug = true;

          if(Wzz1[iii] < 0.0)  debug = true;
          if(Wzz2[iii] < 0.0)  debug = true;

          if(Wx1[iii] < 0.0)  debug = true;
          if(Wy1[iii] < 0.0)  debug = true;
          if(Wz1[iii] < 0.0)  debug = true;
          if(Wx2[iii] < 0.0)  debug = true;
          if(Wy2[iii] < 0.0)  debug = true;
          if(Wz2[iii] < 0.0)  debug = true;
      }

      if(D<=2) debug=false;
        
      if(debug){
        std::cout 
          << "cur xyz: " << "(" << xr << "," << yr << "," << zr << ")"
          << " i1: "  << "(" << i1p << "," << j1p << "," << k1p << ")"
          << " i2: "  << "(" << i2p << "," << j2p << "," << k2p << ")"
          << " i1d: "  << "(" << i1d << "," << j1d << "," << k1d << ")"
          << " i2d: "  << "(" << i2d << "," << j2d << "," << k2d << ")"
          << " x1: " << "(" << x1 << "," << y1 << "," << z1 << ")"
          << " x2: " << "(" << x2 << "," << y2 << "," << z2 << ")"
          << " dx1: " << "(" << dx1 << "," << dy1 << "," << dz1 << ")"
          << " dx2: " << "(" << dx2 << "," << dy2 << "," << dz2 << ")"
          << " dx1d: "<< "(" << dx1d << "," << dy1d << "," << dz1d << ")"
          << " dx2d: "<< "(" << dx2d << "," << dy2d << "," << dz2d << ")"
          << " W1x: " << "(" << Wxx1[0] << "," << Wxx1[1] << "," << Wxx1[2] << /* "," << Wxx1[3] << */ ")"
          << " W2x: " << "(" << Wxx2[0] << "," << Wxx2[1] << "," << Wxx2[2] << /* "," << Wxx2[3] << */ ")"
          << " W1y: " << "(" << Wyy1[0] << "," << Wyy1[1] << "," << Wyy1[2] << /* "," << Wyy1[3] << */ ")"
          << " W2y: " << "(" << Wyy2[0] << "," << Wyy2[1] << "," << Wyy2[2] << /* "," << Wyy2[3] << */ ")"
          << " W1z: " << "(" << Wzz1[0] << "," << Wzz1[1] << "," << Wzz1[2] << /* "," << Wzz1[3] << */ ")"
          << " W2z: " << "(" << Wzz2[0] << "," << Wzz2[1] << "," << Wzz2[2] << /* "," << Wzz2[3] << */ ")"
          << " W1x1: "<< "(" << Wx1[0] << "," << Wx1[1] << "," << Wx1[2]    << /* "," << Wx1[3]  << */ ")"
          << " W1y1: "<< "(" << Wy1[0] << "," << Wy1[1] << "," << Wy1[2]    << /* "," << Wy1[3]  << */ ")"
          << " W1z1: "<< "(" << Wz1[0] << "," << Wz1[1] << "," << Wz1[2]    << /* "," << Wz1[3]  << */ ")"
          << " W2x1: "<< "(" << Wx2[0] << "," << Wx2[1] << "," << Wx2[2]    << /* "," << Wx1[3]  << */ ")"
          << " W2y1: "<< "(" << Wy2[0] << "," << Wy2[1] << "," << Wy2[2]    << /* "," << Wy1[3]  << */ ")"
          << " W2z1: "<< "(" << Wz2[0] << "," << Wz2[1] << "," << Wz2[2]    << /* "," << Wz1[3]  << */ ")"
          << "\n";
      }
#endif

      // TODO: what about original scheme? Does it equal this?
      // TODO: this still lacks -1/2 staggering which is strange
      // TODO: Sokolev correction to interpolator too
      // TODO: is there a more general formula that can be applied to calculate these?

      // NOTE: incrementing loop counter with ++i to enforce at least 1 iteration always
      //jx
      for(int zi=-zlim; zi <=zlim; ++zi)
      for(int yi=-ylim; yi <=ylim; ++yi){

        //if(debug){
        //  std::cout << "jx: injecting into" <<
        //  "(" << i1-1 <<","<< j1+yi <<","<< k1+zi <<") " <<
        //  "(" << i1   <<","<< j1+yi <<","<< k1+zi <<") " <<
        //  "(" << i2-1 <<","<< j2+yi <<","<< k2+zi <<") " <<
        //  "(" << i2   <<","<< j2+yi <<","<< k2+zi <<") " << "\n";
        //}

        //if(D >= 1) atomic_add( gs.jx(i1d  , j1p+yi, k1p+zi), qvx1* Wx1[1]* Wyy1[yi+1]*Wzz1[zi+1] );
        //if(D >= 1) atomic_add( gs.jx(i1d+1, j1p+yi, k1p+zi), qvx1* Wx1[2]* Wyy1[yi+1]*Wzz1[zi+1] );
        //if(D >= 1) atomic_add( gs.jx(i2d  , j2p+yi, k2p+zi), qvx2* Wx2[1]* Wyy2[yi+1]*Wzz2[zi+1] );
        //if(D >= 1) atomic_add( gs.jx(i2d+1, j2p+yi, k2p+zi), qvx2* Wx2[2]* Wyy2[yi+1]*Wzz2[zi+1] );

        if(D >= 1) atomic_add( gs.jx(i1d-1, j1p+yi, k1p+zi), qvx1* Wx1[1]* Wyy1[yi+1]*Wzz1[zi+1] );
        if(D >= 1) atomic_add( gs.jx(i1d  , j1p+yi, k1p+zi), qvx1* Wx1[2]* Wyy1[yi+1]*Wzz1[zi+1] );
        if(D >= 1) atomic_add( gs.jx(i2d-1, j2p+yi, k2p+zi), qvx2* Wx2[1]* Wyy2[yi+1]*Wzz2[zi+1] );
        if(D >= 1) atomic_add( gs.jx(i2d  , j2p+yi, k2p+zi), qvx2* Wx2[2]* Wyy2[yi+1]*Wzz2[zi+1] );

      }

      //jy
      for(int zi=-zlim; zi <=zlim; ++zi)
      for(int xi=-xlim; xi <=xlim; ++xi){
        //std::cout << "jy: injecting into" <<
        //"(" << i1+xi <<","<< j1-1 <<","<< k1+zi <<") " <<
        //"(" << i1+xi <<","<< j1   <<","<< k1+zi <<") " <<
        //"(" << i2+xi <<","<< j2-1 <<","<< k2+zi <<") " <<
        //"(" << i2+xi <<","<< j2   <<","<< k2+zi <<") " << "\n";

        //if(D >= 2) atomic_add( gs.jy(i1p+xi, j1d  , k1p+zi), qvy1* Wy1[1] *Wxx1[xi+1]*Wzz1[zi+1] );
        //if(D >= 1) atomic_add( gs.jy(i1p+xi, j1d+1, k1p+zi), qvy1* Wy1[2] *Wxx1[xi+1]*Wzz1[zi+1] );
        //if(D >= 2) atomic_add( gs.jy(i2p+xi, j2d  , k2p+zi), qvy2* Wy2[1] *Wxx2[xi+1]*Wzz2[zi+1] );
        //if(D >= 1) atomic_add( gs.jy(i2p+xi, j2d+1, k2p+zi), qvy2* Wy2[2] *Wxx2[xi+1]*Wzz2[zi+1] );

        if(D >= 2) atomic_add( gs.jy(i1p+xi, j1d-1, k1p+zi), qvy1* Wy1[1] *Wxx1[xi+1]*Wzz1[zi+1] );
        if(D >= 1) atomic_add( gs.jy(i1p+xi, j1d  , k1p+zi), qvy1* Wy1[2] *Wxx1[xi+1]*Wzz1[zi+1] );
        if(D >= 2) atomic_add( gs.jy(i2p+xi, j2d-1, k2p+zi), qvy2* Wy2[1] *Wxx2[xi+1]*Wzz2[zi+1] );
        if(D >= 1) atomic_add( gs.jy(i2p+xi, j2d  , k2p+zi), qvy2* Wy2[2] *Wxx2[xi+1]*Wzz2[zi+1] );
      }                                                                                
                                                                                       
      //jz                                                                             
      for(int yi=-ylim; yi <=ylim; ++yi)                                                     
      for(int xi=-xlim; xi <=xlim; ++xi){                                                    
        //std::cout << "jz: injecting into" <<                                         
        //"(" << i1+xi <<","<< j1+yi <<","<< k1-1 <<") " <<                            
        //"(" << i1+xi <<","<< j1+yi <<","<< k1   <<") " <<                            
        //"(" << i2+xi <<","<< j2+yi <<","<< k2-1 <<") " <<                            
        //"(" << i2+xi <<","<< j2+yi <<","<< k2   <<") " << "\n";                      
                                                                                       
        //if(D >= 3) atomic_add( gs.jz(i1p+xi, j1p+yi, k1d  ), qvz1* Wz1[1] *Wxx1[xi+1]*Wyy1[yi+1] );
        //if(D >= 1) atomic_add( gs.jz(i1p+xi, j1p+yi, k1d+1), qvz1* Wz1[2] *Wxx1[xi+1]*Wyy1[yi+1] );
        //if(D >= 3) atomic_add( gs.jz(i2p+xi, j2p+yi, k2d  ), qvz2* Wz2[1] *Wxx2[xi+1]*Wyy2[yi+1] );
        //if(D >= 1) atomic_add( gs.jz(i2p+xi, j2p+yi, k2d+1), qvz2* Wz2[2] *Wxx2[xi+1]*Wyy2[yi+1] );

        if(D >= 3) atomic_add( gs.jz(i1p+xi, j1p+yi, k1d-1), qvz1* Wz1[1] *Wxx1[xi+1]*Wyy1[yi+1] );
        if(D >= 1) atomic_add( gs.jz(i1p+xi, j1p+yi, k1d  ), qvz1* Wz1[2] *Wxx1[xi+1]*Wyy1[yi+1] );
        if(D >= 3) atomic_add( gs.jz(i2p+xi, j2p+yi, k2d-1), qvz2* Wz2[1] *Wxx2[xi+1]*Wyy2[yi+1] );
        if(D >= 1) atomic_add( gs.jz(i2p+xi, j2p+yi, k2d  ), qvz2* Wz2[2] *Wxx2[xi+1]*Wyy2[yi+1] );
      }

    }
    //}, con.size(), gs, con);

  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag_2nd<1,3>; // 1D3V
template class pic::ZigZag_2nd<2,3>; // 2D3V
template class pic::ZigZag_2nd<3,3>; // 3D3V

