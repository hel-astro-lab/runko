#include "driftcurrent.h"
#include "../../em-fields/tile.h"

#include <cmath>

template<>
void ffe::DriftCurrent<2>::interpolate_b(
    fields::YeeLattice& yee
    )
{
  double dx=0.0, dy=0.0, dz=0.0;
  double f,g;

  //int iz = 1;
  int iz = 0; // for D<=2; flip switch for making array queries 2D

  int k = 0;
  for(int j=0; j<static_cast<int>(yee.Ny); j++) 
  for(int i=0; i<static_cast<int>(yee.Nx); i++) {

    // XXX: lengths
    dx=+0.5;
    dy=-0.5;
    //dz=-0.5;
    f = yee.bx(i,j-1,k)  +yee.bx(i,j-1,k-iz )   +dz*(yee.bx(i,j-1,k+iz)   - yee.bx(i,j-1,k-iz));
    f = yee.bx(i,j,k)    +yee.bx(i,j,k-iz)      +dz*(yee.bx(i,j,k+iz)     - yee.bx(i,j,k-iz))+f + 
     dy*(yee.bx(i,j+1,k) +yee.bx(i,j+1,k-iz)    +dz*(yee.bx(i,j+1,k+iz)   - yee.bx(i,j+1,k-iz))-f);
    g = yee.bx(i+1,j-1,k)+yee.bx(i+1,j-1,k-iz)  +dz*(yee.bx(i+1,j-1,k+iz) - yee.bx(i+1,j-1,k-iz));
    g = yee.bx(i+1,j,k)  +yee.bx(i+1,j,k-iz)    +dz*(yee.bx(i+1,j,k+iz)   - yee.bx(i+1,j,k-iz))
                                         +g     +dy*(yee.bx(i+1,j+1,k)    + yee.bx(i+1,j+1,k-iz)
                                                +dz*(yee.bx(i+1,j+1,k+iz) - yee.bx(i+1,j+1,k-iz))-g);
    bxf(i,j,k)=(f+dx*(g-f))*(.25);

    dx=-0.5;
    dy=+0.5;
    //dz=-0.5;
    f = yee.by(i,j,k-iz)+yee.by(i-1,j,k-iz)     +dx*(yee.by(i+1,j,k-iz)   - yee.by(i-1,j,k-iz));
    f = yee.by(i,j,k)+yee.by(i-1,j,k)           +dx*(yee.by(i+1,j,k)      - yee.by(i-1,j,k))+f+dz 
      *(yee.by(i,j,k+iz)+yee.by(i-1,j,k+iz)     +dx*(yee.by(i+1,j,k+iz)   - yee.by(i-1,j,k+iz))-f);
    g = yee.by(i,j+1,k-iz)+yee.by(i-1,j+1,k-iz) +dx*(yee.by(i+1,j+1,k-iz) - yee.by(i-1,j+1,k-iz));
    g = yee.by(i,j+1,k)+yee.by(i-1,j+1,k)       +dx*(yee.by(i+1,j+1,k)    - yee.by(i-1,j+1,k))
                                             +g +dz*(yee.by(i,j+1,k+iz)   + yee.by(i-1,j+1,k+iz)
                                                +dx*(yee.by(i+1,j+1,k+iz) - yee.by(i-1,j+1,k+iz))-g);
    byf(i,j,k)=(f+dy*(g-f))*(.25);

    dx=-0.5;
    dy=-0.5;
    //dz=+0.5;
    f = yee.bz(i-1,j,k)+yee.bz(i-1,j-1,k )      +dy*(yee.bz(i-1,j+1,k)    - yee.bz(i-1,j-1,k));
    f = yee.bz(i,j,k)+yee.bz(i,j-1,k)           +dy*(yee.bz(i,j+1,k)      - yee.bz(i,j-1,k))+f+dx 
      * (yee.bz(i+1,j,k)+yee.bz(i+1,j-1,k)      +dy*(yee.bz(i+1,j+1,k)    - yee.bz(i+1,j-1,k))-f);
    g = yee.bz(i-1,j, k+iz)+yee.bz(i-1,j-1,k+iz)+dy*(yee.bz(i-1,j+1,k+iz) - yee.bz(i-1,j-1,k+iz));
    g = yee.bz(i,j,k+iz)+yee.bz(i,j-1,k+iz )    +dy*(yee.bz(i,j+1,k+iz)   - yee.bz(i,j-1,k+iz))
                                             +g +dx*(yee.bz(i+1,j,k+iz)   + yee.bz(i+1,j-1,k+iz)
                                                +dy*(yee.bz(i+1,j+1,k+iz) - yee.bz(i+1,j-1,k+iz))-g);
    bzf(i,j,k)=(f+dz*(g-f))*(.25);
  }
}


/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::comp_drift_cur(ffe::Tile<2>& tile)
{
  fields::YeeLattice& mesh = tile.get_yee();

  // half stagger B; stored in (bxf, byf, bzf)
  interpolate_b(mesh);


  //Realf C = 1.0 * tile.cfl;

  double dive=0;
  double crossx=0, crossy=0, crossz=0;
  double b2=0;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // div E
    //
    // XXX: norm fac
    dive = 1.0*(
        mesh.ex(i+1,j,k) - mesh.ex(i-1,j,k) 
      + mesh.ey(i,j+1,k) - mesh.ey(i,j-1,k)
      // + mesh.ez(i,j,k+1) - mesh.ez(i,j,k-1)
      );

    // E x B
    crossx = mesh.ey(i,j,k)*bzf(i,j,k) - mesh.ez(i,j,k)*byf(i,j,k);
    crossy =-mesh.ex(i,j,k)*bzf(i,j,k) - mesh.ez(i,j,k)*bxf(i,j,k);
    crossz = mesh.ex(i,j,k)*byf(i,j,k) - mesh.ey(i,j,k)*bxf(i,j,k);

    // B^2
    b2 = bxf(i,j,k)*bxf(i,j,k) 
       + byf(i,j,k)*byf(i,j,k) 
       + bzf(i,j,k)*bzf(i,j,k);

    // perp current
    mesh.jx(i,j,k) = dive*crossx/b2;
    mesh.jy(i,j,k) = dive*crossy/b2;
    mesh.jz(i,j,k) = dive*crossz/b2;

    /*
    std::cout << " dive  " << dive 
              << " crosx " << crossx
              << " crosy " << crossy
              << " crosz " << crossz
              << " b2 " << b2
              << " bxf " << bxf(i,j,k)
              << " byf " << byf(i,j,k)
              << " bzf " << bzf(i,j,k)
              << " jx  " << mesh.jx(i,j,k)
              << " jy  " << mesh.jy(i,j,k)
              << " jz  " << mesh.jz(i,j,k)
              << "\n";
    */

  }
}


/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::comp_parallel_cur(
    ffe::Tile<2>& tile
    )
{
  fields::YeeLattice& mesh = tile.get_yee();

  // half stagger B; stored in (bxf, byf, bzf)
  interpolate_b(mesh);

  double spara;
  double jparax, jparay, jparaz;
  double b2, e2;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // E.B
    spara = 
        mesh.ex(i,j,k)*bxf(i,j,k)
      + mesh.ey(i,j,k)*byf(i,j,k)
      + mesh.ez(i,j,k)*bzf(i,j,k);

    // B^2
    b2 = bxf(i,j,k)*bxf(i,j,k) 
       + byf(i,j,k)*byf(i,j,k) 
       + bzf(i,j,k)*bzf(i,j,k);

    // J_parallel
    jparax = spara*bxf(i,j,k)/b2;
    jparay = spara*byf(i,j,k)/b2;
    jparaz = spara*bzf(i,j,k)/b2;

    // Add J_parallel to E
    mesh.ex(i,j,k) -= jparax;
    mesh.ey(i,j,k) -= jparay;
    mesh.ez(i,j,k) -= jparaz;

    // store J_parallel
    mesh.jx1(i,j,k) = jparax;
    mesh.jy1(i,j,k) = jparay;
    mesh.jz1(i,j,k) = jparaz;

    // E^2
    e2 = mesh.ex(i,j,k)*mesh.ex(i,j,k) 
       + mesh.ey(i,j,k)*mesh.ey(i,j,k) 
       + mesh.ez(i,j,k)*mesh.ez(i,j,k);

    // enforce E < B; E -> sqrt(B^2/E^2) E
    mesh.ex(i,j,k) *= sqrt(b2/e2);
    mesh.ey(i,j,k) *= sqrt(b2/e2);
    mesh.ez(i,j,k) *= sqrt(b2/e2);

  }
}







