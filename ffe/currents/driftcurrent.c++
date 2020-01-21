#include "driftcurrent.h"
#include "../../em-fields/tile.h"
#include "../../tools/signum.h"

#include <cmath>



template<>
void ffe::DriftCurrent<2>::interpolate_b(
    fields::YeeLattice& yee
    )
{
  double f,g;

  int k = 0;
  for(int j=0; j<static_cast<int>(yee.Ny); j++) 
  for(int i=0; i<static_cast<int>(yee.Nx); i++) {

    bxf(i,j,k)=0.25*(
      yee.bx(i,   j,   k) + 
      yee.bx(i,   j-1, k) +
      yee.bx(i+1, j,   k) +
      yee.bx(i+1, j-1, k) 
      );

    byf(i,j,k) = 0.25*(
        yee.by(i,   j,   k) + 
        yee.by(i,   j+1, k) + 
        yee.by(i-1, j,   k) + 
        yee.by(i-1, j+1, k)
        );

    bzf(i,j,k) = 0.25*(
        yee.bz(i,   j,   k) + 
        yee.bz(i,   j-1, k) + 
        yee.bz(i-1, j,   k) + 
        yee.bz(i-1, j-1, k)
        );
  }
}


template<>
void ffe::DriftCurrent<3>::interpolate_b(
    fields::YeeLattice& yee
    )
{
  double f,g;

  for(int k=0; k<static_cast<int>(yee.Nz); k++) 
  for(int j=0; j<static_cast<int>(yee.Ny); j++) 
  for(int i=0; i<static_cast<int>(yee.Nx); i++) {

  bxf(i,j,k) = 0.125*(
    yee.bx(i,   j,   k  ) + 
    yee.bx(i,   j,   k-1) + 
    yee.bx(i,   j-1, k  ) + 
    yee.bx(i,   j-1, k-1) + 
    yee.bx(i+1, j,   k  ) + 
    yee.bx(i+1, j,   k-1) + 
    yee.bx(i+1, j-1, k  ) + 
    yee.bx(i+1, j-1, k-1)
  );

  byf(i,j,k) = 0.125*(
    yee.by(i,   j,   k)   + 
    yee.by(i,   j,   k-1) + 
    yee.by(i,   j+1, k)   + 
    yee.by(i,   j+1, k-1) + 
    yee.by(i-1, j,   k)   + 
    yee.by(i-1, j,   k-1) + 
    yee.by(i-1, j+1, k)   + 
    yee.by(i-1, j+1, k-1)
  );

  bzf(i,j,k) = 0.125*(
    yee.bz(i,   j,   k)   + 
    yee.bz(i,   j,   k+1) + 
    yee.bz(i,   j-1, k)   + 
    yee.bz(i,   j-1, k+1) + 
    yee.bz(i-1, j,   k)   + 
    yee.bz(i-1, j,   k+1) + 
    yee.bz(i-1, j-1, k)   + 
    yee.bz(i-1, j-1, k+1)
  );

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
  double dt = tile.cfl; 

  double dive;
  double crossx, crossy, crossz;
  double b2;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // div E
    // nabla E = dEx/dx + dEy/dx + dEz/dx
    //
    // NOTE: dx = 1.0 by definition
    // NOTE: dEz = 0 in 2D
    dive = 0.5*( mesh.ex(i+1,j,k) - mesh.ex(i-1,j,k) 
               + mesh.ey(i,j+1,k) - mesh.ey(i,j-1,k)
            // + mesh.ez(i,j,k+1) - mesh.ez(i,j,k-1)
               );

    // E x B
    crossx = mesh.ey(i,j,k)*bzf(i,j,k) - mesh.ez(i,j,k)*byf(i,j,k);
    crossy =-mesh.ex(i,j,k)*bzf(i,j,k) + mesh.ez(i,j,k)*bxf(i,j,k);
    crossz = mesh.ex(i,j,k)*byf(i,j,k) - mesh.ey(i,j,k)*bxf(i,j,k);

    // B^2
    b2 = bxf(i,j,k)*bxf(i,j,k) 
       + byf(i,j,k)*byf(i,j,k) 
       + bzf(i,j,k)*bzf(i,j,k);

    // perp current
    mesh.jx(i,j,k) += dive*crossx/b2;
    mesh.jy(i,j,k) += dive*crossy/b2;
    mesh.jz(i,j,k) += dive*crossz/b2;

    mesh.ex(i,j,k) -= mesh.jx(i,j,k)*dt;
    mesh.ey(i,j,k) -= mesh.jy(i,j,k)*dt;
    mesh.ez(i,j,k) -= mesh.jz(i,j,k)*dt;


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
  double b2;

  double dt = tile.cfl; // time step is basically = cfl

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

    // J_parallel x dt
    jparax = spara*bxf(i,j,k)/b2;
    jparay = spara*byf(i,j,k)/b2;
    jparaz = spara*bzf(i,j,k)/b2;

    // Add J_parallel to E
    mesh.ex(i,j,k) -= jparax;
    mesh.ey(i,j,k) -= jparay;
    mesh.ez(i,j,k) -= jparaz;

    // store J_parallel
    mesh.jx(i,j,k) += jparax/dt;
    mesh.jy(i,j,k) += jparay/dt;
    mesh.jz(i,j,k) += jparaz/dt;

  }
}


/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::limiter(
    ffe::Tile<2>& tile
    )
{
  fields::YeeLattice& mesh = tile.get_yee();
  interpolate_b(mesh);

  double diss;
  double e2, b2;

  double dt = tile.cfl;

  double ex0, ey0, ez0;
  double ex1, ey1, ez1;
  double jxd, jyd, jzd;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    ex0 = mesh.ex(i,j,k);
    ey0 = mesh.ey(i,j,k);
    ez0 = mesh.ez(i,j,k);

    // component-wise resistivity
    //if( abs(ex0) > abs(bxf(i,j,k))) ex1 = toolbox::sign(ex0)*abs(bxf(i,j,k));
    //if( abs(ey0) > abs(byf(i,j,k))) ey1 = toolbox::sign(ey0)*abs(byf(i,j,k));
    //if( abs(ez0) > abs(bzf(i,j,k))) ez1 = toolbox::sign(ez0)*abs(bzf(i,j,k));

    //// E^2 and B^2
    e2 = ex0*ex0 + ey0*ey0 + ez0*ez0;
    b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k);

    // limiter by sqrt(b^2/e^2)
    diss = sqrt(b2/e2);
    if(diss < 1.0) {
      ex1 = diss*ex0;
      ey1 = diss*ey0;
      ez1 = diss*ez0;

      // dissipated current
      jxd = ex0 - ex1;
      jyd = ey0 - ey1;
      jzd = ez0 - ez1;

      // deposit dissipating current
      mesh.ex(i,j,k) -= jxd;
      mesh.ey(i,j,k) -= jyd;
      mesh.ez(i,j,k) -= jzd;

    } else {
      jxd = 0.0;
      jyd = 0.0;
      jzd = 0.0;
    }

    // store into mesh
    mesh.jx1(i,j,k) = jxd/dt;
    mesh.jy1(i,j,k) = jyd/dt;
    mesh.jz1(i,j,k) = jzd/dt;

    // E^2
    //e2 = mesh.ex(i,j,k)*mesh.ex(i,j,k) 
    //   + mesh.ey(i,j,k)*mesh.ey(i,j,k) 
    //   + mesh.ez(i,j,k)*mesh.ez(i,j,k);

    //// B^2
    //b2 = bxf(i,j,k)*bxf(i,j,k) 
    //   + byf(i,j,k)*byf(i,j,k) 
    //   + bzf(i,j,k)*bzf(i,j,k);

    //if(e2 > b2) {
    //  mesh.ex(i,j,k) *= sqrt(b2/e2);
    //  mesh.ey(i,j,k) *= sqrt(b2/e2);
    //  mesh.ez(i,j,k) *= sqrt(b2/e2);
    //}


  }
}







