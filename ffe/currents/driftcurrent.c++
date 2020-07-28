#include "driftcurrent.h"
#include "../../em-fields/tile.h"
#include "../../tools/signum.h"

#include <cmath>

#include <nvtx3/nvToolsExt.h> 


template<>
void ffe::DriftCurrent<2>::interpolate_b(
    fields::YeeLattice& yee
    )
{
  nvtxRangePush(__FUNCTION__);

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
  nvtxRangePop();

}


template<>
void ffe::DriftCurrent<3>::interpolate_b(
    fields::YeeLattice& yee
    )
{
  nvtxRangePush(__FUNCTION__);

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
  nvtxRangePop();

}


/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::comp_drift_cur(ffe::Tile<2>& tile)
{
  nvtxRangePush(__FUNCTION__);

  fields::YeeLattice& mesh = tile.get_yee();

  // half stagger B; stored in (bxf, byf, bzf)
  interpolate_b(mesh);

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
    divE = 0.5*( mesh.ex(i+1,j,k) - mesh.ex(i-1,j,k) 
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

    if(b2 == 0.0) continue;

    // perp current
    mesh.jx(i,j,k) +=   divE*crossx/b2;
    mesh.jy(i,j,k) +=   divE*crossy/b2;
    mesh.jz(i,j,k) +=   divE*crossz/b2;

    mesh.ex(i,j,k) -=dt*divE*crossx/b2;
    mesh.ey(i,j,k) -=dt*divE*crossy/b2;
    mesh.ez(i,j,k) -=dt*divE*crossz/b2;

    //mesh.ex(i,j,k) -= mesh.jx(i,j,k)*dt;
    //mesh.ey(i,j,k) -= mesh.jy(i,j,k)*dt;
    //mesh.ez(i,j,k) -= mesh.jz(i,j,k)*dt;

  }
  nvtxRangePop();

}


/// 3D Drift current solver
template<>
void ffe::DriftCurrent<3>::comp_drift_cur(ffe::Tile<3>& tile)
{
  nvtxRangePush(__FUNCTION__);

  fields::YeeLattice& mesh = tile.get_yee();

  // half stagger B; stored in (bxf, byf, bzf)
  interpolate_b(mesh);

  double dt = tile.cfl; 

  double dive;
  double crossx, crossy, crossz;
  double b2;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) 
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // div E
    // nabla E = dEx/dx + dEy/dx + dEz/dx
    //
    // NOTE: dx = 1.0 by definition
    dive = 0.5*( mesh.ex(i+1,j,k) - mesh.ex(i-1,j,k) 
               + mesh.ey(i,j+1,k) - mesh.ey(i,j-1,k)
               + mesh.ez(i,j,k+1) - mesh.ez(i,j,k-1)
               );

    // E x B
    crossx = mesh.ey(i,j,k)*bzf(i,j,k) - mesh.ez(i,j,k)*byf(i,j,k);
    crossy =-mesh.ex(i,j,k)*bzf(i,j,k) + mesh.ez(i,j,k)*bxf(i,j,k);
    crossz = mesh.ex(i,j,k)*byf(i,j,k) - mesh.ey(i,j,k)*bxf(i,j,k);

    // B^2
    b2 = bxf(i,j,k)*bxf(i,j,k) 
       + byf(i,j,k)*byf(i,j,k) 
       + bzf(i,j,k)*bzf(i,j,k);

    if(b2 == 0.0) continue;

    // perp current
    mesh.jx(i,j,k) += dive*crossx/b2;
    mesh.jy(i,j,k) += dive*crossy/b2;
    mesh.jz(i,j,k) += dive*crossz/b2;

    mesh.ex(i,j,k) -=dt*dive*crossx/b2;
    mesh.ey(i,j,k) -=dt*dive*crossy/b2;
    mesh.ez(i,j,k) -=dt*dive*crossz/b2;

  }
  nvtxRangePop();

}

/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::comp_parallel_cur(
    ffe::Tile<2>& tile
    )
{
  nvtxRangePush(__FUNCTION__);

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

    if(b2 == 0.0) continue;

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
  nvtxRangePop();

}

/// 3D Drift current solver
template<>
void ffe::DriftCurrent<3>::comp_parallel_cur(
    ffe::Tile<3>& tile
    )
{
  nvtxRangePush(__FUNCTION__);

  fields::YeeLattice& mesh = tile.get_yee();

  // half stagger B; stored in (bxf, byf, bzf)
  interpolate_b(mesh);

  double spara;
  double jparax, jparay, jparaz;
  double b2;

  double dt = tile.cfl; // time step is basically = cfl

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) 
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

    if(b2 == 0.0) continue;

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
  nvtxRangePop();

}

/// 2D Drift current solver
template<>
void ffe::DriftCurrent<2>::limiter(
    ffe::Tile<2>& tile
    )
{
  nvtxRangePush(__FUNCTION__);

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

      mesh.jx(i,j,k) += jxd/dt;
      mesh.jy(i,j,k) += jyd/dt;
      mesh.jz(i,j,k) += jzd/dt;

    } else {
      jxd = 0.0;
      jyd = 0.0;
      jzd = 0.0;
    }

  }
  nvtxRangePop();

}



/// 3D Drift current solver
template<>
void ffe::DriftCurrent<3>::limiter(
    ffe::Tile<3>& tile
    )
{
  nvtxRangePush(__FUNCTION__);

  fields::YeeLattice& mesh = tile.get_yee();
  interpolate_b(mesh);

  double diss;
  double e2, b2;

  double dt = tile.cfl;

  double ex0, ey0, ez0;
  double ex1, ey1, ez1;
  double jxd, jyd, jzd;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) 
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    ex0 = mesh.ex(i,j,k);
    ey0 = mesh.ey(i,j,k);
    ez0 = mesh.ez(i,j,k);

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

      mesh.jx(i,j,k) += jxd/dt;
      mesh.jy(i,j,k) += jyd/dt;
      mesh.jz(i,j,k) += jzd/dt;

    } else {
      jxd = 0.0;
      jyd = 0.0;
      jzd = 0.0;
    }

  }
  nvtxRangePop();

}


//--------------------------------------------------
// explicit template instantiation
template class ffe::DriftCurrent<2>; // 2D
template class ffe::DriftCurrent<3>; // 3D
