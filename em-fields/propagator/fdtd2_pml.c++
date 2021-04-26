#include "fdtd2_pml.h"

#include <cmath>
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

/// Radial PML damping coefficient profile
template<>
DEVCALLABLE float_m fields::FDTD2_pml<3>::lambda( float_m sx, float_m sy, float_m sz )
{
  // normalized unit radius 
  float_m r;

  // radial 
  if (mode == 0) {
    r = std::sqrt( 
            pow( (sx - cenx)/radx, 2) 
          + pow( (sy - ceny)/rady, 2) 
          + pow( (sz - cenz)/radz, 2) 
               );

  // cylindrical in x
  } else if(mode == 1) {
    r = std::sqrt( 
          + pow( (sy - ceny)/rady, 2) 
          + pow( (sz - cenz)/radz, 2) 
               );

  // cylindrical in z
  } else if(mode == 2) {
    r = std::sqrt( 
            pow( (sx - cenx)/radx, 2) 
          + pow( (sy - ceny)/rady, 2) 
               );
  }

  if(r > rad_lim) {
    return -std::min(
            norm_abs*pow( (r - rad_lim)/(1.0 - rad_lim), 3),
            norm_abs
            );
  } else {
    return 0.0;
  }
}



/// 3D E pusher
template<>
void fields::FDTD2_pml<3>::push_e(fields::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  YeeLattice& mesh = tile.get_yee();
  const float_m C = 1.0 * tile.cfl * dt * corr;
  const auto mins = tile.mins;

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, YeeLattice &mesh)
  {
   
    //-----------
    // global grid coordinates
    float_m iglob = i + mins[0];
    float_m jglob = j + mins[1];
    float_m kglob = k + mins[2];
    float_m lam;

    // dE = dt*curl B  (1/2)
	lam = 0.5*lambda(iglob+0.5, jglob, kglob);
    mesh.ex(i,j,k) = (
            mesh.ex(i,j,k)*(1+lam) + 
            C*(mesh.by(i,j,k-1)-mesh.by(i,j,k)  - mesh.bz(i,j-1,k)+mesh.bz(i,j,k)) 
            )/(1-lam);

	lam = 0.5*lambda(iglob, jglob+0.5, kglob);
    mesh.ey(i,j,k) = (
            mesh.ey(i,j,k)*(1+lam) + 
            C*(mesh.bz(i-1,j,k)-mesh.bz(i,j,k)  - mesh.bx(i,j,k-1)+mesh.bx(i,j,k)) 
            )/(1-lam);

	lam = 0.5*lambda(iglob, jglob, kglob+0.5);
    mesh.ez(i,j,k) = (
            mesh.ez(i,j,k)*(1+lam) + 
            C*(mesh.bx(i,j-1,k)-mesh.bx(i,j,k)  - mesh.by(i-1,j,k)+mesh.by(i,j,k)) 
            )/(1-lam);

  }, 
    tile.mesh_lengths[0], 
    tile.mesh_lengths[1], 
    tile.mesh_lengths[2], 
    mesh);

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------

/// 3D B pusher
template<>
void fields::FDTD2_pml<3>::push_half_b(fields::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  YeeLattice& mesh = tile.get_yee();
  const float_m C = 0.5 * tile.cfl * dt * corr;
  const auto mins = tile.mins;

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, YeeLattice &mesh)
  {
    // global grid coordinates
    float_m iglob = i + mins[0];
    float_m jglob = j + mins[1];
    float_m kglob = k + mins[2];
    float_m lam;

    // dB = -dt*curl E  (1/4)
	lam = 0.25*lambda(iglob, jglob+0.5, kglob+0.5);
	mesh.bx(i,j,k) = ( mesh.bx(i,j,k)*(1+lam)
            + C*(mesh.ey(i,j,k+1)-mesh.ey(i,j,k)  - mesh.ez(i,j+1,k)+mesh.ez(i,j,k)) 
            )/(1-lam);

	lam = 0.25*lambda(iglob+0.5, jglob, kglob+0.5);
	mesh.by(i,j,k) = ( mesh.by(i,j,k)*(1+lam)
            + C*(mesh.ez(i+1,j,k)-mesh.ez(i,j,k)  - mesh.ex(i,j,k+1)+mesh.ex(i,j,k)) 
            )/(1-lam);

	lam = 0.25*lambda(iglob+0.5, jglob+0.5, kglob);
	mesh.bz(i,j,k) = ( mesh.bz(i,j,k)*(1+lam)
            + C*(mesh.ex(i,j+1,k)-mesh.ex(i,j,k)  - mesh.ey(i+1,j,k)+mesh.ey(i,j,k)) 
            )/(1-lam);
      
  }, 
    tile.mesh_lengths[0], 
    tile.mesh_lengths[1], 
    tile.mesh_lengths[2], 
    mesh);

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}


/// 3D E & B coupled pusher (for FFE)
template<>
void fields::FDTD2_pml<3>::push_eb(::ffe::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // refs to storages
  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 
    
  // Tile limits
  const auto mins = tile.mins;
  //const auto maxs = tile.maxs;

  const float_m c = tile.cfl; // dt / dx
  const float_m C1 =  c; // high-order curl operator coefficients

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, 
      ffe::SkinnyYeeLattice& dm,
      fields::YeeLattice& m
      )
  {
    //-----------
    // global grid coordinates
    float_m iglob = i + mins[0];
    float_m jglob = j + mins[1];
    float_m kglob = k + mins[2];
    float_m lam;

    //-----------
    // simultaneous E & B update with perfectly matched layer damping

    // dB = -dt*curl E  (1/4)
	lam = 0.5*lambda(iglob, jglob+0.5, kglob+0.5);
	dm.bx(i,j,k) = (m.bx(i,j,k)*lam + C1*(m.ey(i,j,k+1)-m.ey(i,j,k)  - m.ez(i,j+1,k)+m.ez(i,j,k)) )/(1-lam);

	lam = 0.5*lambda(iglob+0.5, jglob, kglob+0.5);
	dm.by(i,j,k) = (m.by(i,j,k)*lam + C1*(m.ez(i+1,j,k)-m.ez(i,j,k)  - m.ex(i,j,k+1)+m.ex(i,j,k)) )/(1-lam);

	lam = 0.5*lambda(iglob+0.5, jglob+0.5, kglob);
	dm.bz(i,j,k) = (m.bz(i,j,k)*lam + C1*(m.ex(i,j+1,k)-m.ex(i,j,k)  - m.ey(i+1,j,k)+m.ey(i,j,k)) )/(1-lam);


    // dE = dt*curl B  (1/2)
	lam = 0.5*lambda(iglob+0.5, jglob, kglob);
    dm.ex(i,j,k) = (m.ex(i,j,k)*lam + C1*(m.by(i,j,k-1)-m.by(i,j,k)  - m.bz(i,j-1,k)+m.bz(i,j,k)) )/(1-lam);

	lam = 0.5*lambda(iglob, jglob+0.5, kglob);
    dm.ey(i,j,k) = (m.ey(i,j,k)*lam + C1*(m.bz(i-1,j,k)-m.bz(i,j,k)  - m.bx(i,j,k-1)+m.bx(i,j,k)) )/(1-lam);

	lam = 0.5*lambda(iglob, jglob, kglob+0.5);
    dm.ez(i,j,k) = (m.ez(i,j,k)*lam + C1*(m.bx(i,j-1,k)-m.bx(i,j,k)  - m.by(i-1,j,k)+m.by(i,j,k)) )/(1-lam);

  }, 
    tile.mesh_lengths[0], 
    tile.mesh_lengths[1], 
    tile.mesh_lengths[2], 
    dm, m
    );

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}


//template class fields::FDTD2_pml<1>;
//template class fields::FDTD2_pml<2>;
template class fields::FDTD2_pml<3>; 
