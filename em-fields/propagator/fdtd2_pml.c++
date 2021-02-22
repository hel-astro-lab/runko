#include "fdtd2_pml.h"

#include <cmath>


/// PML damping coefficient profile
template<>
real_short fields::FDTD2_pml<3>::lambda( real_short sx, real_short sy, real_short sz )
{
  // normalized unit radius 
  real_short r = std::sqrt( 
            pow( (sx - cenx)/radx, 2) 
          + pow( (sy - ceny)/rady, 2) 
          + pow( (sz - cenz)/radz, 2) 
               );

  if(r > rad_lim) {
    return -std::min(
            norm_abs*pow( (r - rad_lim)/(1.0 - rad_lim), 3),
            norm_abs
            );
  } else {
    return 0.0;
  }
}


/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */


/// 2D E pusher
//template<>
//void fields::FDTD2_pml<2>::push_e(fields::Tile<2>& tile)
//{
//  YeeLattice& mesh = tile.get_yee();
//
//  Realf C = 1.0 * tile.cfl * dt * corr;
//
//  int k = 0;
//  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
//  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
//
//    // Ex
//    mesh.ex(i,j,k) += 
//      + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));
//
//    // Ey
//    mesh.ey(i,j,k) += 
//      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));
//
//    // Ez
//    mesh.ez(i,j,k) += 
//      + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k) 
//           -mesh.by(i-1,j,   k) + mesh.by(i,j,k));
//
//  }
//}

/// 3D E pusher
template<>
void fields::FDTD2_pml<3>::push_e(fields::Tile<3>& tile)
{
  YeeLattice& mesh = tile.get_yee();
  Realf C = 1.0 * tile.cfl * dt * corr;

  auto mins = tile.mins;
  real_short lam, iglob, jglob, kglob;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++)
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++)
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Ex
    //mesh.ex(i,j,k) += 
    //  + C*( mesh.by(i,j,  k-1) - mesh.by(i,j,k))
    //  + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

    // Ey
    //mesh.ey(i,j,k) += 
    //  + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k))
    //  + C*(-mesh.bx(i,  j, k-1) + mesh.bx(i,j,k));

    // Ez
    //mesh.ez(i,j,k) += 
    //  + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k))
    //  + C*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k));

    //-----------
    // global grid coordinates
    iglob = static_cast<real_short>(i) + mins[0];
    jglob = static_cast<real_short>(j) + mins[1];
    kglob = static_cast<real_short>(k) + mins[2];

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

  }
}


//--------------------------------------------------

/// 2D B pusher
//template<>
//void fields::FDTD2_pml<2>::push_half_b(fields::Tile<2>& tile)
//{
//  YeeLattice& mesh = tile.get_yee();
//
//  Realf C = 0.5 * tile.cfl * dt * corr;
//
//  int k = 0;
//  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++)
//  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
//
//    // Bx
//    mesh.bx(i,j,k) += 
//      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));
//
//    // By
//    mesh.by(i,j,k) += 
//      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k));
//
//    // Bz
//    mesh.bz(i,j,k) += 
//      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)
//           -mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));
//
//  }
//
//}

/// 3D B pusher
template<>
void fields::FDTD2_pml<3>::push_half_b(fields::Tile<3>& tile)
{
  YeeLattice& mesh = tile.get_yee();
  Realf C = 0.5 * tile.cfl * dt * corr;

  auto mins = tile.mins;
  real_short lam, iglob, jglob, kglob;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) 
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Bx
    //mesh.bx(i,j,k) += 
    //  + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k))
    //  + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    //mesh.by(i,j,k) += 
    //  + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k))
    //  + C*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k));

    // Bz
    //mesh.bz(i,j,k) += 
    //  + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k))
    //  + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

    //-----------
    // global grid coordinates
    iglob = static_cast<real_short>(i) + mins[0];
    jglob = static_cast<real_short>(j) + mins[1];
    kglob = static_cast<real_short>(k) + mins[2];

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
      
  }
}


/// 3D E & B coupled pusher (for FFE)
template<>
void fields::FDTD2_pml<3>::push_eb(::ffe::Tile<3>& tile)
{
  // refs to storages
  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 
    
  // Tile limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;

  // refs to fields for easier access
  auto& ex  = m.ex;
  auto& ey  = m.ey;
  auto& ez  = m.ez;

  auto& bx  = m.bx;
  auto& by  = m.by;
  auto& bz  = m.bz;

  // dt / dx
  real_short c = tile.cfl;

  // high-order curl operator coefficients
  real_short C1 =  c;

  // pml parameters
  real_short lam, iglob, jglob, kglob;


  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        // original non-damped reference version
		//dm.bx(i,j,k) = C1*(ey(i,j,k+1)-ey(i,j,k)  - ez(i,j+1,k)+ez(i,j,k));
		//dm.by(i,j,k) = C1*(ez(i+1,j,k)-ez(i,j,k)  - ex(i,j,k+1)+ex(i,j,k));
		//dm.bz(i,j,k) = C1*(ex(i,j+1,k)-ex(i,j,k)  - ey(i+1,j,k)+ey(i,j,k));

        //dm.ex(i,j,k) = C1*(by(i,j,k-1)-by(i,j,k)  - bz(i,j-1,k)+bz(i,j,k));
        //dm.ey(i,j,k) = C1*(bz(i-1,j,k)-bz(i,j,k)  - bx(i,j,k-1)+bx(i,j,k));
        //dm.ez(i,j,k) = C1*(bx(i,j-1,k)-bx(i,j,k)  - by(i-1,j,k)+by(i,j,k));

        //-----------
        // global grid coordinates
        iglob = static_cast<real_short>(i) + mins[0];
        jglob = static_cast<real_short>(j) + mins[1];
        kglob = static_cast<real_short>(k) + mins[2];

        //-----------
        // simultaneous E & B update with perfectly matched layer damping

        // dB = -dt*curl E  (1/4)
		lam = 0.5*lambda(iglob, jglob+0.5, kglob+0.5);
		dm.bx(i,j,k) = (bx(i,j,k)*lam + C1*(ey(i,j,k+1)-ey(i,j,k)  - ez(i,j+1,k)+ez(i,j,k)) )/(1-lam);

	    lam = 0.5*lambda(iglob+0.5, jglob, kglob+0.5);
		dm.by(i,j,k) = (by(i,j,k)*lam + C1*(ez(i+1,j,k)-ez(i,j,k)  - ex(i,j,k+1)+ex(i,j,k)) )/(1-lam);

	    lam = 0.5*lambda(iglob+0.5, jglob+0.5, kglob);
		dm.bz(i,j,k) = (bz(i,j,k)*lam + C1*(ex(i,j+1,k)-ex(i,j,k)  - ey(i+1,j,k)+ey(i,j,k)) )/(1-lam);


        // dE = dt*curl B  (1/2)
		lam = 0.5*lambda(iglob+0.5, jglob, kglob);
        dm.ex(i,j,k) = (ex(i,j,k)*lam + C1*(by(i,j,k-1)-by(i,j,k)  - bz(i,j-1,k)+bz(i,j,k)) )/(1-lam);

		lam = 0.5*lambda(iglob, jglob+0.5, kglob);
        dm.ey(i,j,k) = (ey(i,j,k)*lam + C1*(bz(i-1,j,k)-bz(i,j,k)  - bx(i,j,k-1)+bx(i,j,k)) )/(1-lam);

		lam = 0.5*lambda(iglob, jglob, kglob+0.5);
        dm.ez(i,j,k) = (ez(i,j,k)*lam + C1*(bx(i,j-1,k)-bx(i,j,k)  - by(i-1,j,k)+by(i,j,k)) )/(1-lam);
      }
    }
  }
}


//template class fields::FDTD2_pml<1>;
//template class fields::FDTD2_pml<2>;
template class fields::FDTD2_pml<3>;
  
