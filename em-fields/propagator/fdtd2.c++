#include "fdtd2.h"

#include <cmath>

/* 
 * 1D version:
 *	ey(i,j,k)=ey(i,j,k) + c *(bz(im1,j,k)-bz(i,j,k))  
 *	ez(i,j,k)=ez(i,j,k) + c *(by(i,j,k)-by(im1,j,k)) 
 *
 * 2D version
 * ex(i,j,k)=ex(i,j,k)+const*(-bz(i,jm1,k)+bz(i,j,k))
 * ey(i,j,k)=ey(i,j,k)+const*(bz(im1,j,k)-bz(i,j,k))
 * ez(i,j,k)=ez(i,j,k)+const*(bx(i,jm1,k)-bx(i,j,k)-  by(im1,j,k)+by(i,j,k))
*/

/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */

#include <nvtx3/nvToolsExt.h> 

#include "../../tools/iter/iter.h"

/// 1D E pusher
template<>
void fields::FDTD2<1>::push_e(fields::Tile<1>& tile)
{
  YeeLattice& mesh = tile.get_yee();
  Realf C = 1.0 * tile.cfl * dt * corr;

  int k = 0;
  int j = 0;
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Ex
    // NONE

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k));

  }

}

/// 2D E pusher
template<>
void fields::FDTD2<2>::push_e(fields::Tile<2>& tile)
{
  YeeLattice& mesh = tile.get_yee();

  Realf C = 1.0 * tile.cfl * dt * corr;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Ex
    mesh.ex(i,j,k) += 
      + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k) 
           -mesh.by(i-1,j,   k) + mesh.by(i,j,k));

  }
}


/// 3D E pusher
template<>
void fields::FDTD2<3>::push_e(fields::Tile<3>& tile)
{
nvtxRangePush(__PRETTY_FUNCTION__);

  YeeLattice& mesh = tile.get_yee();
  Realf C = 1.0 * tile.cfl * dt * corr;


    UniIter::iterate3D(
    [=] DEVCALLABLE (int i, int j, int k, YeeLattice &mesh)
    {

    // Ex
    mesh.ex(i,j,k) += 
      + C*( mesh.by(i,j,  k-1) - mesh.by(i,j,k))
      + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k))
      + C*(-mesh.bx(i,  j, k-1) + mesh.bx(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k))
      + C*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k));
    }, static_cast<int>(tile.mesh_lengths[2]), static_cast<int>(tile.mesh_lengths[1]), static_cast<int>(tile.mesh_lengths[0]), mesh);
UniIter::sync();
    
nvtxRangePop();

}


//--------------------------------------------------

/*
 * 1D version:
  by(i,j,k)=by(i,j,k)+const*(ez(ip1,j,k)-ez(i,j,k)) !-0*ex(i,j,k+1)+0*ex(i,j,k))
	bz(i,j,k)=bz(i,j,k)+const*(ey(i,j,k)-ey(ip1,j,k)) !+0*ex(i,j+1,k)-0*ex(i,j,k))
  
 * 2D version:
 * bx(i,j,k)=bx(i,j,k)+const*(-ez(i,jp1,k)+ez(i,j,k))
 * by(i,j,k)=by(i,j,k)+const*(ez(ip1,j,k)-ez(i,j,k))
 * bz(i,j,k)=bz(i,j,k)+const*(ex(i,jp1,k)-ex(i,j,k) -ey(ip1,j,k)+ey(i,j,k))
*/

/// Update B field with a half step


/// 1D B pusher
template<>
void fields::FDTD2<1>::push_half_b(fields::Tile<1>& tile)
{
  YeeLattice& mesh = tile.get_yee();
  Realf C = 0.5 * tile.cfl * dt * corr;

  int k = 0;
  int j = 0;
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Bx
    // NONE

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));
  }

}

/// 2D B pusher
template<>
void fields::FDTD2<2>::push_half_b(fields::Tile<2>& tile)
{
  YeeLattice& mesh = tile.get_yee();

  Realf C = 0.5 * tile.cfl * dt * corr;

  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++)
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    // Bx
    mesh.bx(i,j,k) += 
      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)
           -mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

  }

}


/// 3D B pusher
template<>
void fields::FDTD2<3>::push_half_b(fields::Tile<3>& tile)
{
nvtxRangePush(__PRETTY_FUNCTION__);

  YeeLattice& mesh = tile.get_yee();
  Realf C = 0.5 * tile.cfl * dt * corr;

    UniIter::iterate3D(
    [=] DEVCALLABLE (int i, int j, int k, YeeLattice &mesh)
    {

    // Bx
    mesh.bx(i,j,k) += 
      + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k))
      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k))
      + C*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k))
      + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

      }, static_cast<int>(tile.mesh_lengths[2]), static_cast<int>(tile.mesh_lengths[1]), static_cast<int>(tile.mesh_lengths[0]), mesh);
UniIter::sync();
    
nvtxRangePop();

}



template class fields::FDTD2<1>;
template class fields::FDTD2<2>;
template class fields::FDTD2<3>;
  
