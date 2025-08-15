#include <cmath>

#include "core/emf/propagators/fdtd4.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */



/// 1D E pusher
template<>
void emf::FDTD4<1>::push_e(emf::Tile<1>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();

  const float C1 = coeff1*corr*tile.cfl;
  const float C2 = coeff2*corr*tile.cfl;

  UniIter::iterate(
  [=] DEVCALLABLE (int i, Grids &mesh)
  {
  // Ex NONE
	mesh.ey(i,0,0) +=   C1*( mesh.bz(i-1, 0, 0) - mesh.bz(i,   0, 0)) 
		                 +C2*( mesh.bz(i-2, 0, 0) - mesh.bz(i+1, 0, 0));

	mesh.ez(i,0,0) +=   C1*(-mesh.by(i-1, 0, 0) + mesh.by(i,   0, 0))
                     +C2*(-mesh.by(i-2, 0, 0) + mesh.by(i+1, 0, 0));
  }, 
    tile.mesh_lengths[0], 
    mesh);

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}


/// 2D E pusher
template<>
void emf::FDTD4<2>::push_e(emf::Tile<2>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();

  const float C1 = coeff1*corr*tile.cfl;
  const float C2 = coeff2*corr*tile.cfl;
  const int k = 0;

  UniIter::iterate2D(
  [=] DEVCALLABLE (int i, int j, Grids &mesh)
  {
	mesh.ex(i,j,k) += C1*(-mesh.bz(i,  j-1,k) + mesh.bz(i,  j,  k)) 
		               +C2*(-mesh.bz(i,  j-2,k) + mesh.bz(i,  j+1,k));

	mesh.ey(i,j,k) += C1*( mesh.bz(i-1,j,  k) - mesh.bz(i,  j,  k)) 
		               +C2*( mesh.bz(i-2,j,  k) - mesh.bz(i+1,j,  k));

	mesh.ez(i,j,k) += C1*( mesh.bx(i,  j-1,k) - mesh.bx(i,  j,  k) 
                       - mesh.by(i-1,j,  k) + mesh.by(i,  j,  k))
                   +C2*( mesh.bx(i,  j-2,k) - mesh.bx(i,  j+1,k)
                        -mesh.by(i-2,j,  k) + mesh.by(i+1,j,  k));
  }, 
    tile.mesh_lengths[0], 
    tile.mesh_lengths[1], 
    mesh);

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}

/// 3D E pusher
template<>
void emf::FDTD4<3>::push_e(emf::Tile<3>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C1 = coeff1*corr*tile.cfl;
  const float C2 = coeff2*corr*tile.cfl;

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, Grids &mesh)
  {
    mesh.ex(i,j,k)+= C1*(mesh.by(i,  j,  k-1) - mesh.by(i,  j,  k  )  
                       - mesh.bz(i,  j-1,k  ) + mesh.bz(i,  j,  k  ))
                    +C2*(mesh.by(i,  j,  k-2) - mesh.by(i,  j,  k+1)
                       - mesh.bz(i,  j-2,k  ) + mesh.bz(i,  j+1,k  ));

    mesh.ey(i,j,k)+= C1*(mesh.bz(i-1,j,  k  ) - mesh.bz(i,  j,  k  )  
                       - mesh.bx(i,  j,  k-1) + mesh.bx(i,  j,  k  ))
                    +C2*(mesh.bz(i-2,j,  k  ) - mesh.bz(i+1,j,  k  )
                       - mesh.bx(i,  j,  k-2) + mesh.bx(i,  j,  k+1));

    mesh.ez(i,j,k)+= C1*(mesh.bx(i,  j-1,k  ) - mesh.bx(i,  j,  k  )  
                       - mesh.by(i-1,j,  k  ) + mesh.by(i,  j,  k  ))
                    +C2*(mesh.bx(i,  j-2,k  ) - mesh.bx(i,  j+1,k  )
                       - mesh.by(i-2,j,  k  ) + mesh.by(i+1,j,  k  ));
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

/// Update B field with a half step


/// 1D B pusher
template<>
void emf::FDTD4<1>::push_half_b(emf::Tile<1>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C1 = 0.5*coeff1*corr*tile.cfl;
  const float C2 = 0.5*coeff2*corr*tile.cfl;

  UniIter::iterate(
  [=] DEVCALLABLE (int i, Grids &mesh)
  {
    // Bx NONE
	   mesh.by(i,0,0) += C1*( mesh.ez(i+1, 0, 0) - mesh.ez(i,   0, 0))
	                    +C2*( mesh.ez(i+2, 0, 0) - mesh.ez(i-1, 0, 0));
	   mesh.bz(i,0,0) += C1*(-mesh.ey(i+1, 0, 0) + mesh.ey(i,   0, 0))
                      +C2*(-mesh.ey(i+2, 0, 0) + mesh.ey(i-1, 0, 0));

  }, 
    tile.mesh_lengths[0], 
    mesh);

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}



/// 2D B pusher
template<>
void emf::FDTD4<2>::push_half_b(emf::Tile<2>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C1 = 0.5*coeff1*corr*tile.cfl;
  const float C2 = 0.5*coeff2*corr*tile.cfl;
  const int k = 0;

  UniIter::iterate2D(
  [=] DEVCALLABLE (int i, int j, Grids &mesh)
  {
     mesh.bx(i,j,k) += C1*(-mesh.ez(i,  j+1,k) + mesh.ez(i,  j,  k))
	                  +C2*(-mesh.ez(i,  j+2,k) + mesh.ez(i,  j-1,k));
	 mesh.by(i,j,k) += C1*( mesh.ez(i+1,j,  k) - mesh.ez(i,  j,  k))
	                  +C2*( mesh.ez(i+2,j,  k) - mesh.ez(i-1,j,  k));
	 mesh.bz(i,j,k) += C1*( mesh.ex(i,  j+1,k) - mesh.ex(i,  j,  k)  
                           -mesh.ey(i+1,j,  k) + mesh.ey(i,  j,  k))
                      +C2*( mesh.ex(i,  j+2,k) - mesh.ex(i,  j-1,k)
                           -mesh.ey(i+2,j,  k) + mesh.ey(i-1,j,  k));

  }, 
    tile.mesh_lengths[0], 
    tile.mesh_lengths[1], 
    mesh);

  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}


/// 3D B pusher
template<>
void emf::FDTD4<3>::push_half_b(emf::Tile<3>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C1 = 0.5*coeff1*corr*tile.cfl;
  const float C2 = 0.5*coeff2*corr*tile.cfl;

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, Grids &mesh)
  {
    mesh.bx(i,j,k)+= C1*(mesh.ey(i,  j,  k+1) - mesh.ey(i,  j,  k  )  
                       - mesh.ez(i,  j+1,k  ) + mesh.ez(i,  j,  k  ))
                    +C2*(mesh.ey(i,  j,  k+2) - mesh.ey(i,  j,  k-1)
                       - mesh.ez(i,  j+2,k  ) + mesh.ez(i,  j-1,k  ));

	mesh.by(i,j,k)+= C1*(mesh.ez(i+1,j,  k  ) - mesh.ez(i,  j,  k  )  
                       - mesh.ex(i,  j,  k+1) + mesh.ex(i,  j,  k  ))
                    +C2*(mesh.ez(i+2,j,  k  ) - mesh.ez(i-1,j,  k  )
                        -mesh.ex(i,  j,  k+2) + mesh.ex(i,  j,  k-1));

	mesh.bz(i,j,k)+= C1*(mesh.ex(i,  j+1,k  ) - mesh.ex(i,  j,  k  )  
                       - mesh.ey(i+1,j,  k  ) + mesh.ey(i,  j,  k  ))
                    +C2*(mesh.ex(i,  j+2,k  ) - mesh.ex(i,  j-1,k  )
                       - mesh.ey(i+2,j,  k  ) + mesh.ey(i-1,j,  k  ));

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



template class emf::FDTD4<1>;
template class emf::FDTD4<2>;
template class emf::FDTD4<3>;
  
