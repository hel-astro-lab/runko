#include <cmath>

#include "core/emf/propagators/fdtd2_pml.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


/// Radial PML damping coefficient profile
template<size_t D>
DEVCALLABLE float emf::FDTD2_pml<D>::lambda( 
    float sx, 
    float sy, 
    float sz )
{
  // normalized unit radius 
  auto drx = (sx - cenx)/radx;
  auto dry = (sy - ceny)/rady;
  auto drz = (sz - cenz)/radz;

  // radial 
  //if(mode == 0) 
  auto r = std::sqrt( drx*drx + dry*dry + drz*drz );

  // cylindrical in x
  //if(mode == 1) r = std::sqrt( dry*dry + drz*drz );

  // cylindrical in z
  //if(mode == 2) r = std::sqrt( drx*drx + dry*dry );

  // 2D plane
  //r = std::sqrt( drx*drx + dry*dry );

  //std::cout << "r" << r;

  // cylindrical region (left and right sides)
  //float rcyl = std::sqrt( drx*drx + drz*drz );
  //float lam = pow( (rcyl - rad_lim)/(1.0f - rad_lim), 3);
  //if( rcyl < rad_lim ) lam = 0.0;
    
  // NOTE turned off
  //float lam = 0.0f;
  //float z    = std::sqrt( dry*dry );
  //float lamz = pow( (z - rad_lim)/(1.0f - rad_lim), 3);
  //if( z < rad_lim ) lamz = 0.0;

  //return -norm_abs*std::min(lam + lamz, 1.0f - EPS);

  // radial
  if(r > rad_lim) {
    return -std::min(
            norm_abs*pow( (r - rad_lim)/(1.0f - rad_lim), 3),
            norm_abs
            );
  } else {
    return 0.0;
  }
}


/// 1D E pusher
template<>
void emf::FDTD2_pml<1>::push_e(emf::Tile<1>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C = tile.cfl;
  const auto mins = tile.mins;

  UniIter::iterate(
  [=] DEVCALLABLE (int i, Grids &mesh)
  {
    //-----------
    // global grid coordinates
    float iglob = i + float( mins[0] );

    // dE = dt*curl B  (1/2)
    float dex =  0.0f; // Ex NONE
    float dey =  C*( mesh.bz(i-1,0,0) - mesh.bz(i,0,0));
    float dez =  C*(-mesh.by(i-1,0,0) + mesh.by(i,0,0));

    float lamx = 0.5f*lambda(iglob+0.5f, 0, 0);
    float lamy = 0.5f*lambda(iglob,      0, 0);
    float lamz = 0.5f*lambda(iglob,      0, 0);

    //mesh.ex(i,0,0) = ( mesh.ex(i,0,0)*(1.0f+lamx) + dex )/(1.0f-lamx);
    mesh.ey(i,0,0) = ( mesh.ey(i,0,0)*(1.0f+lamy) + dey )/(1.0f-lamy);
    mesh.ez(i,0,0) = ( mesh.ez(i,0,0)*(1.0f+lamz) + dez )/(1.0f-lamz);

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
void emf::FDTD2_pml<2>::push_e(emf::Tile<2>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C = tile.cfl;
  const auto mins = tile.mins;

  UniIter::iterate2D(
  [=] DEVCALLABLE (int i, int j, Grids &mesh)
  {
    //-----------
    // global grid coordinates
    float iglob = i + float( mins[0] );
    float jglob = j + float( mins[1] );

    // dE = dt*curl B  (1/2)

    float dex = + C*(-mesh.bz(i,  j-1,0) + mesh.bz(i,j,0)); 
    float dey = + C*( mesh.bz(i-1,j,  0) - mesh.bz(i,j,0));
    float dez = + C*( mesh.bx(i,  j-1,0) - mesh.bx(i,j,0))
                  + C*(-mesh.by(i-1,j,  0) + mesh.by(i,j,0));

	  float lamx = 0.5f*lambda(iglob+0.5f, jglob,      0);
	  float lamy = 0.5f*lambda(iglob,      jglob+0.5f, 0);
	  float lamz = 0.5f*lambda(iglob,      jglob,      0);

    mesh.ex(i,j,0) = ( mesh.ex(i,j,0)*(1.0f+lamx) + dex )/(1.0f-lamx);
    mesh.ey(i,j,0) = ( mesh.ey(i,j,0)*(1.0f+lamy) + dey )/(1.0f-lamy);
    mesh.ez(i,j,0) = ( mesh.ez(i,j,0)*(1.0f+lamz) + dez )/(1.0f-lamz);

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
void emf::FDTD2_pml<3>::push_e(emf::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  //const float C = 1.0 * tile.cfl * dt * corr;
  const float C = tile.cfl;
  const auto mins = tile.mins;

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, Grids &mesh)
  {
    //-----------
    // global grid coordinates
    float iglob = i + float( mins[0] );
    float jglob = j + float( mins[1] );
    float kglob = k + float( mins[2] );
    //float lam;

    // dE = dt*curl B  (1/2)
	  //lam = 0.5*lambda(iglob+0.5, jglob, kglob);
    //mesh.ex(i,j,k) = (
    //        mesh.ex(i,j,k)*(1+lam) + 
    //        C*(mesh.by(i,j,k-1)-mesh.by(i,j,k)  - mesh.bz(i,j-1,k)+mesh.bz(i,j,k)) 
    //        )/(1-lam);

	  //lam = 0.5*lambda(iglob, jglob+0.5, kglob);
    //mesh.ey(i,j,k) = (
    //        mesh.ey(i,j,k)*(1+lam) + 
    //        C*(mesh.bz(i-1,j,k)-mesh.bz(i,j,k)  - mesh.bx(i,j,k-1)+mesh.bx(i,j,k)) 
    //        )/(1-lam);

	  //lam = 0.5*lambda(iglob, jglob, kglob+0.5);
    //mesh.ez(i,j,k) = (
    //        mesh.ez(i,j,k)*(1+lam) + 
    //        C*(mesh.bx(i,j-1,k)-mesh.bx(i,j,k)  - mesh.by(i-1,j,k)+mesh.by(i,j,k)) 
    //        )/(1-lam);

    float dex = + C*( mesh.by(i,  j,  k-1) - mesh.by(i,j,k))
                  + C*(-mesh.bz(i,  j-1,k  ) + mesh.bz(i,j,k)); 
    float dey = + C*( mesh.bz(i-1,j,  k  ) - mesh.bz(i,j,k))
                  + C*(-mesh.bx(i,  j,  k-1) + mesh.bx(i,j,k));
    float dez = + C*( mesh.bx(i,  j-1,k  ) - mesh.bx(i,j,k))
                  + C*(-mesh.by(i-1,j,  k  ) + mesh.by(i,j,k));

	  float lamx = 0.5f*lambda(iglob+0.5f, jglob,      kglob);
	  float lamy = 0.5f*lambda(iglob,      jglob+0.5f, kglob);
	  float lamz = 0.5f*lambda(iglob,      jglob,      kglob+0.5f);

    mesh.ex(i,j,k) = ( mesh.ex(i,j,k)*(1.0f+lamx) + dex )/(1.0f-lamx);
    mesh.ey(i,j,k) = ( mesh.ey(i,j,k)*(1.0f+lamy) + dey )/(1.0f-lamy);
    mesh.ez(i,j,k) = ( mesh.ez(i,j,k)*(1.0f+lamz) + dez )/(1.0f-lamz);


    //mesh.ex(i,j,k) += + C*( mesh.by(i,  j,  k-1) - mesh.by(i,j,k))
    //                  + C*(-mesh.bz(i,  j-1,k  ) + mesh.bz(i,j,k)); 
    //mesh.ey(i,j,k) += + C*( mesh.bz(i-1,j,  k  ) - mesh.bz(i,j,k))
    //                  + C*(-mesh.bx(i,  j,  k-1) + mesh.bx(i,j,k));
    //mesh.ez(i,j,k) += + C*( mesh.bx(i,  j-1,k  ) - mesh.bx(i,j,k))
    //                  + C*(-mesh.by(i-1,j,  k  ) + mesh.by(i,j,k));
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
   

/// 1D B pusher
template<>
void emf::FDTD2_pml<1>::push_half_b(emf::Tile<1>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C = 0.5*tile.cfl;
  const auto mins = tile.mins;

  UniIter::iterate(
  [=] DEVCALLABLE (int i, Grids &mesh)
  {

    // global grid coordinates
    float iglob = i + mins[0];

	  float lamx = 0.25f*lambda(iglob,      0, 0);
	  float lamy = 0.25f*lambda(iglob+0.5f, 0, 0);
	  float lamz = 0.25f*lambda(iglob+0.5f, 0, 0);
      
    float dbx = 0.0f; // Bx NONE
    float dby = C*( mesh.ez(i+1,0,0) - mesh.ez(i,0,0));
    float dbz = C*(-mesh.ey(i+1,0,0) + mesh.ey(i,0,0));

	  //mesh.bx(i,0,0) = ( mesh.bx(i,0,0)*(1.0f+lamx) + dbx )/(1.0f-lamx); // no bx update
	  mesh.by(i,0,0) = ( mesh.by(i,0,0)*(1.0f+lamy) + dby )/(1.0f-lamy);
	  mesh.bz(i,0,0) = ( mesh.bz(i,0,0)*(1.0f+lamz) + dbz )/(1.0f-lamz);

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
void emf::FDTD2_pml<2>::push_half_b(emf::Tile<2>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  const float C = 0.5*tile.cfl;
  const auto mins = tile.mins;

  UniIter::iterate2D(
  [=] DEVCALLABLE (int i, int j, Grids &mesh)
  {

    // global grid coordinates
    float iglob = i + mins[0];
    float jglob = j + mins[1];

	  float lamx = 0.25f*lambda(iglob,      jglob+0.5f, 0);
	  float lamy = 0.25f*lambda(iglob+0.5f, jglob,      0);
	  float lamz = 0.25f*lambda(iglob+0.5f, jglob+0.5f, 0);
      
    float dbx = + C*(-mesh.ez(i,  j+1,0) + mesh.ez(i,j,0));
    float dby = + C*( mesh.ez(i+1,j,  0) - mesh.ez(i,j,0));
    float dbz = + C*( mesh.ex(i,  j+1,0) - mesh.ex(i,j,0))
                  + C*(-mesh.ey(i+1,j,  0) + mesh.ey(i,j,0));

	  mesh.bx(i,j,0) = ( mesh.bx(i,j,0)*(1.0f+lamx) + dbx )/(1.0f-lamx);
	  mesh.by(i,j,0) = ( mesh.by(i,j,0)*(1.0f+lamy) + dby )/(1.0f-lamy);
	  mesh.bz(i,j,0) = ( mesh.bz(i,j,0)*(1.0f+lamz) + dbz )/(1.0f-lamz);

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
void emf::FDTD2_pml<3>::push_half_b(emf::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();
  //const float C = 0.5 * tile.cfl * dt * corr;
  const float C = 0.5*tile.cfl;
  const auto mins = tile.mins;

  //for(int k=0; k<tile.mesh_lengths[2]; k++)
  //for(int j=0; j<tile.mesh_lengths[1]; j++)
  //for(int i=0; i<tile.mesh_lengths[0]; i++) {

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, Grids &mesh)
  {

    // global grid coordinates
    float iglob = i + mins[0];
    float jglob = j + mins[1];
    float kglob = k + mins[2];
    //float lam;

    // dB = -dt*curl E  (1/4)
	  //lam = 0.25*lambda(iglob, jglob+0.5, kglob+0.5);
	  //mesh.bx(i,j,k) = ( mesh.bx(i,j,k)*(1+lam)
    //        + C*(mesh.ey(i,j,k+1)-mesh.ey(i,j,k)  - mesh.ez(i,j+1,k)+mesh.ez(i,j,k)) 
    //        )/(1-lam);

	  //lam = 0.25*lambda(iglob+0.5, jglob, kglob+0.5);
	  //mesh.by(i,j,k) = ( mesh.by(i,j,k)*(1+lam)
    //        + C*(mesh.ez(i+1,j,k)-mesh.ez(i,j,k)  - mesh.ex(i,j,k+1)+mesh.ex(i,j,k)) 
    //        )/(1-lam);

	  //lam = 0.25*lambda(iglob+0.5, jglob+0.5, kglob);
	  //mesh.bz(i,j,k) = ( mesh.bz(i,j,k)*(1+lam)
    //        + C*(mesh.ex(i,j+1,k)-mesh.ex(i,j,k)  - mesh.ey(i+1,j,k)+mesh.ey(i,j,k)) 
    //        )/(1-lam);

	  float lamx = 0.25f*lambda(iglob,      jglob+0.5f, kglob+0.5f);
	  float lamy = 0.25f*lambda(iglob+0.5f, jglob,      kglob+0.5f);
	  float lamz = 0.25f*lambda(iglob+0.5f, jglob+0.5f, kglob);
      
    float dbx = + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k))
                  + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));
    float dby = + C*( mesh.ez(i+1,j,  k  ) - mesh.ez(i,j,k))
                  + C*(-mesh.ex(i,  j,  k+1) + mesh.ex(i,j,k));
    float dbz = + C*( mesh.ex(i,  j+1,k  ) - mesh.ex(i,j,k))
                  + C*(-mesh.ey(i+1,j,  k  ) + mesh.ey(i,j,k));

	  mesh.bx(i,j,k) = ( mesh.bx(i,j,k)*(1.0f+lamx) + dbx )/(1.0f-lamx);
	  mesh.by(i,j,k) = ( mesh.by(i,j,k)*(1.0f+lamy) + dby )/(1.0f-lamy);
	  mesh.bz(i,j,k) = ( mesh.bz(i,j,k)*(1.0f+lamz) + dbz )/(1.0f-lamz);

    //mesh.bx(i,j,k) += + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k))
    //                  + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));
    //mesh.by(i,j,k) += + C*( mesh.ez(i+1,j,  k  ) - mesh.ez(i,j,k))
    //                  + C*(-mesh.ex(i,  j,  k+1) + mesh.ex(i,j,k));
    //mesh.bz(i,j,k) += + C*( mesh.ex(i,  j+1,k  ) - mesh.ex(i,j,k))
    //                  + C*(-mesh.ey(i+1,j,  k  ) + mesh.ey(i,j,k));
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
void emf::FDTD2_pml<1>::push_eb(::ffe::Tile<1>& tile)
{
  assert(false);
}

template<>
void emf::FDTD2_pml<2>::push_eb(::ffe::Tile<2>& tile)
{
  assert(false);
}


/// 3D E & B coupled pusher (for FFE)
template<>
void emf::FDTD2_pml<3>::push_eb(::ffe::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // refs to storages
  emf::Grids&     m = tile.get_grids();
  ffe::SlimGrids& dm = tile.dF; 
    
  // Tile limits
  const auto mins = tile.mins;
  //const auto maxs = tile.maxs;

  const float c = tile.cfl; // dt / dx
  const float C1 =  c; // high-order curl operator coefficients

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, 
      ffe::SlimGrids& dm,
      emf::Grids& m
      )
  {
    //-----------
    // global grid coordinates
    float iglob = i + float( mins[0] );
    float jglob = j + float( mins[1] );
    float kglob = k + float( mins[2] );
    float lam;

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


template class emf::FDTD2_pml<1>;
template class emf::FDTD2_pml<2>;
template class emf::FDTD2_pml<3>; 
