#include "rffe2.h"
#include "../../tools/signum.h"
#include "../../em-fields/tile.h"


#include <cmath>

#include "../../tools/iter/iter.h"

//#include <nvtx3/nvToolsExt.h> 

// general bilinear interpolation
//template<>
//void ffe::rFFE<2>::interpolate(
//  toolbox::Mesh<real_short,3>& f,
//  int i, int j,
//  std::array<int,2>& in,
//  std::array<int,2>& out
//    )
//{
//  int im = in[0] == out[0] ? 0 :  -out[0];
//  int ip = in[0] == out[0] ? 0 : 1-out[0];
//  int jm = in[1] == out[1] ? 0 :  -out[1];
//  int jp = in[1] == out[1] ? 0 : 1-out[1];
//
//  real_short f1 = f[i+ip, j+jp, k] + f[i+ip, j+jm, k]
//  real_short f0 = f[i+im, j+jp, k] + f[i+im, j+jm, k] 
//  return 0.25*(f1 + f0);
//}

// general trilinear interpolation
template<>
void ffe::rFFE2<3>::interpolate( 
        toolbox::Mesh<real_short,3>& f,
        toolbox::Mesh<real_short,0>& fi,
        const std::array<int,3>& in,
        const std::array<int,3>& out
      )
{
  //nvtxRangePush(__FUNCTION__);

  int im = in[2] == out[2] ? 0 :  -out[2];
  int ip = in[2] == out[2] ? 0 : 1-out[2];

  int jm = in[1] == out[1] ? 0 :  -out[1];
  int jp = in[1] == out[1] ? 0 : 1-out[1];

  int km = in[0] == out[0] ? 0 :  -out[0];
  int kp = in[0] == out[0] ? 0 : 1-out[0];

  UniIter::iterate3D(
    [=] DEVCALLABLE (int i, int j, int k, toolbox::Mesh<real_short,3>& f, toolbox::Mesh<real_short,0>& fi)
    {
      real_short f11, f10, f01, f00, f1, f0;
      f11 = f(i+ip, j+jp, k+km) + f(i+ip, j+jp, k+kp);
      f10 = f(i+ip, j+jm, k+km) + f(i+ip, j+jm, k+kp);
      f01 = f(i+im, j+jp, k+km) + f(i+im, j+jp, k+kp);
      f00 = f(i+im, j+jm, k+km) + f(i+im, j+jm, k+kp);
      f1  = f11 + f10;
      f0  = f01 + f00;

      fi(i,j,k) = 0.125*(f1 + f0);
    }, {std::tuple{8,8,4}},
    f.Nx, f.Ny, f.Nz, f, fi);
    
  //nvtxRangePop();
}



/// 3D divE = rho
template<>
void ffe::rFFE2<3>::comp_rho(ffe::Tile<3>& tile)
{
  //nvtxRangePush(__FUNCTION__);
  
  fields::YeeLattice& mesh = tile.get_yee();

  // NOTE: compute rho from -1 to +1 because later on we re-stagger it 
  // and need the guard zones for interpolation
  //for(int k=-1; k<static_cast<int>(tile.mesh_lengths[2]+1); k++) {
  //  for(int j=-1; j<static_cast<int>(tile.mesh_lengths[1]+1); j++) {
  //    for(int i=-1; i<static_cast<int>(tile.mesh_lengths[0]+1); i++) {
  //      rho(i,j,k) = 
  //        (ex(i,j,k) - ex(i-1,j,  k  )) +
  //        (ey(i,j,k) - ey(i  ,j-1,k  )) + 
  //        (ez(i,j,k) - ez(i  ,j,  k-1));
  //    }
  //  }
  //}
    
  // TODO check vs above
  UniIter::iterate3D(
    [=] DEVCALLABLE (int i, int j, int k, fields::YeeLattice& m)
    {
        m.rho(i-1,j-1,k-1) = 
          (m.ex(i-1,j-1,k-1) - m.ex(i-1-1,j-1,  k-1  )) +
          (m.ey(i-1,j-1,k-1) - m.ey(i-1  ,j-1-1,k-1  )) + 
          (m.ez(i-1,j-1,k-1) - m.ez(i-1  ,j-1,  k-1-1));
    }, 
    static_cast<int>(tile.mesh_lengths[2])+2,
    static_cast<int>(tile.mesh_lengths[1])+2,
    static_cast<int>(tile.mesh_lengths[0])+2,
    mesh);

  //nvtxRangePop();
  UniIter::sync();

}


/// 3D 
template<>
void ffe::rFFE2<3>::push_eb(ffe::Tile<3>& tile)
{
  //nvtxRangePush(__FUNCTION__);

  // refs to storages
  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  // refs to fields for easier access

  real_short c = tile.cfl;

  // dt / dx
  real_short cx = c;
  real_short cy = c;
  real_short cz = c;

  //for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
  //  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
  //    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

  //      // dB = dt*curl E
  //      dm.bx(i,j,k) = cz*( ey(i,  j,  k+1) - ey(i,j,k) ) - cy*( ez(i,  j+1,k) - ez(i,j,k) );
  //      dm.by(i,j,k) = cx*( ez(i+1,j,  k  ) - ez(i,j,k) ) - cz*( ex(i,  j,k+1) - ex(i,j,k) );
  //      dm.bz(i,j,k) = cy*( ex(i,  j+1,k  ) - ex(i,j,k) ) - cx*( ey(i+1,j,k  ) - ey(i,j,k) );

  //      // dE = dt*curl B 
  //      dm.ex(i,j,k) = cz*( by(i,  j,  k-1) - by(i,j,k) ) - cy*( bz( i,  j-1,k  ) - bz(i,j,k) );
  //      dm.ey(i,j,k) = cx*( bz(i-1,j,  k  ) - bz(i,j,k) ) - cz*( bx( i,  j,  k-1) - bx(i,j,k) );
  //      dm.ez(i,j,k) = cy*( bx(i,  j-1,k  ) - bx(i,j,k) ) - cx*( by( i-1,j,  k  ) - by(i,j,k) );
  //    }
  //  }
  //}

  // TODO check vs above
  UniIter::iterate3D(
    [=] DEVCALLABLE (int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice &m)
    {
      // dB = dt*curl E
      dm.bx(i,j,k) = cz*( m.ey(i,  j,  k+1) - m.ey(i,j,k) ) - cy*( m.ez(i,  j+1,k) - m.ez(i,j,k) );
      dm.by(i,j,k) = cx*( m.ez(i+1,j,  k  ) - m.ez(i,j,k) ) - cz*( m.ex(i,  j,k+1) - m.ex(i,j,k) );
      dm.bz(i,j,k) = cy*( m.ex(i,  j+1,k  ) - m.ex(i,j,k) ) - cx*( m.ey(i+1,j,k  ) - m.ey(i,j,k) );

      // dE = dt*curl B 
      dm.ex(i,j,k) = cz*( m.by(i,  j,  k-1) - m.by(i,j,k) ) - cy*( m.bz( i,  j-1,k  ) - m.bz(i,j,k) );
      dm.ey(i,j,k) = cx*( m.bz(i-1,j,  k  ) - m.bz(i,j,k) ) - cz*( m.bx( i,  j,  k-1) - m.bx(i,j,k) );
      dm.ez(i,j,k) = cy*( m.bx(i,  j-1,k  ) - m.bx(i,j,k) ) - cx*( m.by( i-1,j,  k  ) - m.by(i,j,k) );
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m);
    //nvtxRangePop();
    UniIter::sync();
}



template<>
void ffe::rFFE2<3>::stagger_x_eb(fields::YeeLattice& m)
{
  //nvtxRangePush(__FUNCTION__);

  interpolate(m.ex, exf, {{1,1,0}}, {{1,1,0}} ); //x
  interpolate(m.ey, eyf, {{1,0,1}}, {{1,1,0}} );
  interpolate(m.ez, ezf, {{0,1,1}}, {{1,1,0}} );
  interpolate(m.bx, bxf, {{0,0,1}}, {{1,1,0}} );
  interpolate(m.by, byf, {{0,1,0}}, {{1,1,0}} );
  interpolate(m.bz, bzf, {{1,0,0}}, {{1,1,0}} );
  //nvtxRangePop();

}

template<>
void ffe::rFFE2<3>::stagger_y_eb(fields::YeeLattice& m)
{
  //nvtxRangePush(__FUNCTION__);

  interpolate(m.ex, exf, {{1,1,0}}, {{1,0,1}} );
  interpolate(m.ey, eyf, {{1,0,1}}, {{1,0,1}} ); //y
  interpolate(m.ez, ezf, {{0,1,1}}, {{1,0,1}} );
  interpolate(m.bx, bxf, {{0,0,1}}, {{1,0,1}} );
  interpolate(m.by, byf, {{0,1,0}}, {{1,0,1}} );
  interpolate(m.bz, bzf, {{1,0,0}}, {{1,0,1}} );
  //nvtxRangePop();

}

template<>
void ffe::rFFE2<3>::stagger_z_eb(fields::YeeLattice& m)
{
  //nvtxRangePush(__FUNCTION__);

  interpolate(m.ex, exf, {{1,1,0}}, {{0,1,1}} );
  interpolate(m.ey, eyf, {{1,0,1}}, {{0,1,1}} );
  interpolate(m.ez, ezf, {{0,1,1}}, {{0,1,1}} ); //z
  interpolate(m.bx, bxf, {{0,0,1}}, {{0,1,1}} );
  interpolate(m.by, byf, {{0,1,0}}, {{0,1,1}} );
  interpolate(m.bz, bzf, {{1,0,0}}, {{0,1,1}} );
  //nvtxRangePop();

}


/// 3D 
template<>
void ffe::rFFE2<3>::add_jperp(ffe::Tile<3>& tile)
{
  //nvtxRangePush(__FUNCTION__);

  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  real_short dt = tile.cfl;

  interpolate(m.rho, rhf, { { 1, 1, 1 } }, { { 1, 1, 0 } });
  stagger_x_eb(m);
  UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {

      real_short b2, cur;
      b2 =
        (bxf(i, j, k) * bxf(i, j, k) + byf(i, j, k) * byf(i, j, k) +
         bzf(i, j, k) * bzf(i, j, k) + EPS);

      cur =
        rhf(i, j, k) * (eyf(i, j, k) * bzf(i, j, k) - byf(i, j, k) * ezf(i, j, k)) / b2;
      m.jx(i, j, k) = cur;
      dm.ex(i, j, k) -= dt * cur;
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);

  interpolate(m.rho, rhf, {{1,1,1}}, {{1,0,1}} );
  stagger_y_eb(m);

  UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
        real_short b2, cur;
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);

        cur = rhf(i,j,k) * (ezf(i,j,k)*bxf(i,j,k) - exf(i,j,k)*bzf(i,j,k))/b2;
        m.jy(i,j,k) = cur;
        dm.ey(i,j,k) -= dt*cur;
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);


  interpolate(m.rho, rhf, {{1,1,1}}, {{0,1,1}} );
  stagger_z_eb(m);

    UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
        real_short b2, cur;
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);

        cur = rhf(i,j,k) * (exf(i,j,k)*byf(i,j,k) - bxf(i,j,k)*eyf(i,j,k))/b2;
        m.jz(i,j,k) = cur;
        dm.ez(i,j,k) -= dt*cur;
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);


//nvtxRangePop();
  UniIter::sync();

 }


//TODO move to tile
template<>
void ffe::rFFE2<3>::update_eb(
    ffe::Tile<3>& tile,
    real_short c1, 
    real_short c2, 
    real_short c3
    )
{
  //nvtxRangePush(__FUNCTION__);
  fields::YeeLattice&    m = tile.get_yee();
  ffe::SkinnyYeeLattice& n = tile.Fn; 
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  UniIter::iterate3D(
    [=] DEVCALLABLE(
      int i,
      int j,
      int k,
      ffe::SkinnyYeeLattice& dm,
      fields::YeeLattice& m,
      ffe::SkinnyYeeLattice& n)
      {
        // RK3 E update
        m.ex(i,j,k) = c1*n.ex(i,j,k) + c2*m.ex(i,j,k) + c3*dm.ex(i,j,k);
        m.ey(i,j,k) = c1*n.ey(i,j,k) + c2*m.ey(i,j,k) + c3*dm.ey(i,j,k);
        m.ez(i,j,k) = c1*n.ez(i,j,k) + c2*m.ez(i,j,k) + c3*dm.ez(i,j,k);

        // RK3 B update
        m.bx(i,j,k) = c1*n.bx(i,j,k) + c2*m.bx(i,j,k) + c3*dm.bx(i,j,k);
        m.by(i,j,k) = c1*n.by(i,j,k) + c2*m.by(i,j,k) + c3*dm.by(i,j,k);
        m.bz(i,j,k) = c1*n.bz(i,j,k) + c2*m.bz(i,j,k) + c3*dm.bz(i,j,k);

        // variable switch for 1) e > b and 2) j_par calcs.
        // Enables to calculate both of the above as independent
        // corrections because interpolation is done via m.ex
        // meshes and results are stored in dm.ex meshes:
        dm.ex(i,j,k) = m.ex(i,j,k);
        dm.ey(i,j,k) = m.ey(i,j,k);
        dm.ez(i,j,k) = m.ez(i,j,k);
      },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, n);

//nvtxRangePop();
    UniIter::sync();

  }

template<>
void ffe::rFFE2<3>::remove_jpar(ffe::Tile<3>& tile)
{
UniIter::sync();
  //nvtxRangePush(__FUNCTION__);

  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  real_short cur, b2;
  real_short dt = tile.cfl;

  // NOTE: updates done via dm array to avoid cross contamination between x/y/z diretions

  stagger_x_eb(m);
  
  UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
        real_short cur, b2;
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*bxf(i,j,k) /b2/dt;

        m.jx(i,j,k) += cur;
        //dm.ex(i,j,k) -= cur;
        dm.ex(i,j,k) = m.ex(i,j,k) - cur*dt;
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);

  stagger_y_eb(m);
  
  UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
        real_short cur, b2;
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*byf(i,j,k) /b2/dt;

        m.jy(i,j,k) += cur;
        //dm.ey(i,j,k) -= cur;
        dm.ey(i,j,k) = m.ey(i,j,k) - cur*dt;
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);

  stagger_z_eb(m);
  
  UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
        real_short cur, b2;
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*bzf(i,j,k) /b2/dt;

        m.jz(i,j,k) += cur;
        //dm.ez(i,j,k) -= cur;
        dm.ez(i,j,k) = m.ez(i,j,k) - cur*dt;
    },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);

//nvtxRangePop();
    UniIter::sync();

  //NOTE: only done at the end of the loop so it does not affect previous calculations: dm used as temporary array
  //for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
  //  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
  //    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
  //      m.ex(i,j,k) = dm.ex(i,j,k);
  //      m.ey(i,j,k) = dm.ey(i,j,k);
  //      m.ez(i,j,k) = dm.ez(i,j,k);
  //    }
  //  }
  //}


}


template<>
void ffe::rFFE2<3>::limit_e(ffe::Tile<3>& tile)
{
  //nvtxRangePush(__FUNCTION__);

  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  real_short dt = tile.cfl;


  stagger_x_eb(m);

    UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
        real_short e2, b2, diss, cur;

        e2 = exf(i,j,k)*exf(i,j,k) + eyf(i,j,k)*eyf(i,j,k) + ezf(i,j,k)*ezf(i,j,k) + EPS;
        b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k) + EPS;

        diss = 1.0;
        if (e2 > b2) diss = std::sqrt(b2/e2); 

        cur = (1.-diss)*dm.ex(i,j,k)/dt;
        m.jx(i,j,k) += cur;
        //dm.ex(i,j,k) = diss*m.ex(i,j,k);
        m.ex(i,j,k) = diss*dm.ex(i,j,k);
  },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);

  stagger_y_eb(m);
    UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {
      
        real_short e2, b2, diss, cur;

        e2 = exf(i,j,k)*exf(i,j,k) + eyf(i,j,k)*eyf(i,j,k) + ezf(i,j,k)*ezf(i,j,k) + EPS;
        b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k) + EPS;

        diss = 1.0;
        if (e2 > b2) diss = std::sqrt(b2/e2);

        cur = (1.-diss)*dm.ey(i,j,k)/dt;
        m.jy(i,j,k) += cur;
        //dm.ey(i,j,k) = diss*m.ey(i,j,k);
        m.ey(i,j,k) = diss*dm.ey(i,j,k);

  },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);

  stagger_z_eb(m);
    UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, ffe::SkinnyYeeLattice& dm, fields::YeeLattice& m, 
    toolbox::Mesh<real_short, 0>& bxf,
    toolbox::Mesh<real_short, 0>& byf, toolbox::Mesh<real_short, 0>& bzf,
    toolbox::Mesh<real_short, 0>& exf, toolbox::Mesh<real_short, 0>& eyf,
    toolbox::Mesh<real_short, 0>& ezf, toolbox::Mesh<real_short, 0>& rhf) {

        real_short e2, b2, diss, cur;
        
        e2 = exf(i,j,k)*exf(i,j,k) + eyf(i,j,k)*eyf(i,j,k) + ezf(i,j,k)*ezf(i,j,k) + EPS;
        b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k) + EPS;

        diss = 1.0;
        if (e2 > b2) diss = std::sqrt(b2/e2);

        cur = (1.-diss)*dm.ez(i,j,k)/dt;
        m.jz(i,j,k) += cur;
        //dm.ez(i,j,k) = diss*m.ez(i,j,k);
        //m.ez(i,j,k) -= cur*dt;
        m.ez(i,j,k) = diss*dm.ez(i,j,k); //new
  },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    dm, m, bxf, byf, bzf, exf, eyf, ezf, rhf);


//nvtxRangePop();
  UniIter::sync();

  //NOTE: only done at the end of the loop so it does not affect previous calculations: dm used as temporary array
  //for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
  //  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
  //    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
  //      m.ex(i,j,k) = dm.ex(i,j,k);
  //      m.ey(i,j,k) = dm.ey(i,j,k);
  //      m.ez(i,j,k) = dm.ez(i,j,k);
  //    }
  //  }
  //}

}


//TODO move to tile
template<>
void ffe::rFFE2<3>::copy_eb( ffe::Tile<3>& tile)
{
  //nvtxRangePush(__FUNCTION__);
  fields::YeeLattice&    m = tile.get_yee();
  ffe::SkinnyYeeLattice& n = tile.Fn; 
  //ffe::SkinnyYeeLattice& dm = tile.dF; 

      UniIter::iterate3D(
    [=] DEVCALLABLE( int i, int j, int k, fields::YeeLattice&    m, ffe::SkinnyYeeLattice& n)
    {
        n.ex(i,j,k) = m.ex(i,j,k);
        n.ey(i,j,k) = m.ey(i,j,k);
        n.ez(i,j,k) = m.ez(i,j,k);

        n.bx(i,j,k) = m.bx(i,j,k);
        n.by(i,j,k) = m.by(i,j,k);
        n.bz(i,j,k) = m.bz(i,j,k);
  },
    static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]),
    m, n);
//nvtxRangePop();
  UniIter::sync();
}


//--------------------------------------------------
// explicit template instantiation
//template class ffe::rFFE2<2>; // 2D
template class ffe::rFFE2<3>; // 3D
