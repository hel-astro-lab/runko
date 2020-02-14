#include "rffe.h"
#include "../../tools/signum.h"
#include "../../em-fields/tile.h"

#include <cmath>


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
  std::array<int,3> in,
  std::array<int,3> out)
{
  int im = in[0] == out[0] ? 0 :  -out[0];
  int ip = in[0] == out[0] ? 0 : 1-out[0];
  int jm = in[1] == out[1] ? 0 :  -out[1];
  int jp = in[1] == out[1] ? 0 : 1-out[1];
  int km = in[2] == out[2] ? 0 :  -out[2];
  int kp = in[2] == out[2] ? 0 : 1-out[2];

  real_short f11, f10, f01, f00, f1, f0;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        f11 = f[i+ip, j+jp, k+km] + f[i+ip, j+jp, k+kp];
        f10 = f[i+ip, j+jm, k+km] + f[i+ip, j+jm, k+kp];
        f01 = f[i+im, j+jp, k+km] + f[i+im, j+jp, k+kp];
        f00 = f[i+im, j+jm, k+km] + f[i+im, j+jm, k+kp];
        f1  = f11 + f10;
        f0  = f01 + f00;

        fi(i,j,k) = 0.125*(f1 + f0);
      }
    }
  }
}



/// 3D divE = rho
template<>
void ffe::rFFE2<3>::comp_rho(ffe::Tile<3>& tile)
{
  fields::YeeLattice& mesh = tile.get_yee();
  auto& rho = mesh.rho;
  auto& ex  = mesh.ex;
  auto& ey  = mesh.ey;
  auto& ez  = mesh.ez;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        rho(i,j,k) = 
          (ex(i,j,k) - ex(i-1,j,  k  )) +
          (ey(i,j,k) - ey(i  ,j-1,k  )) + 
          (ez(i,j,k) - ez(i  ,j,  k-1));
      }
    }
  }
  return;
}

/// 3D 
template<>
void ffe::rFFE2<3>::push_eb(ffe::Tile<3>& tile)
{
  fields::YeeLattice& mesh = tile.get_yee();
  real_short c = tile.cfl;

  // dt / dx
  real_short cx = c;
  real_short cy = c;
  real_short cz = c;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        // dB = curl E
        dbx(i,j,k) = cz*( ey(i,  j,  k+1) - ey(i,j,k) ) - cy*( ez(i,  j+1,k) - ez(i,j,k) );
        dby(i,j,k) = cx*( ez(i+1,j,  k  ) - ez(i,j,k) ) - cz*( ex(i,  j,k+1) - ex(i,j,k) );
        dbz(i,j,k) = cy*( ex(i,  j+1,k  ) - ex(i,j,k) ) - cx*( ey(i+1,j,k  ) - ey(i,j,k) );

        // dE = curl B - curl B_0
        dex(i,j,k) = cz*( by(i,  j,  k-1) - by(i,j,k) ) - cy*( bz( i,  j-1,k  ) - bz(i,j,k) );
        dey(i,j,k) = cx*( bz(i-1,j,  k  ) - bz(i,j,k) ) - cz*( bx( i,  j,  k-1) - bx(i,j,k) );
        dez(i,j,k) = cy*( bx(i,  j-1,k  ) - bx(i,j,k) ) - cx*( by( i-1,j,  k  ) - by(i,j,k) );
      }
    }
  }

  return;
}


template<>
void ffe::rFFE2<3>::stagger_x_eb(fields::YeeLattice& m)
{
  interpolate(m.ex, this->exf, {{1,1,0}}, {{1,1,0}} );
  interpolate(m.ey, this->eyf, {{1,0,1}}, {{1,1,0}} );
  interpolate(m.ez, this->ezf, {{0,1,1}}, {{1,1,0}} );
  interpolate(m.bx, this->bxf, {{0,0,1}}, {{1,1,0}} );
  interpolate(m.by, this->byf, {{0,1,0}}, {{1,1,0}} );
  interpolate(m.bz, this->bzf, {{1,0,0}}, {{1,1,0}} );
}

template<>
void ffe::rFFE2<3>::stagger_y_eb(fields::YeeLattice& m)
{
  interpolate(m.ex, this->exf, {{1,1,0}}, {{1,0,1}} );
  interpolate(m.ey, this->eyf, {{1,0,1}}, {{1,0,1}} );
  interpolate(m.ez, this->ezf, {{0,1,1}}, {{1,0,1}} );
  interpolate(m.bx, this->bxf, {{0,0,1}}, {{1,0,1}} );
  interpolate(m.by, this->byf, {{0,1,0}}, {{1,0,1}} );
  interpolate(m.bz, this->bzf, {{1,0,0}}, {{1,0,1}} );
}

template<>
void ffe::rFFE2<3>::stagger_z_eb(fields::YeeLattice& m)
{
  interpolate(m.ex, this->exf, {{1,1,0}}, {{0,1,1}} );
  interpolate(m.ey, this->eyf, {{1,0,1}}, {{0,1,1}} );
  interpolate(m.ez, this->ezf, {{0,1,1}}, {{0,1,1}} );
  interpolate(m.bx, this->bxf, {{0,0,1}}, {{0,1,1}} );
  interpolate(m.by, this->byf, {{0,1,0}}, {{0,1,1}} );
  interpolate(m.bz, this->bzf, {{1,0,0}}, {{0,1,1}} );
}


/// 3D 
template<>
void ffe::rFFE2<3>::add_jperp(ffe::Tile<3>& tile)
{
  fields::YeeLattice& m = tile.get_yee();

  interpolate(m.rh, this->rhf, {{1,1,1}}, {{1,1,0}} );
  stagger_x_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        b2 = (bxf*bxf + byf*byf + bzf*bzf + EPS);

        jx(i,j,k) = dt *rhf * (eyf * bzf - byf*ezf)/b2;
        dex(ijk) -= jx(i,j,k);
      }
    }
  }

  interpolate(m.rh, this->rhf, {{1,1,1}}, {{1,0,1}} );
  stagger_y_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        b2 = (bxf*bxf + byf*byf + bzf*bzf + EPS);

        jy(i,j,k) = dt *rhf * (ezf*bxf - exf*bzf)/b2;
        dey(ijk) -= jy(i,j,k);
      }
    }
  }

  interpolate(m.rh, this->rhf, {{1,1,1}}, {{0,1,1}} );
  stagger_z_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        b2 = (bxf*bxf + byf*byf + bzf*bzf + EPS);

        jz(i,j,k) = dt *rhf * (exf*byf - bxf*eyf)/b2;
        dez(ijk) -= jz(i,j,k);
      }
    }
  }

 return;
}


template<>
void ffe::rFFE2<3>::update_eb(
    ffe::Tile<3>& tile,
    real_short c1, real_short c2, real_short c3
    )
{
  fields::YeeLattice& m = tile.get_yee();
  //fields::YeeLattice& n = tile.get_skinny_yee();

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        // RK3 E update
        m.ex(i,j,k) = c1 * n.enx(i,j,k) + c2 * m.ex(i,j,k) + c3 * n.dex(i,j,k);
        m.ey(i,j,k) = c1 * n.eny(i,j,k) + c2 * m.ey(i,j,k) + c3 * n.dey(i,j,k);
        m.ez(i,j,k) = c1 * n.enz(i,j,k) + c2 * m.ez(i,j,k) + c3 * n.dez(i,j,k);

        // active variable switch 
        n.dex(i,j,k) = m.ex(i,j,k);
        n.dey(i,j,k) = m.ey(i,j,k);
        n.dez(i,j,k) = m.ez(i,j,k);

        // RK3 B update
        m.bx(i,j,k) = c1 * n.bnx(i,j,k) + c2 * m.bx(i,j,k) + c3 * n.dbx(i,j,k);
        m.by(i,j,k) = c1 * n.bny(i,j,k) + c2 * m.by(i,j,k) + c3 * n.dby(i,j,k);
        m.bz(i,j,k) = c1 * n.bnz(i,j,k) + c2 * m.bz(i,j,k) + c3 * n.dbz(i,j,k);
      }
    }
  }

  return;
}


template<>
void ffe::rFFE2<3>::remove_epar(ffe::Tile<3>& tile)
{
  fields::YeeLattice& m = tile.get_yee();

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        iex = interpolate(ex, i,j,k, {{1,1,0}}, {{1,1,0}} );
        iey = interpolate(ey, i,j,k, {{1,0,1}}, {{1,1,0}} );
        iez = interpolate(ez, i,j,k, {{0,1,1}}, {{1,1,0}} );
        ibx = interpolate(bx, i,j,k, {{0,0,1}}, {{1,1,0}} );
        iby = interpolate(by, i,j,k, {{0,1,0}}, {{1,1,0}} );
        ibz = interpolate(bz, i,j,k, {{1,0,0}}, {{1,1,0}} );

        b2 = ibx*ibx + iby*iby + ibz*ibz + EPS;
        dex(i,j,k) = ex(i,j,k) - (iex*ibx + iey*iby + iex*ibz)*ibx /b2;


        iex = interpolate(ex, i,j,k, {{1,1,0}}, {{1,0,1}} );
        iey = interpolate(ey, i,j,k, {{1,0,1}}, {{1,0,1}} );
        iez = interpolate(ez, i,j,k, {{0,1,1}}, {{1,0,1}} );
        ibx = interpolate(bx, i,j,k, {{0,0,1}}, {{1,0,1}} );
        iby = interpolate(by, i,j,k, {{0,1,0}}, {{1,0,1}} );
        ibz = interpolate(bz, i,j,k, {{1,0,0}}, {{1,0,1}} );

        b2 = ibx*ibx + iby*iby + ibz*ibz + EPS;
        dey(i,j,k) = ey(i,j,k) - (iex*ibx + iey*iby + iex*ibz)*iby /b2;


        iex = interpolate(ex, i,j,k, {{1,1,0}}, {{0,1,1}} );
        iey = interpolate(ey, i,j,k, {{1,0,1}}, {{0,1,1}} );
        iez = interpolate(ez, i,j,k, {{0,1,1}}, {{0,1,1}} );
        ibx = interpolate(bx, i,j,k, {{0,0,1}}, {{0,1,1}} );
        iby = interpolate(by, i,j,k, {{0,1,0}}, {{0,1,1}} );
        ibz = interpolate(bz, i,j,k, {{1,0,0}}, {{0,1,1}} );

        b2 = ibx*ibx + iby*iby + ibz*ibz + EPS;
        dez(i,j,k) = ez(i,j,k) - (iex*ibx + iey*iby + iex*ibz)*ibz /b2;

      }
    }
  }


  return;
}


template<>
void ffe::rFFE2<3>::limit_e(ffe::Tile<3>& tile)
{
  fields::YeeLattice& m = tile.get_yee();

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        iex = interpolate(ex, i,j,k, {{1,1,0}}, {{1,1,0}} );
        iey = interpolate(ey, i,j,k, {{1,0,1}}, {{1,1,0}} );
        iez = interpolate(ez, i,j,k, {{0,1,1}}, {{1,1,0}} );
        ibx = interpolate(bx, i,j,k, {{0,0,1}}, {{1,1,0}} );
        iby = interpolate(by, i,j,k, {{0,1,0}}, {{1,1,0}} );
        ibz = interpolate(bz, i,j,k, {{1,0,0}}, {{1,1,0}} );

        e2 = iex*iex + iey*iey + iez*iez;
        b2 = ibx*ibx + iby*iby + ibz*ibz + EPS;

        diss = 1.0
        if (e2 > b2) diss = sqrt(b2/e2);
        // TODO: dex or ex
        ex(i,j,k) = diss*dex(i,j,k);

        iex = interpolate(ex, i,j,k, {{1,1,0}}, {{1,0,1}} );
        iey = interpolate(ey, i,j,k, {{1,0,1}}, {{1,0,1}} );
        iez = interpolate(ez, i,j,k, {{0,1,1}}, {{1,0,1}} );
        ibx = interpolate(bx, i,j,k, {{0,0,1}}, {{1,0,1}} );
        iby = interpolate(by, i,j,k, {{0,1,0}}, {{1,0,1}} );
        ibz = interpolate(bz, i,j,k, {{1,0,0}}, {{1,0,1}} );

        e2 = iex*iex + iey*iey + iez*iez;
        b2 = ibx*ibx + iby*iby + ibz*ibz + EPS;
        if (e2 > b2) diss = sqrt(b2/e2);
        ey(i,j,k) = diss*dey(i,j,k);


        iex = interpolate(ex, i,j,k, {{1,1,0}}, {{0,1,1}} );
        iey = interpolate(ey, i,j,k, {{1,0,1}}, {{0,1,1}} );
        iez = interpolate(ez, i,j,k, {{0,1,1}}, {{0,1,1}} );
        ibx = interpolate(bx, i,j,k, {{0,0,1}}, {{0,1,1}} );
        iby = interpolate(by, i,j,k, {{0,1,0}}, {{0,1,1}} );
        ibz = interpolate(bz, i,j,k, {{1,0,0}}, {{0,1,1}} );

        e2 = iex*iex + iey*iey + iez*iez;
        b2 = ibx*ibx + iby*iby + ibz*ibz + EPS;
        if (e2 > b2) diss = sqrt(b2/e2);
        ez(i,j,k) = diss*dez(i,j,k);


      }
    }
  }


  return;
}
