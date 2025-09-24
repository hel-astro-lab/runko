#include <cmath>

#include "core/ffe/currents/rffe4.h"
#include "tools/signum.h"
#include "core/emf/tile.h"



// general trilinear interpolation
template<>
void ffe::rFFE4<3>::interpolate( 
        toolbox::Mesh<float,3>& f,
        toolbox::Mesh<float,0>& fi,
        const std::array<int,3>& in,
        const std::array<int,3>& out
      )
{
  int im = in[2] == out[2] ? 0 :  -out[2];
  int ip = in[2] == out[2] ? 0 : 1-out[2];

  int jm = in[1] == out[1] ? 0 :  -out[1];
  int jp = in[1] == out[1] ? 0 : 1-out[1];

  int km = in[0] == out[0] ? 0 :  -out[0];
  int kp = in[0] == out[0] ? 0 : 1-out[0];

  float f11, f10, f01, f00, f1, f0;

  for(int k=0; k<f.Nz; k++) {
    for(int j=0; j<f.Ny; j++) {
      for(int i=0; i<f.Nx; i++) {
        f11 = f(i+ip, j+jp, k+km) + f(i+ip, j+jp, k+kp);
        f10 = f(i+ip, j+jm, k+km) + f(i+ip, j+jm, k+kp);
        f01 = f(i+im, j+jp, k+km) + f(i+im, j+jp, k+kp);
        f00 = f(i+im, j+jm, k+km) + f(i+im, j+jm, k+kp);
        f1  = f11 + f10;
        f0  = f01 + f00;

        fi(i,j,k) = 0.125*(f1 + f0);
      }
    }
  }
}

template<>
void ffe::rFFE4<3>::stagger_x_eb(emf::Grids& m)
{
  interpolate(m.ex, exf, {{1,1,0}}, {{1,1,0}} ); //x
  interpolate(m.ey, eyf, {{1,0,1}}, {{1,1,0}} );
  interpolate(m.ez, ezf, {{0,1,1}}, {{1,1,0}} );
  interpolate(m.bx, bxf, {{0,0,1}}, {{1,1,0}} );
  interpolate(m.by, byf, {{0,1,0}}, {{1,1,0}} );
  interpolate(m.bz, bzf, {{1,0,0}}, {{1,1,0}} );
}

template<>
void ffe::rFFE4<3>::stagger_y_eb(emf::Grids& m)
{
  interpolate(m.ex, exf, {{1,1,0}}, {{1,0,1}} );
  interpolate(m.ey, eyf, {{1,0,1}}, {{1,0,1}} ); //y
  interpolate(m.ez, ezf, {{0,1,1}}, {{1,0,1}} );
  interpolate(m.bx, bxf, {{0,0,1}}, {{1,0,1}} );
  interpolate(m.by, byf, {{0,1,0}}, {{1,0,1}} );
  interpolate(m.bz, bzf, {{1,0,0}}, {{1,0,1}} );
}

template<>
void ffe::rFFE4<3>::stagger_z_eb(emf::Grids& m)
{
  interpolate(m.ex, exf, {{1,1,0}}, {{0,1,1}} );
  interpolate(m.ey, eyf, {{1,0,1}}, {{0,1,1}} );
  interpolate(m.ez, ezf, {{0,1,1}}, {{0,1,1}} ); //z
  interpolate(m.bx, bxf, {{0,0,1}}, {{0,1,1}} );
  interpolate(m.by, byf, {{0,1,0}}, {{0,1,1}} );
  interpolate(m.bz, bzf, {{1,0,0}}, {{0,1,1}} );
}



/// 3D divE = rho
template<>
void ffe::rFFE4<3>::comp_rho(ffe::Tile<3>& tile)
{
  emf::Grids& mesh = tile.get_grids();
  auto& rho = mesh.rho;
  auto& ex  = mesh.ex;
  auto& ey  = mesh.ey;
  auto& ez  = mesh.ez;

  // high-order curl operator coefficients
  float C1 = 9.0/8.0;
  float C2 = 1.0/24.0;

  // NOTE: compute rho from -1 to +1 because later on we re-stagger it 
  // and need the guard zones for interpolation
  for(int k=-1; k<static_cast<int>(tile.mesh_lengths[2]+1); k++) {
    for(int j=-1; j<static_cast<int>(tile.mesh_lengths[1]+1); j++) {
      for(int i=-1; i<static_cast<int>(tile.mesh_lengths[0]+1); i++) {

        rho(i,j,k) = 
          C1*(ex(i,j,k) - ex(i-1,j,k)) + C2*(ex(i-2,j,k) - ex(i+1,j,k)) +
          C1*(ey(i,j,k) - ey(i,j-1,k)) + C2*(ey(i,j-2,k) - ey(i,j+1,k)) +
          C1*(ez(i,j,k) - ez(i,j,k-1)) + C2*(ez(i,j,k-2) - ez(i,j,k+1));
      }
    }
  }
  }

/// 3D 
template<>
void ffe::rFFE4<3>::push_eb(ffe::Tile<3>& tile)
{
  // refs to storages
  emf::Grids&     m = tile.get_grids();
  ffe::SlimGrids& dm = tile.dF; 

  // refs to emf for easier access
  auto& ex  = m.ex;
  auto& ey  = m.ey;
  auto& ez  = m.ez;

  auto& bx  = m.bx;
  auto& by  = m.by;
  auto& bz  = m.bz;

  float c = tile.cfl;

  // dt / dx
  // high-order curl operator coefficients
  float C1 = c*9.0/8.0;
  float C2 = c*1.0/24.0;

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        // dB = dt*curl E
		dm.bx(i,j,k) =
            C1*(ey(i,j,k+1)-ey(i,j,k)  - ez(i,j+1,k)+ez(i,j,k)) -
            C2*(ey(i,j,k+2)-ey(i,j,k-1)- ez(i,j+2,k)+ez(i,j-1,k));

		dm.by(i,j,k) =
            C1*(ez(i+1,j,k)-ez(i,j,k)  - ex(i,j,k+1)+ex(i,j,k)) -
            C2*(ez(i+2,j,k)-ez(i-1,j,k)- ex(i,j,k+2)+ex(i,j,k-1));

		dm.bz(i,j,k) =
            C1*(ex(i,j+1,k)-ex(i,j,k)  - ey(i+1,j,k)+ey(i,j,k)) -
            C2*(ex(i,j+2,k)-ex(i,j-1,k)- ey(i+2,j,k)+ey(i-1,j,k));

        // dE = dt*curl B 
        dm.ex(i,j,k) =
          C1*(by(i,j,k-1)-by(i,j,k)  - bz(i,j-1,k)+bz(i,j,k)) -
          C2*(by(i,j,k-2)-by(i,j,k+1)- bz(i,j-2,k)+bz(i,j+1,k));

        dm.ey(i,j,k) =
          C1*(bz(i-1,j,k)-bz(i,j,k)  - bx(i,j,k-1)+bx(i,j,k)) -
          C2*(bz(i-2,j,k)-bz(i+1,j,k)- bx(i,j,k-2)+bx(i,j,k+1));

        dm.ez(i,j,k) =
          C1*(bx(i,j-1,k)-bx(i,j,k)  - by(i-1,j,k)+by(i,j,k)) -
          C2*(bx(i,j-2,k)-bx(i,j+1,k)- by(i-2,j,k)+by(i+1,j,k));

      }
    }
  }

}


/// 3D 
template<>
void ffe::rFFE4<3>::add_jperp(ffe::Tile<3>& tile)
{
  emf::Grids&     m = tile.get_grids();
  ffe::SlimGrids& dm = tile.dF; 

  auto& jx  = m.jx;
  auto& jy  = m.jy;
  auto& jz  = m.jz;

  float dt = tile.cfl;
  float b2, cur;

  interpolate(m.rho, rhf, {{1,1,1}}, {{1,1,0}} );
  stagger_x_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        b2 = (
            bxf(i,j,k)*bxf(i,j,k) 
          + byf(i,j,k)*byf(i,j,k) 
          + bzf(i,j,k)*bzf(i,j,k) 
          + EPS);

        cur = rhf(i,j,k) * (eyf(i,j,k)*bzf(i,j,k) - byf(i,j,k)*ezf(i,j,k) )/b2;
        jx(i,j,k) = cur;
        dm.ex(i,j,k) -= dt*cur;
      }
    }
  }

  interpolate(m.rho, rhf, {{1,1,1}}, {{1,0,1}} );
  stagger_y_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
          
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);

        cur = rhf(i,j,k) * (ezf(i,j,k)*bxf(i,j,k) - exf(i,j,k)*bzf(i,j,k))/b2;
        jy(i,j,k) = cur;
        dm.ey(i,j,k) -= dt*cur;
      }
    }
  }

  interpolate(m.rho, rhf, {{1,1,1}}, {{0,1,1}} );
  stagger_z_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);

        cur = rhf(i,j,k) * (exf(i,j,k)*byf(i,j,k) - bxf(i,j,k)*eyf(i,j,k))/b2;
        jz(i,j,k) = cur;
        dm.ez(i,j,k) -= dt*cur;
      }
    }
  }

 }



template<>
void ffe::rFFE4<3>::remove_jpar(ffe::Tile<3>& tile)
{
  emf::Grids&     m = tile.get_grids();
  ffe::SlimGrids& dm = tile.dF; 

  float cur, b2;
  float dt = tile.cfl;

  // NOTE: updates done via dm array to avoid cross contamination between x/y/z diretions

  stagger_x_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*bxf(i,j,k) /b2/dt;

        m.jx(i,j,k) += cur;
        //dm.ex(i,j,k) -= cur;
        dm.ex(i,j,k) = m.ex(i,j,k) - cur*dt;
      }
    }
  }

  stagger_y_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*byf(i,j,k) /b2/dt;

        m.jy(i,j,k) += cur;
        //dm.ey(i,j,k) -= cur;
        dm.ey(i,j,k) = m.ey(i,j,k) - cur*dt;
      }
    }
  }

  stagger_z_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*bzf(i,j,k) /b2/dt;

        m.jz(i,j,k) += cur;
        //dm.ez(i,j,k) -= cur;
        dm.ez(i,j,k) = m.ez(i,j,k) - cur*dt;
      }
    }
  }


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
void ffe::rFFE4<3>::limit_e(ffe::Tile<3>& tile)
{
  emf::Grids&     m = tile.get_grids();
  ffe::SlimGrids& dm = tile.dF; 

  float dt = tile.cfl;
  float e2, b2, diss, cur;


  stagger_x_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        e2 = exf(i,j,k)*exf(i,j,k) + eyf(i,j,k)*eyf(i,j,k) + ezf(i,j,k)*ezf(i,j,k) + EPS;
        b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k) + EPS;

        diss = 1.0;
        if (e2 > b2) diss = std::sqrt(b2/e2); 

        cur = (1.-diss)*m.ex(i,j,k)/dt;
        m.jx(i,j,k) += cur;

        //dm.ex(i,j,k) = diss*m.ex(i,j,k); old
        m.ex(i,j,k) = diss*dm.ex(i,j,k); // new
      }
    }
  }

  stagger_y_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        e2 = exf(i,j,k)*exf(i,j,k) + eyf(i,j,k)*eyf(i,j,k) + ezf(i,j,k)*ezf(i,j,k) + EPS;
        b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k) + EPS;

        diss = 1.0;
        if (e2 > b2) diss = std::sqrt(b2/e2);

        cur = (1.-diss)*m.ey(i,j,k)/dt;
        m.jy(i,j,k) += cur;
        //dm.ey(i,j,k) = diss*m.ey(i,j,k);
        m.ey(i,j,k) = diss*dm.ey(i,j,k); //new
      }
    }
  }

  stagger_z_eb(m);

  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

        e2 = exf(i,j,k)*exf(i,j,k) + eyf(i,j,k)*eyf(i,j,k) + ezf(i,j,k)*ezf(i,j,k) + EPS;
        b2 = bxf(i,j,k)*bxf(i,j,k) + byf(i,j,k)*byf(i,j,k) + bzf(i,j,k)*bzf(i,j,k) + EPS;

        diss = 1.0;
        if (e2 > b2) diss = std::sqrt(b2/e2);

        cur = (1.-diss)*m.ez(i,j,k)/dt;
        m.jz(i,j,k) += cur;
        //dm.ez(i,j,k) = diss*m.ez(i,j,k);
        m.ez(i,j,k) = diss*dm.ez(i,j,k); //new
      }
    }
  }

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




//--------------------------------------------------
// explicit template instantiation
template class ffe::rFFE4<3>; // 3D
