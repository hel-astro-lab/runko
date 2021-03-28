#include "fdtd_general.h"

#include <cmath>


inline float_m Dm_x( toolbox::Mesh<float_m, 3>& f, int i, int j, int k, int ai, int /*bi*/, int bj, int bk) {
    return f(i+ai-1, j+bj, k+bk) -f(i-ai, j+bj, k+bk);
}

inline float_m Dm_y( toolbox::Mesh<float_m, 3>& f, int i, int j, int k, int ai, int bi, int /*bj*/, int bk) {
    return f(i+bi, j+ai-1, k+bk) -f(i+bi, j-ai, k+bk);
}

inline float_m Dm_z( toolbox::Mesh<float_m, 3>& f, int i, int j, int k, int ai, int bi, int bj, int /*bk*/) {
    return f(i+bi, j+bj, k+ai-1) -f(i+bi, j+bj, k-ai);
}

//-------------------------------------------------- 
inline float_m Dp_x( toolbox::Mesh<float_m, 3>& f, int i, int j, int k, int ai, int /*bi*/, int bj, int bk) {
    return f(i+ai, j+bj, k+bk) - f(i-ai+1, j+bj, k+bk);
}

inline float_m Dp_y( toolbox::Mesh<float_m, 3>& f, int i, int j, int k, int ai, int bi, int /*bj*/, int bk) {
    return f(i+bi, j+ai, k+bk) - f(i+bi, j-ai+1, k+bk);
}

inline float_m Dp_z( toolbox::Mesh<float_m, 3>& f, int i, int j, int k, int ai, int bi, int bj, int /*bk*/) {
    return f(i+bi, j+bj, k+ai) - f(i+bi, j+bj, k-ai+1);
}


/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */

/// 3D E pusher
template<>
void fields::FDTDGen<3>::push_e(fields::Tile<3>& tile)
{
  YeeLattice& mesh = tile.get_yee();
  
  float_m Cx, Cy, Cz;

  const int Nx = tile.mesh_lengths[0];
  const int Ny = tile.mesh_lengths[1];
  const int Nz = tile.mesh_lengths[2];

  // dE/dt = +curlB
  for(int ai=1; ai<=3; ai++) { //alphas
    Cx = CXs(ai, 0, 0)*corr*tile.cfl;
    Cy = CYs(0, ai, 0)*corr*tile.cfl;
    Cz = CZs(0,  0,ai)*corr*tile.cfl;

    for(int k=0; k<Nz; k++) 
    for(int j=0; j<Ny; j++) 
    for(int i=0; i<Nx; i++){
        //legs with 0-offset
        mesh.ex(i,j,k) += +Cy*Dm_y(mesh.bz,i,j,k, ai, 0,0,0); 
        mesh.ex(i,j,k) += -Cz*Dm_z(mesh.by,i,j,k, ai, 0,0,0);
        mesh.ey(i,j,k) += -Cx*Dm_x(mesh.bz,i,j,k, ai, 0,0,0);
        mesh.ey(i,j,k) += +Cz*Dm_z(mesh.bx,i,j,k, ai, 0,0,0);
        mesh.ez(i,j,k) += +Cx*Dm_x(mesh.by,i,j,k, ai, 0,0,0);
        mesh.ez(i,j,k) += -Cy*Dm_y(mesh.bx,i,j,k, ai, 0,0,0);
    }

    //legs with n-offset
    //for(int bi : {-1,+1}) {
    for(int bi : {-3,-2,-1,+1,+2,+3}) {

        // hand-coded beta fetching from array
        Cx = CXs(ai, 1, 1)*corr*tile.cfl;
        Cy = CYs(1, ai, 1)*corr*tile.cfl;
        Cz = CZs(1,  1,ai)*corr*tile.cfl;

        for(int k=0; k<Nz; k++) 
        for(int j=0; j<Ny; j++) 
        for(int i=0; i<Nx; i++){
            // curl_x
            mesh.ex(i,j,k) += +Cy*Dm_y(mesh.bz,i,j,k, ai, bi,0,0); //beta_yx
            mesh.ex(i,j,k) += +Cy*Dm_y(mesh.bz,i,j,k, ai, 0,0,bi); //beta_yz

            mesh.ex(i,j,k) += -Cz*Dm_z(mesh.by,i,j,k, ai, bi,0,0);
            mesh.ex(i,j,k) += -Cz*Dm_z(mesh.by,i,j,k, ai, 0,bi,0);

            // curl_y
            mesh.ey(i,j,k) += -Cx*Dm_x(mesh.bz,i,j,k, ai, 0,bi,0);
            mesh.ey(i,j,k) += -Cx*Dm_x(mesh.bz,i,j,k, ai, 0,0,bi);

            mesh.ey(i,j,k) += +Cz*Dm_z(mesh.bx,i,j,k, ai, bi,0,0);
            mesh.ey(i,j,k) += +Cz*Dm_z(mesh.bx,i,j,k, ai, 0,bi,0);

            // curl_z
            mesh.ez(i,j,k) += +Cx*Dm_x(mesh.by,i,j,k, ai, 0,bi,0);
            mesh.ez(i,j,k) += +Cx*Dm_x(mesh.by,i,j,k, ai, 0,0,bi);

            mesh.ez(i,j,k) += -Cy*Dm_y(mesh.bx,i,j,k, ai, bi,0,0);
            mesh.ez(i,j,k) += -Cy*Dm_y(mesh.bx,i,j,k, ai, 0,0,bi);
        }
    }
  }

}


//--------------------------------------------------
/// Update B field with a half step

/// 3D B pusher
template<>
void fields::FDTDGen<3>::push_half_b(fields::Tile<3>& tile)
{
  YeeLattice& mesh = tile.get_yee();

  float_m Cx, Cy, Cz;

  const int Nx = tile.mesh_lengths[0];
  const int Ny = tile.mesh_lengths[1];
  const int Nz = tile.mesh_lengths[2];

  // dB/dt = -curlE
  for(int ai=1; ai<=3; ai++) { //alphas
    Cx = 0.5*CXs(ai, 0, 0)*corr*tile.cfl;
    Cy = 0.5*CYs(0, ai, 0)*corr*tile.cfl;
    Cz = 0.5*CZs(0,  0,ai)*corr*tile.cfl;

    for(int k=0; k<Nz; k++) 
    for(int j=0; j<Ny; j++) 
    for(int i=0; i<Nx; i++){
        //legs with 0-offset
        mesh.bx(i,j,k) -= +Cy*Dp_y(mesh.ez,i,j,k, ai, 0,0,0); 
        mesh.bx(i,j,k) -= -Cz*Dp_z(mesh.ey,i,j,k, ai, 0,0,0);
                        
        mesh.by(i,j,k) -= -Cx*Dp_x(mesh.ez,i,j,k, ai, 0,0,0);
        mesh.by(i,j,k) -= +Cz*Dp_z(mesh.ex,i,j,k, ai, 0,0,0);
                        
        mesh.bz(i,j,k) -= +Cx*Dp_x(mesh.ey,i,j,k, ai, 0,0,0);
        mesh.bz(i,j,k) -= -Cy*Dp_y(mesh.ex,i,j,k, ai, 0,0,0);
    }

    //legs with n-offset
    //for(int bi : {-1,+1}) {
    for(int bi : {-3,-2,-1,+1,+2,+3}) {

        // hand-coded beta fetching from array
        Cx = 0.5*CXs(ai, 1, 1)*corr*tile.cfl;
        Cy = 0.5*CYs(1, ai, 1)*corr*tile.cfl;
        Cz = 0.5*CZs(1,  1,ai)*corr*tile.cfl;

        for(int k=0; k<Nz; k++) 
        for(int j=0; j<Ny; j++) 
        for(int i=0; i<Nx; i++){
            // curl_x
            mesh.bx(i,j,k) -= +Cy*Dp_y(mesh.ez,i,j,k, ai, bi,0,0); //beta_yx
            mesh.bx(i,j,k) -= +Cy*Dp_y(mesh.ez,i,j,k, ai, 0,0,bi); //beta_yz
                            
            mesh.bx(i,j,k) -= -Cz*Dp_z(mesh.ey,i,j,k, ai, bi,0,0);
            mesh.bx(i,j,k) -= -Cz*Dp_z(mesh.ey,i,j,k, ai, 0,bi,0);
                            
            // curl_y       
            mesh.by(i,j,k) -= -Cx*Dp_x(mesh.ez,i,j,k, ai, 0,bi,0);
            mesh.by(i,j,k) -= -Cx*Dp_x(mesh.ez,i,j,k, ai, 0,0,bi);
                            
            mesh.by(i,j,k) -= +Cz*Dp_z(mesh.ex,i,j,k, ai, bi,0,0);
            mesh.by(i,j,k) -= +Cz*Dp_z(mesh.ex,i,j,k, ai, 0,bi,0);
                            
            // curl_z       
            mesh.bz(i,j,k) -= +Cx*Dp_x(mesh.ey,i,j,k, ai, 0,bi,0);
            mesh.bz(i,j,k) -= +Cx*Dp_x(mesh.ey,i,j,k, ai, 0,0,bi);
                            
            mesh.bz(i,j,k) -= -Cy*Dp_y(mesh.ex,i,j,k, ai, bi,0,0);
            mesh.bz(i,j,k) -= -Cy*Dp_y(mesh.ex,i,j,k, ai, 0,0,bi);
        }
    }
  }

}


//template class fields::FDTDGen<1>;
//template class fields::FDTDGen<2>;
template class fields::FDTDGen<3>;
  
