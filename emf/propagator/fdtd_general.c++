#include "fdtd_general.h"

#include <cmath>
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


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
void emf::FDTDGen<3>::push_e(emf::Tile<3>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, Grids &mesh)
  {

    // dE/dt = +curlB
    for(int ai=1; ai<=3; ai++) { //alphas
      const float_m Cx = CXs(ai, 0, 0)*corr*tile.cfl;
      const float_m Cy = CYs(0, ai, 0)*corr*tile.cfl;
      const float_m Cz = CZs(0,  0,ai)*corr*tile.cfl;

      //legs with 0-offset
      mesh.ex(i,j,k) += +Cy*Dm_y(mesh.bz,i,j,k, ai, 0,0,0); 
      mesh.ex(i,j,k) += -Cz*Dm_z(mesh.by,i,j,k, ai, 0,0,0);
      mesh.ey(i,j,k) += -Cx*Dm_x(mesh.bz,i,j,k, ai, 0,0,0);
      mesh.ey(i,j,k) += +Cz*Dm_z(mesh.bx,i,j,k, ai, 0,0,0);
      mesh.ez(i,j,k) += +Cx*Dm_x(mesh.by,i,j,k, ai, 0,0,0);
      mesh.ez(i,j,k) += -Cy*Dm_y(mesh.bx,i,j,k, ai, 0,0,0);

      //legs with n-offset
      //for(int bi : {-1,+1}) {
      for(int bi : {-3,-2,-1,+1,+2,+3}) {

          // hand-coded beta fetching from array
          const float_m Cx2 = CXs(ai, 1, 1)*corr*tile.cfl;
          const float_m Cy2 = CYs(1, ai, 1)*corr*tile.cfl;
          const float_m Cz2 = CZs(1,  1,ai)*corr*tile.cfl;

          // curl_x
          mesh.ex(i,j,k) += +Cy2*Dm_y(mesh.bz,i,j,k, ai, bi,0,0); //beta_yx
          mesh.ex(i,j,k) += +Cy2*Dm_y(mesh.bz,i,j,k, ai, 0,0,bi); //beta_yz

          mesh.ex(i,j,k) += -Cz2*Dm_z(mesh.by,i,j,k, ai, bi,0,0);
          mesh.ex(i,j,k) += -Cz2*Dm_z(mesh.by,i,j,k, ai, 0,bi,0);

          // curl_y
          mesh.ey(i,j,k) += -Cx2*Dm_x(mesh.bz,i,j,k, ai, 0,bi,0);
          mesh.ey(i,j,k) += -Cx2*Dm_x(mesh.bz,i,j,k, ai, 0,0,bi);

          mesh.ey(i,j,k) += +Cz2*Dm_z(mesh.bx,i,j,k, ai, bi,0,0);
          mesh.ey(i,j,k) += +Cz2*Dm_z(mesh.bx,i,j,k, ai, 0,bi,0);

          // curl_z
          mesh.ez(i,j,k) += +Cx2*Dm_x(mesh.by,i,j,k, ai, 0,bi,0);
          mesh.ez(i,j,k) += +Cx2*Dm_x(mesh.by,i,j,k, ai, 0,0,bi);

          mesh.ez(i,j,k) += -Cy2*Dm_y(mesh.bx,i,j,k, ai, bi,0,0);
          mesh.ez(i,j,k) += -Cy2*Dm_y(mesh.bx,i,j,k, ai, 0,0,bi);

      }// end of bi
    } //end of ai

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

/// 3D B pusher
template<>
void emf::FDTDGen<3>::push_half_b(emf::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  Grids& mesh = tile.get_grids();

  UniIter::iterate3D(
  [=] DEVCALLABLE (int i, int j, int k, Grids &mesh)
  {
    // dB/dt = -curlE
    for(int ai=1; ai<=3; ai++) { //alphas
      const float_m Cx = 0.5*CXs(ai, 0, 0)*corr*tile.cfl;
      const float_m Cy = 0.5*CYs(0, ai, 0)*corr*tile.cfl;
      const float_m Cz = 0.5*CZs(0,  0,ai)*corr*tile.cfl;

      //legs with 0-offset
      mesh.bx(i,j,k) -= +Cy*Dp_y(mesh.ez,i,j,k, ai, 0,0,0); 
      mesh.bx(i,j,k) -= -Cz*Dp_z(mesh.ey,i,j,k, ai, 0,0,0);
                      
      mesh.by(i,j,k) -= -Cx*Dp_x(mesh.ez,i,j,k, ai, 0,0,0);
      mesh.by(i,j,k) -= +Cz*Dp_z(mesh.ex,i,j,k, ai, 0,0,0);
                      
      mesh.bz(i,j,k) -= +Cx*Dp_x(mesh.ey,i,j,k, ai, 0,0,0);
      mesh.bz(i,j,k) -= -Cy*Dp_y(mesh.ex,i,j,k, ai, 0,0,0);


      //legs with n-offset
      //for(int bi : {-1,+1}) {
      for(int bi : {-3,-2,-1,+1,+2,+3}) {

        // hand-coded beta fetching from array
        const float_m Cx2 = 0.5*CXs(ai, 1, 1)*corr*tile.cfl;
        const float_m Cy2 = 0.5*CYs(1, ai, 1)*corr*tile.cfl;
        const float_m Cz2 = 0.5*CZs(1,  1,ai)*corr*tile.cfl;

        // curl_x
        mesh.bx(i,j,k) -= +Cy2*Dp_y(mesh.ez,i,j,k, ai, bi,0,0); //beta_yx
        mesh.bx(i,j,k) -= +Cy2*Dp_y(mesh.ez,i,j,k, ai, 0,0,bi); //beta_yz
                        
        mesh.bx(i,j,k) -= -Cz2*Dp_z(mesh.ey,i,j,k, ai, bi,0,0);
        mesh.bx(i,j,k) -= -Cz2*Dp_z(mesh.ey,i,j,k, ai, 0,bi,0);
                        
        // curl_y       
        mesh.by(i,j,k) -= -Cx2*Dp_x(mesh.ez,i,j,k, ai, 0,bi,0);
        mesh.by(i,j,k) -= -Cx2*Dp_x(mesh.ez,i,j,k, ai, 0,0,bi);
                        
        mesh.by(i,j,k) -= +Cz2*Dp_z(mesh.ex,i,j,k, ai, bi,0,0);
        mesh.by(i,j,k) -= +Cz2*Dp_z(mesh.ex,i,j,k, ai, 0,bi,0);
                        
        // curl_z       
        mesh.bz(i,j,k) -= +Cx2*Dp_x(mesh.ey,i,j,k, ai, 0,bi,0);
        mesh.bz(i,j,k) -= +Cx2*Dp_x(mesh.ey,i,j,k, ai, 0,0,bi);
                        
        mesh.bz(i,j,k) -= -Cy2*Dp_y(mesh.ex,i,j,k, ai, bi,0,0);
        mesh.bz(i,j,k) -= -Cy2*Dp_y(mesh.ex,i,j,k, ai, 0,0,bi);
      } //end of bi
    } //end of ai
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


//template class emf::FDTDGen<1>;
//template class emf::FDTDGen<2>;
template class emf::FDTDGen<3>;
  
