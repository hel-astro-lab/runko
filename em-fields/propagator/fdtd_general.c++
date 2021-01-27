#include "fdtd_general.h"

#include <cmath>


//push_e
//C1*( mesh.by(i,j,k-1)-mesh.by(i,j,k) 
//   -(mesh.bz(i,j-1,k)-mesh.bz(i,j,k)) 
//
//C2*( mesh.by(i,j,k-2)-mesh.by(i,j,k+1)
//   -(mesh.bz(i,j-2,k)+mesh.bz(i,j+1,k))

//     -f(i-1 , j   , k     + f(i     , j   , k   );
//     -f(i-2 , j   , k     + f(i+1   , j   , k   );
inline real_short Dm_x( toolbox::Mesh<real_short, 3>& f, int i, int j, int k, int ai, int /*bi*/, int bj, int bk) {
    return f(i+ai-1, j+bj, k+bk) -f(i-ai, j+bj, k+bk);
}

inline real_short Dm_y( toolbox::Mesh<real_short, 3>& f, int i, int j, int k, int ai, int bi, int /*bj*/, int bk) {
    return f(i+bi, j+ai-1, k+bk) -f(i+bi, j-ai, k+bk);
}

inline real_short Dm_z( toolbox::Mesh<real_short, 3>& f, int i, int j, int k, int ai, int bi, int bj, int /*bk*/) {
    return f(i+bi, j+bj, k+ai-1) -f(i+bi, j+bj, k-ai);
}

//push half b
//C1*( mesh.ey(i,j,k+1)-mesh.ey(i,j,k)  
//  - (mesh.ez(i,j+1,k)-mesh.ez(i,j,k))
//
//C2*( mesh.ey(i,j,k+2)-mesh.ey(i,j,k-1)
//   - mesh.ez(i,j+2,k)+mesh.ez(i,j-1,k));

//     -f(i+1 , j   , k     + f(i     , j   , k   );
//     -f(i+2 , j   , k     + f(i-1   , j   , k   );
inline real_short Dp_x( toolbox::Mesh<real_short, 3>& f, int i, int j, int k, int ai, int /*bi*/, int bj, int bk) {
    return f(i+ai, j+bj, k+bk) - f(i-ai+1, j+bj, k+bk);
}

inline real_short Dp_y( toolbox::Mesh<real_short, 3>& f, int i, int j, int k, int ai, int bi, int /*bj*/, int bk) {
    return f(i+bi, j+ai, k+bk) - f(i+bi, j-ai+1, k+bk);
}

inline real_short Dp_z( toolbox::Mesh<real_short, 3>& f, int i, int j, int k, int ai, int bi, int bj, int /*bk*/) {
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
  
  real_short Cx, Cy, Cz;

  const int Nx = tile.mesh_lengths[0];
  const int Ny = tile.mesh_lengths[1];
  const int Nz = tile.mesh_lengths[2];


  //for(int ks=0; ks<=1; ks++) {
  //for(int js=0; js<=1; js++) {
  //for(int is=0; is<=1; is++) {
  //  //if(is == 0 && js == 0 && ks == 0) continue;

  //  Cx = CXs(is, js, ks)*corr*tile.cfl;
  //  Cy = CYs(is, js, ks)*corr*tile.cfl;
  //  Cz = CZs(is, js, ks)*corr*tile.cfl;

  //  //if(Cx == 0.0 && Cy == 0.0 && Cz == 0.0) continue;

  //  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
  //  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
  //  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

  //    // dE/dt = +curlB
  //    // NOTE: correcting for the derivate sign here so -= instead of +=
  //    mesh.ex(i,j,k) -= +Cy*( mesh.bz(i-is, j-js, k-ks) - mesh.bz(i+is, j+js-1, k+ks) );
  //    mesh.ex(i,j,k) -= -Cz*( mesh.by(i-is, j-js, k-ks) - mesh.by(i+is, j+js, k+ks-1) );
  //                                     
  //    mesh.ey(i,j,k) -= -Cx*( mesh.bz(i-is, j-js, k-ks) - mesh.bz(i+is-1, j+js, k+ks) );
  //    mesh.ey(i,j,k) -= +Cz*( mesh.bx(i-is, j-js, k-ks) - mesh.bx(i+is, j+js, k+ks-1) );
  //                                     
  //    mesh.ez(i,j,k) -= -Cx*( mesh.by(i-is, j-js, k-ks) - mesh.by(i+is-1, j+js, k+ks) );
  //    mesh.ez(i,j,k) -= +Cy*( mesh.bx(i-is, j-js, k-ks) - mesh.bx(i+is, j+js-1, k+ks) );
  //  }}}
  //}}}


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
    //for(int bi : {-3,-2,-1,+1,+2,+3}) {
    for(int bi : {-1,+1}) {

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

  real_short Cx, Cy, Cz;

  const int Nx = tile.mesh_lengths[0];
  const int Ny = tile.mesh_lengths[1];
  const int Nz = tile.mesh_lengths[2];


  //if(i == 0 && j == 0 && k == 0) {

  //        real_short de1 = mesh.ez(i-is, j-js+1, k-ks) - mesh.ez(i+is, j+js, k+ks);
  //        real_short de2 = mesh.ey(i-is, j-js, k-ks+1) - mesh.ey(i+is, j+js, k+ks); 
  //        real_short de3 = mesh.ez(i-is+1, j-js, k-ks) - mesh.ez(i+is, j+js, k+ks); 
  //        real_short de4 = mesh.ex(i-is, j-js, k-ks+1) - mesh.ex(i+is, j+js, k+ks); 
  //        real_short de5 = mesh.ey(i-is+1, j-js, k-ks) - mesh.ey(i+is, j+js, k+ks); 
  //        real_short de6 = mesh.ex(i-is, j-js+1, k-ks) - mesh.ex(i+is, j+js, k+ks); 

  //    std::cout << " Bupd>> :"
  //        << " ijks: (" << is << "," << js << "," << ks << ")"
  //        << " Cs: [" << Cx << "," << Cy << "," << Cz << "] \n" 

  //        << "bx: +Cy ez: (" << i-is   << "," << j-js+1  << "," << k-ks  << ")" << "Delta:" << de1 << "\n"
  //        << "bx: -Cz ey: (" << i-is   << "," << j-js    << "," << k-ks+1<< ")" << "Delta:" << de2 << "\n"

  //        << "by: -Cx ez: (" << i-is+1 << "," << j-js    << "," << k-ks  << ")" << "Delta:" << de3 << "\n"
  //        << "by: +Cz ex: (" << i-is   << "," << j-js    << "," << k-ks+1<< ")" << "Delta:" << de4 << "\n"

  //        << "bz: +Cx ey: (" << i-is+1 << "," << j-js    << "," << k-ks  << ")" << "Delta:" << de5 << "\n"
  //        << "bz: -Cy ex: (" << i-is   << "," << j-js+1  << "," << k-ks  << ")" << "Delta:" << de6 << "\n";

  //}

  // dB/dt = -curlE
  // NOTE: taking the derivative sign into account here so += instead of -=

  //Realf C = 0.5 * tile.cfl * dt * corr;

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


        //real_short dbx = 0.0;
        //dbx -= +Cy*Dp_y(mesh.ez,i,j,k, ai, 0,0,0); 
        //dbx -= -Cz*Dp_z(mesh.ey,i,j,k, ai, 0,0,0);

        //real_short dby = 0.0;
        //dby -= -Cx*Dp_x(mesh.ez,i,j,k, ai, 0,0,0);
        //dby -= +Cz*Dp_z(mesh.ex,i,j,k, ai, 0,0,0);

        //real_short dbz = 0.0;
        //dbz -= +Cx*Dp_x(mesh.ey,i,j,k, ai, 0,0,0);
        //dbz -= -Cy*Dp_y(mesh.ex,i,j,k, ai, 0,0,0);

        ////--------------------------------------------------
        //real_short dbx1 = 0.0;
        //real_short dby1 = 0.0;
        //real_short dbz1 = 0.0;
        //dbx1 += + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k)) + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));
        //dby1 += + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k)) + C*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k));
        //dbz1 += + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)) + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));


        //if(i == 0 && j == 0 && k == 0) {
        //    std::cout <<
        //        " dbx " << dbx << " vs " << dbx1 << "\n"
        //        " dby " << dby << " vs " << dby1 << "\n"
        //        " dbz " << dbz << " vs " << dbz1 << "\n";
        //}

  //if(i == 0 && j == 0 && k == 0) {
  //        real_short de1 = +Cy*Dp_y(mesh.ez,i,j,k, ai, 0,0,0); 
  //        real_short de2 = -Cz*Dp_z(mesh.ey,i,j,k, ai, 0,0,0);  
  //                                                             
  //        real_short de3 = -Cx*Dp_x(mesh.ez,i,j,k, ai, 0,0,0);  
  //        real_short de4 = +Cz*Dp_z(mesh.ex,i,j,k, ai, 0,0,0);  
  //                                                             
  //        real_short de5 = +Cx*Dp_x(mesh.ey,i,j,k, ai, 0,0,0);  
  //        real_short de6 = -Cy*Dp_y(mesh.ex,i,j,k, ai, 0,0,0);  

  //        int bj = 0, bk = 0;
  //    std::cout << " Bupd>> :"
  //        //<< " ijks: (" << is << "," << js << "," << ks << ")"
  //        << " ijks: (" << ai << "," << ai << "," << ai << ")"
  //        << " Cs: [" << Cx << "," << Cy << "," << Cz << "] \n" 

  //        << " ijk: (" << i+ai <<","<< j+bj   <<","<< k+bk   <<") ("<<  i-ai+1<<","<< j+bj     <<","<< k+bk     <<"\n"
  //        << " ijk: (" << i    <<","<< j+bj+ai<<","<< k+bk   <<") ("<<  i     <<","<< j+bj-ai+1<<","<< k+bk     <<"\n"
  //        << " ijk: (" << i    <<","<< j+bj   <<","<< k+bk+ai<<") ("<<  i     <<","<< j+bj     <<","<< k+bk-ai+1<<"\n"

  //        << "bx: +Cy ez:  Delta:" << de1 << "\n"
  //        << "bx: -Cz ey:  Delta:" << de2 << "\n"

  //        << "by: -Cx ez:  Delta:" << de3 << "\n"
  //        << "by: +Cz ex:  Delta:" << de4 << "\n"

  //        << "bz: +Cx ey:  Delta:" << de5 << "\n"
  //        << "bz: -Cy ex:  Delta:" << de6 << "\n";

  //}

    }

    //legs with n-offset
    //for(int bi : {-3,-2,-1,+1,+2,+3}) {
    for(int bi : {-1,+1}) {

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
  
