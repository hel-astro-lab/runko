#include "fdtd_general.h"

#include <cmath>


/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */

/// 3D E pusher
template<>
void fields::FDTDGen<3>::push_e(fields::Tile<3>& tile)
{
  YeeLattice& mesh = tile.get_yee();

//	double C1 = coeff1*corr*tile.cfl;
//	double C2 = coeff2*corr*tile.cfl;
  
  real_short Cx, Cy, Cz;

  for(int is=0; is<=3; is++) 
  for(int js=0; js<=3; js++) 
  for(int ks=0; ks<=3; ks++) {
    if(is == 0 && js == 0 && ks == 0) continue;

    Cx = CXs(is, js, ks)*corr*tile.cfl;
    Cy = CYs(is, js, ks)*corr*tile.cfl;
    Cz = CZs(is, js, ks)*corr*tile.cfl;

    if(Cx == 0.0 && Cy == 0.0 && Cz == 0.0) continue;

    //Cx(1,0,0) = C1;
    //Cy(0,1,0) = C1;
    //Cz(0,0,1) = C1;
    //
    //Cx(2,0,0) = C2;
    //Cy(0,2,0) = C2;
    //Cz(0,0,2) = C2;

    // dE/dt = +curlB
    for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++)
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++)
    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

      // DONE inverse signs
        //-Cy*((mesh.bz(i,j-1,k)                  -mesh.bz(i,j,k)) )
         //+Cz*( mesh.by(i,j,k-1)                 -mesh.by(i,j,k)  )

      //mesh.ex(i,j,k)+=
      //    -Cy*(mesh.bz(i,j-2,k) - mesh.bz(i,j+1,k) )
      //    +Cz*(mesh.by(i,j,k-2) - mesh.by(i,j,k+1) )

      // NOTE: correcting for the derivate sign here so -= instead of +=
      mesh.ex(i,j,k) -= 
           +Cy*( mesh.bz(i-is, j-js, k-ks) - mesh.bz(i+is, j+js-1, k+ks) )
           -Cz*( mesh.by(i-is, j-js, k-ks) - mesh.by(i+is, j+js, k+ks-1) );
                                       
      mesh.ey(i,j,k) -=                
           -Cx*( mesh.bz(i-is, j-js, k-ks) - mesh.bz(i+is-1, j+js, k+ks) )
           +Cz*( mesh.bx(i-is, j-js, k-ks) - mesh.bx(i+is, j+js, k+ks-1) );
                                       
      mesh.ez(i,j,k) -=                
           -Cx*( mesh.by(i-is, j-js, k-ks) - mesh.by(i+is-1, j+js, k+ks) )
           +Cy*( mesh.bx(i-is, j-js, k-ks) - mesh.bx(i+is, j+js-1, k+ks) );
    } 
  }


    // alphas
    //Cx(1,0,0) = C1
    //Cy(0,1,0) = C1
    //Cz(0,0,1) = C1

    ////deltas
    //Cx(2,0,0) = C2
    //Cy(0,2,0) = C2
    //Cz(0,0,2) = C2

    //mesh.ex(i,j,k)+=
    //  C1*(mesh.by(i,j,k-1)-mesh.by(i,j,k)  - mesh.bz(i,j-1,k)+mesh.bz(i,j,k))+
    //  C2*(mesh.by(i,j,k-2)-mesh.by(i,j,k+1)- mesh.bz(i,j-2,k)+mesh.bz(i,j+1,k));

    //mesh.ey(i,j,k)+=
    //  C1*(mesh.bz(i-1,j,k)-mesh.bz(i,j,k)  - mesh.bx(i,j,k-1)+mesh.bx(i,j,k))+
    //  C2*(mesh.bz(i-2,j,k)-mesh.bz(i+1,j,k)- mesh.bx(i,j,k-2)+mesh.bx(i,j,k+1));

    //mesh.ez(i,j,k)+=
    //  C1*(mesh.bx(i,j-1,k)-mesh.bx(i,j,k)  - mesh.by(i-1,j,k)+mesh.by(i,j,k))+
    //  C2*(mesh.bx(i,j-2,k)-mesh.bx(i,j+1,k)- mesh.by(i-2,j,k)+mesh.by(i+1,j,k));
}


//--------------------------------------------------
/// Update B field with a half step

/// 3D B pusher
template<>
void fields::FDTDGen<3>::push_half_b(fields::Tile<3>& tile)
{
  YeeLattice& mesh = tile.get_yee();

  real_short Cx, Cy, Cz;

  //std::cout << "starting reading\n";
  //for(int is=1; is<=3; is++) 
  //for(int js=1; js<=3; js++) 
  //for(int ks=1; ks<=3; ks++) {

  //  Cx = CXs(is-1, js-1, ks-1);
  //  Cy = CYs(is-1, js-1, ks-1);
  //  Cz = CZs(is-1, js-1, ks-1);

  //  std::cout << "CXs " 
  //      << is << ","
  //      << js << ","
  //      << ks << " :"
  //      << Cx << std::endl;
  //      
  //  std::cout << "CYs " 
  //      << is << ","
  //      << js << ","
  //      << ks << " :"
  //      << Cy << std::endl;

  //  std::cout << "CZs " 
  //      << is << ","
  //      << js << ","
  //      << ks << " :"
  //      << Cz << std::endl;
  //}

    // dB/dt = -curlE
    for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++)
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++)
    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

  for(int is=0; is<=3; is++) 
  for(int js=0; js<=3; js++) 
  for(int ks=0; ks<=3; ks++) {
    //if(is == 0 && js == 0 && ks == 0) continue;

    Cx = 0.5*CXs(is, js, ks)*corr*tile.cfl;
    Cy = 0.5*CYs(is, js, ks)*corr*tile.cfl;
    Cz = 0.5*CZs(is, js, ks)*corr*tile.cfl;

    //if(Cx == 0.0 && Cy == 0.0 && Cz == 0.0) continue;


            // DONE
            //-=
            //+Cy*(mesh.ez(i,j+1,k)-mesh.ez(i,j,k) )
            //-Cz*(mesh.ey(i,j,k+1)-mesh.ey(i,j,k) )

      //real_short DT= +Cy*( mesh.ez(i+is, j+js, k+ks) - mesh.ez(i-is+1, j-js+1, k-ks+1) ) -Cz*( mesh.ey(i+is, j+js, k+ks) - mesh.ey(i-is+1, j-js+1, k-ks+1) );

      //if (DT*DT > 1.0e-5) {
      //td::cout << 
      //      " meshbx:" << mesh.bx(i,j,k) 
      //   << " CyT:" << +Cy*( mesh.ez(i+is, j+js, k+ks) - mesh.ez(i-is+1, j-js+1, k-ks+1) ) 
      //   << " CzT:" << -Cz*( mesh.ey(i+is, j+js, k+ks) - mesh.ey(i-is+1, j-js+1, k-ks+1) )
      //   << " ey:"  << mesh.ey(i+is, j+js, k+ks)
      //   << " ez:"  << mesh.ez(i+is, j+js, k+ks)
      //   << " ez1:" << mesh.ez(i-is+1, j-js+1, k-ks+1) 
      //   << " ey1:" << -mesh.ey(i-is+1, j-js+1, k-ks+1) 
      //   << " is:"  << is
      //   << " js:"  << js
      //   << " ks:"  << ks
      //   << " i:"  << i
      //   << " j:"  << j
      //   << " k:"  << k
      //   << std::endl;
      // }
        
      //mesh.bx(i,j,k)+=
      //     +C1*(
      //        - ( mesh.ez(i,j+1,k)-mesh.ez(i,j,k) )
      //          ( mesh.ey(i,j,k+1)-mesh.ey(i,j,k) )
      //        )
      //     +C2*(
      //          mesh.ey(i,j,k+2)-mesh.ey(i,j,k-1)
      //        - mesh.ez(i,j+2,k)+mesh.ez(i,j-1,k)
      //        );
        
      // DONE
      // mesh.bx(i,j,k)+=
      // C1*(
      //     +(mesh.ez(i,j,k)) - mesh.ez(i,j+1,k)
      //     -(mesh.ey(i,j,k) - mesh.ey(i,j,k+1))

      // FDTD2 DONE
      //mesh.bx(i,j,k) += 
      //    + C*( mesh.ez(i,j,k) - mesh.ez(i,  j+1,k  )  )
      //    - C*( mesh.ey(i,j,k) - mesh.ey(i,  j,  k+1)  )

      // NOTE: taking the derivative sign into account here so += instead -=
      mesh.bx(i,j,k) += 
           +Cy*( mesh.ez(i-is, j-js+1, k-ks) - mesh.ez(i+is, j+js, k+ks) )
           -Cz*( mesh.ey(i-is, j-js, k-ks+1) - mesh.ey(i+is, j+js, k+ks) );

      mesh.by(i,j,k) += 
           -Cx*( mesh.ez(i-is+1, j-js, k-ks) - mesh.ez(i+is, j+js, k+ks) )
           +Cz*( mesh.ex(i-is, j-js, k-ks+1) - mesh.ex(i+is, j+js, k+ks) );

      mesh.bz(i,j,k) += 
           +Cx*( mesh.ey(i-is+1, j-js, k-ks) - mesh.ey(i+is, j+js, k+ks) )
           -Cy*( mesh.ex(i-is, j-js+1, k-ks) - mesh.ex(i+is, j+js, k+ks) );

    } 
  }

  //double C1 = 0.5*coeff1*corr*tile.cfl;
  //double C2 = 0.5*coeff2*corr*tile.cfl;

  //for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) 
  //for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  //for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

  //  	mesh.bx(i,j,k)+=
  //    C1*(mesh.ey(i,j,k+1)-mesh.ey(i,j,k)  - mesh.ez(i,j+1,k)+mesh.ez(i,j,k))+
  //    C2*(mesh.ey(i,j,k+2)-mesh.ey(i,j,k-1)- mesh.ez(i,j+2,k)+mesh.ez(i,j-1,k));


  //  	mesh.bz(i,j,k)+=
  //    C1*(mesh.ex(i,j+1,k)-mesh.ex(i,j,k)  - mesh.ey(i+1,j,k)+mesh.ey(i,j,k))+
  //    C2*(mesh.ex(i,j+2,k)-mesh.ex(i,j-1,k)- mesh.ey(i+2,j,k)+mesh.ey(i-1,j,k));

  //}


}



//template class fields::FDTDGen<1>;
//template class fields::FDTDGen<2>;
template class fields::FDTDGen<3>;
  
