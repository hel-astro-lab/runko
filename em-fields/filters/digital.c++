#include "digital.h"

#include <cmath>


/// single 2D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D 1st order binomial coefficients
  double winv=1./16.,
         wtm=4.*winv, //middle
         wts=2.*winv, //side
         wtc=1.*winv; //corner

  auto& mesh = tile.get_yee();

  // using tmp as scratch arrays
  //
  // values are processed in the following order
  //
  // | 7 8 9 |
  // | 4 5 6 |
  // | 1 2 3 |
  // 
  // because this access pattern translates to continuous sweeps of
  // ...,1,2,3,....,4,5,6,....,7,8,9,.....
  // for the target array memory.


  //--------------------------------------------------
  // Jx
  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jx(i-1, j-1, k)*wtc + 
      mesh.jx(i  , j-1, k)*wts + 
      mesh.jx(i+1, j-1, k)*wtc + 

      mesh.jx(i-1, j  , k)*wts + 
      mesh.jx(i  , j  , k)*wtm + 
      mesh.jx(i+1, j  , k)*wts + 

      mesh.jx(i-1, j+1, k)*wtc + 
      mesh.jx(i  , j+1, k)*wts + 
      mesh.jx(i+1, j+1, k)*wtc;
  }
  mesh.jx = tmp; // then copy from scratch to original arrays

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jy(i-1, j-1, k)*wtc + 
      mesh.jy(i  , j-1, k)*wts + 
      mesh.jy(i+1, j-1, k)*wtc + 

      mesh.jy(i-1, j  , k)*wts + 
      mesh.jy(i  , j  , k)*wtm + 
      mesh.jy(i+1, j  , k)*wts + 

      mesh.jy(i-1, j+1, k)*wtc + 
      mesh.jy(i  , j+1, k)*wts + 
      mesh.jy(i+1, j+1, k)*wtc;
  }
  mesh.jy = tmp; // then copy from scratch to original arrays

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jz(i-1, j-1, k)*wtc + 
      mesh.jz(i  , j-1, k)*wts + 
      mesh.jz(i+1, j-1, k)*wtc + 

      mesh.jz(i-1, j  , k)*wts + 
      mesh.jz(i  , j  , k)*wtm + 
      mesh.jz(i+1, j  , k)*wts + 

      mesh.jz(i-1, j+1, k)*wtc + 
      mesh.jz(i  , j+1, k)*wts + 
      mesh.jz(i+1, j+1, k)*wtc;
  }
  mesh.jz = tmp; // then copy from scratch to original arrays

}



/// single 2D 3-point general filter pass
template<>
void fields::General3p<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  double winv=1./4.;                         //normalization
  double wtm=winv * 4.0*alpha*alpha,         //middle
         wts=winv * 2.0*alpha*(1.0-alpha),   //side
         wtc=winv * (1.0-alpha)*(1.0-alpha); //corner

  auto& mesh = tile.get_yee();


  //--------------------------------------------------
  // Jx
  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jx(i-1, j-1, k)*wtc + 
      mesh.jx(i  , j-1, k)*wts + 
      mesh.jx(i+1, j-1, k)*wtc + 

      mesh.jx(i-1, j  , k)*wts + 
      mesh.jx(i  , j  , k)*wtm + 
      mesh.jx(i+1, j  , k)*wts + 

      mesh.jx(i-1, j+1, k)*wtc + 
      mesh.jx(i  , j+1, k)*wts + 
      mesh.jx(i+1, j+1, k)*wtc;
  }
  mesh.jx = tmp; // then copy from scratch to original arrays

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jy(i-1, j-1, k)*wtc + 
      mesh.jy(i  , j-1, k)*wts + 
      mesh.jy(i+1, j-1, k)*wtc + 

      mesh.jy(i-1, j  , k)*wts + 
      mesh.jy(i  , j  , k)*wtm + 
      mesh.jy(i+1, j  , k)*wts + 

      mesh.jy(i-1, j+1, k)*wtc + 
      mesh.jy(i  , j+1, k)*wts + 
      mesh.jy(i+1, j+1, k)*wtc;
  }
  mesh.jy = tmp; // then copy from scratch to original arrays

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jz(i-1, j-1, k)*wtc + 
      mesh.jz(i  , j-1, k)*wts + 
      mesh.jz(i+1, j-1, k)*wtc + 

      mesh.jz(i-1, j  , k)*wts + 
      mesh.jz(i  , j  , k)*wtm + 
      mesh.jz(i+1, j  , k)*wts + 

      mesh.jz(i-1, j+1, k)*wtc + 
      mesh.jz(i  , j+1, k)*wts + 
      mesh.jz(i+1, j+1, k)*wtc;
  }
  mesh.jz = tmp; // then copy from scratch to original arrays

}


/// single 2D 3-point general filter pass
template<>
void fields::General3pStrided<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  double winv=1./4.;                         //normalization
  double wtm=winv * 4.0*alpha*alpha,         //middle
         wts=winv * 2.0*alpha*(1.0-alpha),   //side
         wtc=winv * (1.0-alpha)*(1.0-alpha); //corner

  auto& mesh = tile.get_yee();


  //--------------------------------------------------
  // Jx
  int k = 0;
  for(int istr=stride; istr>0; istr--) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

      tmp(i,j,k) = 
        mesh.jx(i-istr, j-istr, k)*wtc + 
        mesh.jx(i     , j-istr, k)*wts + 
        mesh.jx(i+istr, j-istr, k)*wtc + 

        mesh.jx(i-istr, j     , k)*wts + 
        mesh.jx(i     , j     , k)*wtm + 
        mesh.jx(i+istr, j     , k)*wts + 

        mesh.jx(i-istr, j+istr, k)*wtc + 
        mesh.jx(i     , j+istr, k)*wts + 
        mesh.jx(i+istr, j+istr, k)*wtc;
    }
    mesh.jx = tmp; // then copy from scratch to original arrays
  }

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    //for(int istr=1; istr<=stride; istr++) {
    for(int istr=stride; istr>0; istr--) {

      tmp(i,j,k) = 
        mesh.jy(i-istr, j-istr, k)*wtc + 
        mesh.jy(i     , j-istr, k)*wts + 
        mesh.jy(i+istr, j-istr, k)*wtc + 
               
        mesh.jy(i-istr, j     , k)*wts + 
        mesh.jy(i     , j     , k)*wtm + 
        mesh.jy(i+istr, j     , k)*wts + 
               
        mesh.jy(i-istr, j+istr, k)*wtc + 
        mesh.jy(i     , j+istr, k)*wts + 
        mesh.jy(i+istr, j+istr, k)*wtc;
    }


    //for(int jstr=1; jstr<=stride; jstr++) 
    //for(int istr=1; istr<=stride; istr++) {
    //  tmp(i,j,k) = 
    //    mesh.jy(i-istr, j-jstr, k)*wtc + 
    //    mesh.jy(i     , j-jstr, k)*wts + 
    //    mesh.jy(i+istr, j-jstr, k)*wtc + 
    //           
    //    mesh.jy(i-istr, j     , k)*wts + 
    //    mesh.jy(i     , j     , k)*wtm + 
    //    mesh.jy(i+istr, j     , k)*wts + 
    //           
    //    mesh.jy(i-istr, j+jstr, k)*wtc + 
    //    mesh.jy(i     , j+jstr, k)*wts + 
    //    mesh.jy(i+istr, j+jstr, k)*wtc;
    //}
  }
  mesh.jy = tmp; // then copy from scratch to original arrays

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {

    //for(int istr=1; istr<=stride; istr++) {
    for(int istr=stride; istr>0; istr--) {

      tmp(i,j,k) = 
        mesh.jz(i-istr, j-istr, k)*wtc + 
        mesh.jz(i     , j-istr, k)*wts + 
        mesh.jz(i+istr, j-istr, k)*wtc + 
               
        mesh.jz(i-istr, j     , k)*wts + 
        mesh.jz(i     , j     , k)*wtm + 
        mesh.jz(i+istr, j     , k)*wts + 
               
        mesh.jz(i-istr, j+istr, k)*wtc + 
        mesh.jz(i     , j+istr, k)*wts + 
        mesh.jz(i+istr, j+istr, k)*wtc;
    }

    //for(int jstr=1; jstr<=stride; jstr++) 
    //for(int istr=1; istr<=stride; istr++) {
    //  tmp(i,j,k) = 
    //    mesh.jz(i-istr, j-jstr, k)*wtc + 
    //    mesh.jz(i     , j-jstr, k)*wts + 
    //    mesh.jz(i+istr, j-jstr, k)*wtc + 
    //           
    //    mesh.jz(i-istr, j     , k)*wts + 
    //    mesh.jz(i     , j     , k)*wtm + 
    //    mesh.jz(i+istr, j     , k)*wts + 
    //           
    //    mesh.jz(i-istr, j+jstr, k)*wtc + 
    //    mesh.jz(i     , j+jstr, k)*wts + 
    //    mesh.jz(i+istr, j+jstr, k)*wtc;
    //}
  }
  mesh.jz = tmp; // then copy from scratch to original arrays

}





//template class fields::Binomial2<1>; // 1D
template class fields::Binomial2<2>; // 2D
//template class fields::Binomial2<3>; // 3D


//template class fields::General3p<1>; // 1D
template class fields::General3p<2>; // 2D
//template class fields::General3p<3>; // 3D


//template class fields::General3pStrided<1>; // 1D
template class fields::General3pStrided<2>; // 2D
//template class fields::General3pStrided<3>; // 3D
