#include "digital.h"

#include <cmath>


/// single 2D 2nd order binomial filter pass
template<>
void fields::Binomial2<2>::solve(fields::Tile<2>& tile)
{

  // 2D 1st order binomial coefficients
  double winv=1./16.,
         wtm=4.*winv, //middle
         wts=2.*winv, //side
         wtc=1.*winv; //corner

  YeeLattice& mesh = tile.get_yee();

  // using jx1, jy1, and jz1 as scratch arrays
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
    mesh.jx1(i,j,k) = 
      jx(i-1, j-1, k)*wtc + 
      jx(i  , j-1, k)*wts + 
      jx(i+1, j-1, k)*wtc + 
                          
      jx(i-1, j  , k)*wts + 
      jx(i  , j  , k)*wtm + 
      jx(i+1, j  , k)*wts + 

      jx(i-1, j+1, k)*wtc + 
      jx(i  , j+1, k)*wts + 
      jx(i+1, j+1, k)*wtc;
  }

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    mesh.jy1(i,j,k) = 
      jy(i-1, j-1, k)*wtc + 
      jy(i  , j-1, k)*wts + 
      jy(i+1, j-1, k)*wtc + 
                          
      jy(i-1, j  , k)*wts + 
      jy(i  , j  , k)*wtm + 
      jy(i+1, j  , k)*wts + 
        
      jy(i-1, j+1, k)*wtc + 
      jy(i  , j+1, k)*wts + 
      jy(i+1, j+1, k)*wtc;
  }

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    mesh.jz1(i,j,k) = 
      jz(i-1, j-1, k)*wtc + 
      jz(i  , j-1, k)*wts + 
      jz(i+1, j-1, k)*wtc + 
                          
      jz(i-1, j  , k)*wts + 
      jz(i  , j  , k)*wtm + 
      jz(i+1, j  , k)*wts + 
        
      jz(i-1, j+1, k)*wtc + 
      jz(i  , j+1, k)*wts + 
      jz(i+1, j+1, k)*wtc;
  }

  // then copy from scratch to original arrays

}









