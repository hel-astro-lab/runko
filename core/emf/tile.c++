#include <iostream>
#include <cmath>

#include "core/emf/tile.h"

#include "tools/has_element.h"
#include "external/iter/iter.h"
#include "external/iter/allocator.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


// Optional vectorization of tile boundary updates;
// VEC_FLD controls exchange_boundary routines 
// VEC_CUR controls exchange_currents routines
// Latter has atomic_add operations that can be very slow on some platforms
//#define VEC_FLD2D 
//#define VEC_FLD3D
//#define VEC_CUR2D 
//#define VEC_CUR3D 



namespace emf {
  using namespace mpi4cpp;

// SPHINX emf docs addcur start
// Deposit current into electric field
template<std::size_t D>
void Tile<D>::deposit_current() 
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& mesh = get_grids();

  UniIter::iterate3D(
    [=] DEVCALLABLE (int i, int j, int k, Grids& mesh)
    {
      mesh.ex(i,j,k) -= mesh.jx(i,j,k);
      mesh.ey(i,j,k) -= mesh.jy(i,j,k);
      mesh.ez(i,j,k) -= mesh.jz(i,j,k);
    },mesh.ex.Nx, mesh.ex.Ny, mesh.ex.Nz, mesh);

  UniIter::sync();


#ifdef GPU
  nvtxRangePop();
#endif 
}
// SPHINX emf docs addcur stop



/// Update Yee grid boundaries
template<>
void Tile<1>::update_boundaries(
        corgi::Grid<1>& grid,
        std::vector<int> iarr
        ) 
{
  using Tile_t  = Tile<1>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito=0, ifro=0;
  const int halo = 3; // halo region size for fields
  Tileptr tpr; 

  // target
  auto& lhs = get_grids();

  const int Nx = lhs.Nx;


  for(int in=-1; in <= 1; in++) {
    if (in == 0) continue;

    tpr = std::dynamic_pointer_cast<Tile<1>>(grid.get_tileptr( neighs(in) ));

    if (tpr) {
      auto& rhs = tpr->get_grids();

      // copy from right side to left
      //copy_vert_grids(mesh, mleft, iarr, -1, mleft.Nx-1); 
        
      // copy from left side to right
      //copy_vert_grids(mesh, mright, iarr, mesh.Nx, 0); 

      if (in == +1) { ito = Nx; ifro = 0; }
      if (in == -1) { ito = -1; ifro = Nx-1; }

      if(has_elem(iarr, 0)) {
        UniIter::iterate([=] DEVCALLABLE (int h, Grids &lhs_in, Grids &rhs_in){
          lhs_in.jx(ito+in*h, 0, 0) = rhs_in.jx(ifro+in*h, 0, 0);
          lhs_in.jy(ito+in*h, 0, 0) = rhs_in.jy(ifro+in*h, 0, 0);
          lhs_in.jz(ito+in*h, 0, 0) = rhs_in.jz(ifro+in*h, 0, 0);
        }, halo, lhs, rhs);
      }

      if(has_elem(iarr, 1)) {
        UniIter::iterate([=] DEVCALLABLE (int h, Grids &lhs_in, Grids &rhs_in){
          lhs_in.ex(ito+in*h, 0, 0) = rhs_in.ex(ifro+in*h, 0, 0);
          lhs_in.ey(ito+in*h, 0, 0) = rhs_in.ey(ifro+in*h, 0, 0);
          lhs_in.ez(ito+in*h, 0, 0) = rhs_in.ez(ifro+in*h, 0, 0);
        }, halo, lhs, rhs);
      }

      if(has_elem(iarr, 2)) {
        UniIter::iterate([=] DEVCALLABLE (int h, Grids &lhs_in, Grids &rhs_in){
          lhs_in.bx(ito+in*h, 0, 0) = rhs_in.bx(ifro+in*h, 0, 0);
          lhs_in.by(ito+in*h, 0, 0) = rhs_in.by(ifro+in*h, 0, 0);
          lhs_in.bz(ito+in*h, 0, 0) = rhs_in.bz(ifro+in*h, 0, 0);
        }, halo, lhs, rhs);
      }

    }
  }

}


template<>
void Tile<2>::update_boundaries(
        corgi::Grid<2>& grid,
        std::vector<int> iarr
        ) 
{
  using Tile_t  = Tile<2>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito=0, jto=0, ifro=0, jfro=0;
  Tileptr tpr;

  auto& lhs = get_grids(); // target as a reference to update into

  const int Nx = lhs.Nx;
  const int Ny = lhs.Ny;
  //const int Nz = lhs.Nz;

  const int halo = 3; // halo region size for fields

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      if (in == 0 && jn == 0) continue;

      tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn) ));
      if (tpr) {
        auto& rhs = tpr->get_grids();

        /* diagonal rules are:
        if + then to   n
        if + then from 0

        if - then to   -1
        if - then from n-1
        */

        if (in == +1) { ito = Nx; ifro = 0; }
        if (jn == +1) { jto = Ny; jfro = 0; }

        if (in == -1) { ito = -1; ifro = Nx-1; }
        if (jn == -1) { jto = -1; jfro = Ny-1; }

        // copy (halo = 1) assignment
        //if      (jn == 0) copy_vert_grids(    mesh, mpr, ito, ifro);   // vertical
        //else if (in == 0) copy_horz_grids(    mesh, mpr, jto, jfro);   // horizontal
        //else              copy_z_pencil_grids(mesh, mpr, ito, jto, ifro, jfro); // diagonal
        
        // generalized halo >= 1 loops
        if (jn == 0) { // vertical
          //for(int h=0; h<halo; h++) copy_vert_grids(mesh, mpr, iarr, ito+in*h, ifro+in*h);   

          if(has_elem(iarr, 0)) {

            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int j, int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.jx(ito+in*h, j, 0) = rhs_in.jx(ifro+in*h, j, 0);
              lhs_in.jy(ito+in*h, j, 0) = rhs_in.jy(ifro+in*h, j, 0);
              lhs_in.jz(ito+in*h, j, 0) = rhs_in.jz(ifro+in*h, j, 0);
            }, Ny, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D([=] (int j, int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.jx(ito+in*h, j, 0) = rhs_in.jx(ifro+in*h, j, 0);
              lhs_in.jy(ito+in*h, j, 0) = rhs_in.jy(ifro+in*h, j, 0);
              lhs_in.jz(ito+in*h, j, 0) = rhs_in.jz(ifro+in*h, j, 0);
            }, Ny, halo, lhs, rhs);
            #endif
          }

          // TODO rest of the operations
          if(has_elem(iarr, 1)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int j, int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.ex(ito+in*h, j, 0) = rhs_in.ex(ifro+in*h, j, 0);
              lhs_in.ey(ito+in*h, j, 0) = rhs_in.ey(ifro+in*h, j, 0);
              lhs_in.ez(ito+in*h, j, 0) = rhs_in.ez(ifro+in*h, j, 0);
            }, Ny, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int j, int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.ex(ito+in*h, j, 0) = rhs_in.ex(ifro+in*h, j, 0);
              lhs_in.ey(ito+in*h, j, 0) = rhs_in.ey(ifro+in*h, j, 0);
              lhs_in.ez(ito+in*h, j, 0) = rhs_in.ez(ifro+in*h, j, 0);
            }, Ny, halo, lhs, rhs);
            #endif
          }

          if(has_elem(iarr, 2)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int j, int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.bx(ito+in*h, j, 0) = rhs_in.bx(ifro+in*h, j, 0);
              lhs_in.by(ito+in*h, j, 0) = rhs_in.by(ifro+in*h, j, 0);
              lhs_in.bz(ito+in*h, j, 0) = rhs_in.bz(ifro+in*h, j, 0);
            }, Ny, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int j, int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.bx(ito+in*h, j, 0) = rhs_in.bx(ifro+in*h, j, 0);
              lhs_in.by(ito+in*h, j, 0) = rhs_in.by(ifro+in*h, j, 0);
              lhs_in.bz(ito+in*h, j, 0) = rhs_in.bz(ifro+in*h, j, 0);
            }, Ny, halo, lhs, rhs);
            #endif
          }



        } else if (in == 0) { // horizontal
          //for(int g=0; g<halo; g++) copy_horz_grids(mesh, mpr, iarr, jto+jn*g, jfro+jn*g);   

          if(has_elem(iarr, 0)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int i, int g, Grids &lhs_in, Grids &rhs_in){
              lhs_in.jx(i, jto+jn*g, 0) = rhs_in.jx(i, jfro+jn*g, 0);
              lhs_in.jy(i, jto+jn*g, 0) = rhs_in.jy(i, jfro+jn*g, 0);
              lhs_in.jz(i, jto+jn*g, 0) = rhs_in.jz(i, jfro+jn*g, 0);
            }, Nx, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int i, int g, Grids &lhs_in, Grids &rhs_in){
              lhs_in.jx(i, jto+jn*g, 0) = rhs_in.jx(i, jfro+jn*g, 0);
              lhs_in.jy(i, jto+jn*g, 0) = rhs_in.jy(i, jfro+jn*g, 0);
              lhs_in.jz(i, jto+jn*g, 0) = rhs_in.jz(i, jfro+jn*g, 0);
            }, Nx, halo, lhs, rhs);

            #endif
          }
            
          if(has_elem(iarr, 1)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int i, int g, Grids &lhs_in, Grids &rhs_in){
              lhs_in.ex(i, jto+jn*g, 0) = rhs_in.ex(i, jfro+jn*g, 0);
              lhs_in.ey(i, jto+jn*g, 0) = rhs_in.ey(i, jfro+jn*g, 0);
              lhs_in.ez(i, jto+jn*g, 0) = rhs_in.ez(i, jfro+jn*g, 0);
            }, Nx, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int i, int g, Grids &lhs_in, Grids &rhs_in){
              lhs_in.ex(i, jto+jn*g, 0) = rhs_in.ex(i, jfro+jn*g, 0);
              lhs_in.ey(i, jto+jn*g, 0) = rhs_in.ey(i, jfro+jn*g, 0);
              lhs_in.ez(i, jto+jn*g, 0) = rhs_in.ez(i, jfro+jn*g, 0);
            }, Nx, halo, lhs, rhs);

            #endif
          }

          if(has_elem(iarr, 2)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int i, int g, Grids &lhs_in, Grids &rhs_in){
              lhs_in.bx(i, jto+jn*g, 0) = rhs_in.bx(i, jfro+jn*g, 0);
              lhs_in.by(i, jto+jn*g, 0) = rhs_in.by(i, jfro+jn*g, 0);
              lhs_in.bz(i, jto+jn*g, 0) = rhs_in.bz(i, jfro+jn*g, 0);
            }, Nx, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int i, int g, Grids &lhs_in, Grids &rhs_in){
              lhs_in.bx(i, jto+jn*g, 0) = rhs_in.bx(i, jfro+jn*g, 0);
              lhs_in.by(i, jto+jn*g, 0) = rhs_in.by(i, jfro+jn*g, 0);
              lhs_in.bz(i, jto+jn*g, 0) = rhs_in.bz(i, jfro+jn*g, 0);
            }, Nx, halo, lhs, rhs);

            #endif
          }

        } else { // diagonal
          //for(int h=0; h<halo; h++) {
          //  for(int g=0; g<halo; g++) {
          //    copy_z_pencil_grids(mesh, mpr,iarr, ito+in*h, jto+jn*g, ifro+in*h, jfro+jn*g); 
          //  }
          //}

          if(has_elem(iarr, 0)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int g ,int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.jx(ito+in*h, jto+jn*g, 0) = rhs_in.jx(ifro+in*h, jfro+jn*g, 0);
              lhs_in.jy(ito+in*h, jto+jn*g, 0) = rhs_in.jy(ifro+in*h, jfro+jn*g, 0);
              lhs_in.jz(ito+in*h, jto+jn*g, 0) = rhs_in.jz(ifro+in*h, jfro+jn*g, 0);
            }, halo, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int g ,int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.jx(ito+in*h, jto+jn*g, 0) = rhs_in.jx(ifro+in*h, jfro+jn*g, 0);
              lhs_in.jy(ito+in*h, jto+jn*g, 0) = rhs_in.jy(ifro+in*h, jfro+jn*g, 0);
              lhs_in.jz(ito+in*h, jto+jn*g, 0) = rhs_in.jz(ifro+in*h, jfro+jn*g, 0);
            }, halo, halo, lhs, rhs);

            #endif
          }
            
          if(has_elem(iarr, 1)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int g ,int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.ex(ito+in*h, jto+jn*g, 0) = rhs_in.ex(ifro+in*h, jfro+jn*g, 0);
              lhs_in.ey(ito+in*h, jto+jn*g, 0) = rhs_in.ey(ifro+in*h, jfro+jn*g, 0);
              lhs_in.ez(ito+in*h, jto+jn*g, 0) = rhs_in.ez(ifro+in*h, jfro+jn*g, 0);
            }, halo, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int g ,int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.ex(ito+in*h, jto+jn*g, 0) = rhs_in.ex(ifro+in*h, jfro+jn*g, 0);
              lhs_in.ey(ito+in*h, jto+jn*g, 0) = rhs_in.ey(ifro+in*h, jfro+jn*g, 0);
              lhs_in.ez(ito+in*h, jto+jn*g, 0) = rhs_in.ez(ifro+in*h, jfro+jn*g, 0);
            }, halo, halo, lhs, rhs);

            #endif
          }

          if(has_elem(iarr, 2)) {
            #ifdef VEC_FLD2D
            UniIter::iterate2D([=] DEVCALLABLE (int g ,int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.bx(ito+in*h, jto+jn*g, 0) = rhs_in.bx(ifro+in*h, jfro+jn*g, 0);
              lhs_in.by(ito+in*h, jto+jn*g, 0) = rhs_in.by(ifro+in*h, jfro+jn*g, 0);
              lhs_in.bz(ito+in*h, jto+jn*g, 0) = rhs_in.bz(ifro+in*h, jfro+jn*g, 0);
            }, halo, halo, lhs, rhs);
            #else
            UniIter::UniIterHost::iterate2D_nonvec([=] (int g ,int h, Grids &lhs_in, Grids &rhs_in){
              lhs_in.bx(ito+in*h, jto+jn*g, 0) = rhs_in.bx(ifro+in*h, jfro+jn*g, 0);
              lhs_in.by(ito+in*h, jto+jn*g, 0) = rhs_in.by(ifro+in*h, jfro+jn*g, 0);
              lhs_in.bz(ito+in*h, jto+jn*g, 0) = rhs_in.bz(ifro+in*h, jfro+jn*g, 0);
            }, halo, halo, lhs, rhs);

            #endif
          }

        }
      } // end of if(tpr)
    }
  }
}



template<>
void Tile<3>::update_boundaries(
        corgi::Grid<3>& grid,
        std::vector<int> iarr
        )
{
  //std::cout << "upB: updating boundaries\n";
#ifdef GPU
  nvtxRangePush(__FUNCTION__);
#endif

  using Tile_t  = Tile<3>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito=0, jto=0, kto=0, ifro=0, jfro=0, kfro=0;
  Tileptr tpr;

  auto& lhs = get_grids(); // target as a reference to update into
  const int halo = 3; // halo region size for fields

  const int Nx = lhs.Nx;
  const int Ny = lhs.Ny;
  const int Nz = lhs.Nz;

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      for(int kn=-1; kn <= 1; kn++) {

        if (in == 0 && jn == 0 && kn == 0) continue;

        // continue only if the tile exists (get_tileptr return nullptr otherwise)
        tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn, kn) ));
        if (tpr) {
          auto& rhs = tpr->get_grids();

          /* diagonal rules are:
          if + then to   n
          if + then from 0

          if - then to   -1
          if - then from n-1
          */

          if (in == +1) { ito = Nx; ifro = 0; }
          if (jn == +1) { jto = Ny; jfro = 0; }
          if (kn == +1) { kto = Nz; kfro = 0; }

          if (in == -1) { ito = -1; ifro = Nx-1; }
          if (jn == -1) { jto = -1; jfro = Ny-1; }
          if (kn == -1) { kto = -1; kfro = Nz-1; }

          // generalized halo >= 1 loops

          if (kn == 0) {

            // vertical
            if (jn == 0) { 
              //for(int h=0; h<halo; h++) copy_vert_grids(mesh, mpr, iarr, ito+in*h, ifro+in*h);   
              //copy_vert_grids_halo(mesh, mpr, mesh.ex.Ny, mesh.ex.Nz, halo, ito, ifro, in, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito+in*h, j, k) = rhs_in.jx(ifro+in*h, j, k);
                  lhs_in.jy(ito+in*h, j, k) = rhs_in.jy(ifro+in*h, j, k);
                  lhs_in.jz(ito+in*h, j, k) = rhs_in.jz(ifro+in*h, j, k);
                }, Ny, Nz, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito+in*h, j, k) = rhs_in.jx(ifro+in*h, j, k);
                  lhs_in.jy(ito+in*h, j, k) = rhs_in.jy(ifro+in*h, j, k);
                  lhs_in.jz(ito+in*h, j, k) = rhs_in.jz(ifro+in*h, j, k);
                }, Ny, Nz, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito+in*h, j, k) = rhs_in.ex(ifro+in*h, j, k);
                  lhs_in.ey(ito+in*h, j, k) = rhs_in.ey(ifro+in*h, j, k);
                  lhs_in.ez(ito+in*h, j, k) = rhs_in.ez(ifro+in*h, j, k);
                }, Ny, Nz, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito+in*h, j, k) = rhs_in.ex(ifro+in*h, j, k);
                  lhs_in.ey(ito+in*h, j, k) = rhs_in.ey(ifro+in*h, j, k);
                  lhs_in.ez(ito+in*h, j, k) = rhs_in.ez(ifro+in*h, j, k);
                }, Ny, Nz, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito+in*h, j, k) = rhs_in.bx(ifro+in*h, j, k);
                  lhs_in.by(ito+in*h, j, k) = rhs_in.by(ifro+in*h, j, k);
                  lhs_in.bz(ito+in*h, j, k) = rhs_in.bz(ifro+in*h, j, k);
                }, Ny, Nz, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito+in*h, j, k) = rhs_in.bx(ifro+in*h, j, k);
                  lhs_in.by(ito+in*h, j, k) = rhs_in.by(ifro+in*h, j, k);
                  lhs_in.bz(ito+in*h, j, k) = rhs_in.bz(ifro+in*h, j, k);
                }, Ny, Nz, halo, lhs, rhs);

                #endif
              }

            // horizontal
            } else if (in == 0) { 
              //for(int g=0; g<halo; g++) copy_horz_grids(mesh, mpr,iarr, jto+jn*g, jfro+jn*g);   
              //copy_horz_grids_halo(mesh, mpr, mesh.ex.Nx, mesh.ex.Nz, halo, jto, jfro, jn, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(i, jto+jn*g, k) = rhs_in.jx(i, jfro+jn*g, k);
                  lhs_in.jy(i, jto+jn*g, k) = rhs_in.jy(i, jfro+jn*g, k);
                  lhs_in.jz(i, jto+jn*g, k) = rhs_in.jz(i, jfro+jn*g, k);
                }, Nx, Nz, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(i, jto+jn*g, k) = rhs_in.jx(i, jfro+jn*g, k);
                  lhs_in.jy(i, jto+jn*g, k) = rhs_in.jy(i, jfro+jn*g, k);
                  lhs_in.jz(i, jto+jn*g, k) = rhs_in.jz(i, jfro+jn*g, k);
                }, Nx, Nz, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(i, jto+jn*g, k) = rhs_in.ex(i, jfro+jn*g, k);
                  lhs_in.ey(i, jto+jn*g, k) = rhs_in.ey(i, jfro+jn*g, k);
                  lhs_in.ez(i, jto+jn*g, k) = rhs_in.ez(i, jfro+jn*g, k);
                }, Nx, Nz, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(i, jto+jn*g, k) = rhs_in.ex(i, jfro+jn*g, k);
                  lhs_in.ey(i, jto+jn*g, k) = rhs_in.ey(i, jfro+jn*g, k);
                  lhs_in.ez(i, jto+jn*g, k) = rhs_in.ez(i, jfro+jn*g, k);
                }, Nx, Nz, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(i, jto+jn*g, k) = rhs_in.bx(i, jfro+jn*g, k);
                  lhs_in.by(i, jto+jn*g, k) = rhs_in.by(i, jfro+jn*g, k);
                  lhs_in.bz(i, jto+jn*g, k) = rhs_in.bz(i, jfro+jn*g, k);
                }, Nx, Nz, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(i, jto+jn*g, k) = rhs_in.bx(i, jfro+jn*g, k);
                  lhs_in.by(i, jto+jn*g, k) = rhs_in.by(i, jfro+jn*g, k);
                  lhs_in.bz(i, jto+jn*g, k) = rhs_in.bz(i, jfro+jn*g, k);
                }, Nx, Nz, halo, lhs, rhs);

                #endif
              }

            // diagonal
            } else { 
              //for(int h=0; h<halo; h++) {
              //  for(int g=0; g<halo; g++) {
              //    copy_z_pencil_grids(mesh, mpr, iarr, ito+in*h, jto+jn*g, ifro+in*h, jfro+jn*g); 
              //} }
             //copy_z_pencil_grids_halo(mesh, mpr, mesh.ex.Nz, halo, ito, ifro, jto, jfro, in, jn, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito+in*h, jto+jn*g, k) = rhs_in.jx(ifro+in*h, jfro+jn*g, k);
                  lhs_in.jy(ito+in*h, jto+jn*g, k) = rhs_in.jy(ifro+in*h, jfro+jn*g, k);
                  lhs_in.jz(ito+in*h, jto+jn*g, k) = rhs_in.jz(ifro+in*h, jfro+jn*g, k);
                }, Nz, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito+in*h, jto+jn*g, k) = rhs_in.jx(ifro+in*h, jfro+jn*g, k);
                  lhs_in.jy(ito+in*h, jto+jn*g, k) = rhs_in.jy(ifro+in*h, jfro+jn*g, k);
                  lhs_in.jz(ito+in*h, jto+jn*g, k) = rhs_in.jz(ifro+in*h, jfro+jn*g, k);
                }, Nz, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito+in*h, jto+jn*g, k) = rhs_in.ex(ifro+in*h, jfro+jn*g, k);
                  lhs_in.ey(ito+in*h, jto+jn*g, k) = rhs_in.ey(ifro+in*h, jfro+jn*g, k);
                  lhs_in.ez(ito+in*h, jto+jn*g, k) = rhs_in.ez(ifro+in*h, jfro+jn*g, k);
                }, Nz, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito+in*h, jto+jn*g, k) = rhs_in.ex(ifro+in*h, jfro+jn*g, k);
                  lhs_in.ey(ito+in*h, jto+jn*g, k) = rhs_in.ey(ifro+in*h, jfro+jn*g, k);
                  lhs_in.ez(ito+in*h, jto+jn*g, k) = rhs_in.ez(ifro+in*h, jfro+jn*g, k);
                }, Nz, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito+in*h, jto+jn*g, k) = rhs_in.bx(ifro+in*h, jfro+jn*g, k);
                  lhs_in.by(ito+in*h, jto+jn*g, k) = rhs_in.by(ifro+in*h, jfro+jn*g, k);
                  lhs_in.bz(ito+in*h, jto+jn*g, k) = rhs_in.bz(ifro+in*h, jfro+jn*g, k);
                }, Nz, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito+in*h, jto+jn*g, k) = rhs_in.bx(ifro+in*h, jfro+jn*g, k);
                  lhs_in.by(ito+in*h, jto+jn*g, k) = rhs_in.by(ifro+in*h, jfro+jn*g, k);
                  lhs_in.bz(ito+in*h, jto+jn*g, k) = rhs_in.bz(ifro+in*h, jfro+jn*g, k);
                }, Nz, halo, halo, lhs, rhs);

                #endif
              }

            } 
         
          // 3D case with kn != 0
          } else {
            
            // infront/behind directions
            if (in == 0 && jn == 0 && kn != 0) { 
              //for(int g=0; g<halo; g++) copy_face_grids(mesh, mpr, iarr, kto+kn*g, kfro+kn*g);   
              //copy_face_grids_halo(mesh, mpr, mesh.ex.Nx, mesh.ex.Ny, halo, kto, kfro, kn, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int j ,int f, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(i, j, kto +kn*f) =  rhs_in.jx(i, j, kfro+kn*f);
                  lhs_in.jy(i, j, kto +kn*f) =  rhs_in.jy(i, j, kfro+kn*f);
                  lhs_in.jz(i, j, kto +kn*f) =  rhs_in.jz(i, j, kfro+kn*f);
                }, Nx, Ny, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int j ,int f, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(i, j, kto +kn*f) =  rhs_in.jx(i, j, kfro+kn*f);
                  lhs_in.jy(i, j, kto +kn*f) =  rhs_in.jy(i, j, kfro+kn*f);
                  lhs_in.jz(i, j, kto +kn*f) =  rhs_in.jz(i, j, kfro+kn*f);
                }, Nx, Ny, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int j ,int f, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(i, j, kto +kn*f) =  rhs_in.ex(i, j, kfro+kn*f);
                  lhs_in.ey(i, j, kto +kn*f) =  rhs_in.ey(i, j, kfro+kn*f);
                  lhs_in.ez(i, j, kto +kn*f) =  rhs_in.ez(i, j, kfro+kn*f);
                }, Nx, Ny, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int j ,int f, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(i, j, kto +kn*f) =  rhs_in.ex(i, j, kfro+kn*f);
                  lhs_in.ey(i, j, kto +kn*f) =  rhs_in.ey(i, j, kfro+kn*f);
                  lhs_in.ez(i, j, kto +kn*f) =  rhs_in.ez(i, j, kfro+kn*f);
                }, Nx, Ny, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int j ,int f, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(i, j, kto +kn*f) =  rhs_in.bx(i, j, kfro+kn*f);
                  lhs_in.by(i, j, kto +kn*f) =  rhs_in.by(i, j, kfro+kn*f);
                  lhs_in.bz(i, j, kto +kn*f) =  rhs_in.bz(i, j, kfro+kn*f);
                }, Nx, Ny, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int j ,int f, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(i, j, kto +kn*f) =  rhs_in.bx(i, j, kfro+kn*f);
                  lhs_in.by(i, j, kto +kn*f) =  rhs_in.by(i, j, kfro+kn*f);
                  lhs_in.bz(i, j, kto +kn*f) =  rhs_in.bz(i, j, kfro+kn*f);
                }, Nx, Ny, halo, lhs, rhs);

                #endif
              }


            // 3D generalized diagonal locations
            // If the finite-difference scheme is purely non-diagonal
            // then these can be dropped off.
              
            // vertical wedges
            } else if (in != 0 && jn == 0 && kn != 0) {

              // y pencils
              //for(int h=0; h<halo; h++) {
              //for(int g=0; g<halo; g++) {
              //    copy_y_pencil_grids(mesh, mpr, iarr, ito+in*h, kto+kn*g, ifro+in*h, kfro+kn*g); 
              //}}
              //copy_y_pencil_grids_halo(mesh, mpr, mesh.ex.Ny, halo, ito, ifro, kto, kfro, in, kn, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito+in*h, j, kto+kn*g) = rhs_in.jx(ifro+in*h, j, kfro+kn*g);
                  lhs_in.jy(ito+in*h, j, kto+kn*g) = rhs_in.jy(ifro+in*h, j, kfro+kn*g);
                  lhs_in.jz(ito+in*h, j, kto+kn*g) = rhs_in.jz(ifro+in*h, j, kfro+kn*g);
                }, Ny, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito+in*h, j, kto+kn*g) = rhs_in.jx(ifro+in*h, j, kfro+kn*g);
                  lhs_in.jy(ito+in*h, j, kto+kn*g) = rhs_in.jy(ifro+in*h, j, kfro+kn*g);
                  lhs_in.jz(ito+in*h, j, kto+kn*g) = rhs_in.jz(ifro+in*h, j, kfro+kn*g);
                }, Ny, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito+in*h, j, kto+kn*g) = rhs_in.ex(ifro+in*h, j, kfro+kn*g);
                  lhs_in.ey(ito+in*h, j, kto+kn*g) = rhs_in.ey(ifro+in*h, j, kfro+kn*g);
                  lhs_in.ez(ito+in*h, j, kto+kn*g) = rhs_in.ez(ifro+in*h, j, kfro+kn*g);
                }, Ny, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito+in*h, j, kto+kn*g) = rhs_in.ex(ifro+in*h, j, kfro+kn*g);
                  lhs_in.ey(ito+in*h, j, kto+kn*g) = rhs_in.ey(ifro+in*h, j, kfro+kn*g);
                  lhs_in.ez(ito+in*h, j, kto+kn*g) = rhs_in.ez(ifro+in*h, j, kfro+kn*g);
                }, Ny, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito+in*h, j, kto+kn*g) = rhs_in.bx(ifro+in*h, j, kfro+kn*g);
                  lhs_in.by(ito+in*h, j, kto+kn*g) = rhs_in.by(ifro+in*h, j, kfro+kn*g);
                  lhs_in.bz(ito+in*h, j, kto+kn*g) = rhs_in.bz(ifro+in*h, j, kfro+kn*g);
                }, Ny, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito+in*h, j, kto+kn*g) = rhs_in.bx(ifro+in*h, j, kfro+kn*g);
                  lhs_in.by(ito+in*h, j, kto+kn*g) = rhs_in.by(ifro+in*h, j, kfro+kn*g);
                  lhs_in.bz(ito+in*h, j, kto+kn*g) = rhs_in.bz(ifro+in*h, j, kfro+kn*g);
                }, Ny, halo, halo, lhs, rhs);

                #endif
              }


            // horizontal wedges
            } else if (in == 0 && jn != 0 && kn != 0) {

              // x pencils
              //for(int h=0; h<halo; h++) {
              //for(int g=0; g<halo; g++) {
              //    copy_x_pencil_grids(mesh, mpr, iarr, jto+jn*h, kto+kn*g, jfro+jn*h, kfro+kn*g); 
              //}}
              //copy_x_pencil_grids_halo(mesh, mpr, mesh.ex.Nx, halo, jto, jfro, kto, kfro, jn, kn, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(i, jto+jn*h, kto+kn*g) = rhs_in.jx(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.jy(i, jto+jn*h, kto+kn*g) = rhs_in.jy(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.jz(i, jto+jn*h, kto+kn*g) = rhs_in.jz(i, jfro+jn*h, kfro+kn*g);
                }, Nx, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(i, jto+jn*h, kto+kn*g) = rhs_in.jx(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.jy(i, jto+jn*h, kto+kn*g) = rhs_in.jy(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.jz(i, jto+jn*h, kto+kn*g) = rhs_in.jz(i, jfro+jn*h, kfro+kn*g);
                }, Nx, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(i, jto+jn*h, kto+kn*g) = rhs_in.ex(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.ey(i, jto+jn*h, kto+kn*g) = rhs_in.ey(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.ez(i, jto+jn*h, kto+kn*g) = rhs_in.ez(i, jfro+jn*h, kfro+kn*g);
                }, Nx, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(i, jto+jn*h, kto+kn*g) = rhs_in.ex(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.ey(i, jto+jn*h, kto+kn*g) = rhs_in.ey(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.ez(i, jto+jn*h, kto+kn*g) = rhs_in.ez(i, jfro+jn*h, kfro+kn*g);
                }, Nx, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(i, jto+jn*h, kto+kn*g) = rhs_in.bx(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.by(i, jto+jn*h, kto+kn*g) = rhs_in.by(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.bz(i, jto+jn*h, kto+kn*g) = rhs_in.bz(i, jfro+jn*h, kfro+kn*g);
                }, Nx, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(i, jto+jn*h, kto+kn*g) = rhs_in.bx(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.by(i, jto+jn*h, kto+kn*g) = rhs_in.by(i, jfro+jn*h, kfro+kn*g);
                  lhs_in.bz(i, jto+jn*h, kto+kn*g) = rhs_in.bz(i, jfro+jn*h, kfro+kn*g);
                }, Nx, halo, halo, lhs, rhs);

                #endif
              }

            // corners
            } else if (in != 0 && jn != 0 && kn != 0) {

              //pointwise
              //for(int h=0; h<halo; h++) {
              //for(int g=0; g<halo; g++) {
              //for(int f=0; f<halo; f++) {
              //  copy_point_grids(mesh, mpr, iarr,
              //          ito +in*h, jto +jn*g, kto +kn*f,
              //          ifro+in*h, jfro+jn*g, kfro+kn*f);
              //}}}
              //copy_point_grids_halo(mesh, mpr, halo, ito, ifro, jto, jfro, kto, kfro, in, jn, kn, ind);

              if(has_elem(iarr, 0)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jx(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.jy(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jy(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.jz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jz(ifro+in*h, jfro+jn*g, kfro+kn*f);
                }, halo, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.jx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jx(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.jy(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jy(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.jz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jz(ifro+in*h, jfro+jn*g, kfro+kn*f);
                }, halo, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 1)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ex(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.ey(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ey(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.ez(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ez(ifro+in*h, jfro+jn*g, kfro+kn*f);
                }, halo, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.ex(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ex(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.ey(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ey(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.ez(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ez(ifro+in*h, jfro+jn*g, kfro+kn*f);
                }, halo, halo, halo, lhs, rhs);

                #endif
              }

              if(has_elem(iarr, 2)) {
                #ifdef VEC_FLD3D
                UniIter::iterate3D([=] DEVCALLABLE (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.bx(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.by(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.by(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.bz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.bz(ifro+in*h, jfro+jn*g, kfro+kn*f);
                }, halo, halo, halo, lhs, rhs);
                #else
                UniIter::UniIterHost::iterate3D_nonvec([=] (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                  lhs_in.bx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.bx(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.by(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.by(ifro+in*h, jfro+jn*g, kfro+kn*f);
                  lhs_in.bz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.bz(ifro+in*h, jfro+jn*g, kfro+kn*f);
                }, halo, halo, halo, lhs, rhs);

                #endif
              }


            }

          } // 3D cases with kn != 0
        } // if tpr
      } // kn
    } // jn
  } // in

  UniIter::sync();
  
#ifdef GPU
  nvtxRangePop();
#endif

}


template<>
void Tile<1>::exchange_currents(corgi::Grid<1>& grid) 
{

  using Tile_t  = Tile<1>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito=0, ifro=0;
  Tileptr tpr; 

  const int halo = 3; // halo region size for currents

  auto& lhs = get_grids(); // target as a reference to update into
  const int Nx = lhs.Nx;

  // add from right side to left
  //for(int h=1; h<= halo; h++) add_vert_grids(mesh, mleft, -h, mleft.Nx-h); 
  // add from left side to right
  //for(int h=1; h<= halo; h++) add_vert_grids(mesh, mright, mesh.Nx+h-1, h-1); 

  for(int in=-1; in <= 1; in++) {
    if (in == 0) continue;
    tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in) ));

    if (tpr) {
      auto& rhs = tpr->get_grids();

      //shifted neg's from vals with -1
      if (in == +1) { ito = Nx-1; ifro = -1; }
      if (in == -1) { ito = 0;    ifro = Nx; }

        UniIter::iterate([=] DEVCALLABLE (int h, Grids &lhs_in, Grids &rhs_in){  
          atomic_add( lhs_in.jx(ito-in*h, 0, 0), rhs_in.jx(ifro-in*h, 0, 0) );
          atomic_add( lhs_in.jy(ito-in*h, 0, 0), rhs_in.jy(ifro-in*h, 0, 0) );
          atomic_add( lhs_in.jz(ito-in*h, 0, 0), rhs_in.jz(ifro-in*h, 0, 0) );
        }, halo, lhs, rhs);

    }
  }
}



/// Update currents on Yee grid boundaries
//
// The whole "FROM" -> "TO" index selection depending on neighbor is abstractified.
// The rules are:
//  - if neighbor is - (i.e., left, or bottom) => TO=-1 & FROM=N-1
//  - if neighbor is + (i.e., right or top)    => TO=N  & FROM=0
//  
// Therefore, given a range h=1,2,3,..halo, we need to copy/add values 
//  to an index TO-S*h from index FRO-S*h
// where S is the sign of the neighbor.
//
//--------------------------------------------------
// Here is an example from an addition of top tile (i.e. +1 neighbor)
//
// local index in mesh:  | neighbors index:
//  (outside mesh)  Ny   | 0  (start of neighbor mesh; values below are halo regions)
//                -------|-------
//                  Ny-1 | -1
//                  Ny-2 | -2
//                  Ny-3 | -3
//
// so we need to add:
//  neighbor at j=-1 to j=Ny-1
//  neighbor at j=-2 to j=Ny-2
//  neighbor at j=-3 to j=Ny-3
//
// In a loop over h=1,2,3 this is, given more succinctly:
//  FRO - S*h into TO - S*h, 
//
//  where FRO=Ny, TO=0, and S=-1.
//--------------------------------------------------
//
template<>
void Tile<2>::exchange_currents(corgi::Grid<2>& grid) 
{

  using Tile_t  = Tile<2>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito=0, jto=0, ifro=0, jfro=0;
  Tileptr tpr; 

  const int halo = 3;

  auto& lhs = get_grids(); // target as a reference to update into

  const int Nx = lhs.Nx;
  const int Ny = lhs.Ny;
  //const int Nz = lhs.Nz;

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      if (in == 0 && jn == 0) continue;

      tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn) ));
      if (tpr) {
        auto& rhs = tpr->get_grids();

        /* diagonal rules are:
        if + then to   n
        if + then from 0

        if - then to   -1
        if - then from n-1
        */

        //shifted neg's from vals with -1
        if (in == +1) { ito = Nx-1; ifro = -1; }
        if (jn == +1) { jto = Ny-1; jfro = -1; }

        if (in == -1) { ito = 0;    ifro = Nx; }
        if (jn == -1) { jto = 0;    jfro = Ny; }

        // add
        if (jn == 0) { // vertical
          //for(int h=1; h<=halo; h++) add_vert_grids(mesh, mpr, ito-in*h, ifro-in*h);   
          
          #ifdef VEC_CUR2D
          UniIter::iterate2D([=] DEVCALLABLE (int j,int h, Grids &lhs_in, Grids &rhs_in){  
            atomic_add( lhs_in.jx(ito-in*h, j, 0), rhs_in.jx(ifro-in*h, j, 0) );
            atomic_add( lhs_in.jy(ito-in*h, j, 0), rhs_in.jy(ifro-in*h, j, 0) );
            atomic_add( lhs_in.jz(ito-in*h, j, 0), rhs_in.jz(ifro-in*h, j, 0) );
          }, Ny, halo, lhs, rhs);
          #else
          UniIter::UniIterHost::iterate2D_nonvec([=] (int j,int h, Grids &lhs_in, Grids &rhs_in){  
            lhs_in.jx(ito-in*h, j, 0) += rhs_in.jx(ifro-in*h, j, 0);
            lhs_in.jy(ito-in*h, j, 0) += rhs_in.jy(ifro-in*h, j, 0);
            lhs_in.jz(ito-in*h, j, 0) += rhs_in.jz(ifro-in*h, j, 0);
          }, Ny, halo, lhs, rhs);
          #endif

        } else if (in == 0) { // horizontal
          //for(int g=1; g<=halo; g++) add_horz_grids(mesh, mpr, jto-jn*g, jfro-jn*g);   

          #ifdef VEC_CUR2D
          UniIter::iterate2D([=] DEVCALLABLE (int i, int g, Grids &lhs_in, Grids &rhs_in){
            atomic_add( lhs_in.jx(i, jto-jn*g, 0), rhs_in.jx(i, jfro-jn*g, 0) );
            atomic_add( lhs_in.jy(i, jto-jn*g, 0), rhs_in.jy(i, jfro-jn*g, 0) );
            atomic_add( lhs_in.jz(i, jto-jn*g, 0), rhs_in.jz(i, jfro-jn*g, 0) );
          }, Nx, halo, lhs, rhs);
          #else
          UniIter::UniIterHost::iterate2D_nonvec([=] (int i, int g, Grids &lhs_in, Grids &rhs_in){
            lhs_in.jx(i, jto-jn*g, 0) += rhs_in.jx(i, jfro-jn*g, 0);
            lhs_in.jy(i, jto-jn*g, 0) += rhs_in.jy(i, jfro-jn*g, 0);
            lhs_in.jz(i, jto-jn*g, 0) += rhs_in.jz(i, jfro-jn*g, 0);
          }, Nx, halo, lhs, rhs);
          #endif
        } else { // diagonal
          //for(int h=1; h<=halo; h++) {
          //  for(int g=1; g<=halo; g++) {
          //    add_z_pencil_grids(mesh, mpr, ito-in*h, jto-jn*g, ifro-in*h, jfro-jn*g); 
          //  }
          //}

          #ifdef VEC_CUR2D
          UniIter::iterate2D([=] DEVCALLABLE (int g ,int h, Grids &lhs_in, Grids &rhs_in){
            atomic_add( lhs_in.jx(ito-in*h, jto-jn*g, 0), rhs_in.jx(ifro-in*h, jfro-jn*g, 0));
            atomic_add( lhs_in.jy(ito-in*h, jto-jn*g, 0), rhs_in.jy(ifro-in*h, jfro-jn*g, 0));
            atomic_add( lhs_in.jz(ito-in*h, jto-jn*g, 0), rhs_in.jz(ifro-in*h, jfro-jn*g, 0));
          }, halo, halo, lhs, rhs);
          #else
          UniIter::UniIterHost::iterate2D_nonvec([=] (int g ,int h, Grids &lhs_in, Grids &rhs_in){
            lhs_in.jx(ito-in*h, jto-jn*g, 0) += rhs_in.jx(ifro-in*h, jfro-jn*g, 0);
            lhs_in.jy(ito-in*h, jto-jn*g, 0) += rhs_in.jy(ifro-in*h, jfro-jn*g, 0);
            lhs_in.jz(ito-in*h, jto-jn*g, 0) += rhs_in.jz(ifro-in*h, jfro-jn*g, 0);
          }, halo, halo, lhs, rhs);
          #endif
        }
      } // end of if(tpr)
    }
  }
}

template<>
void Tile<3>::exchange_currents(corgi::Grid<3>& grid) 
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  using Tile_t  = Tile<3>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito=0, jto=0, kto=0, ifro=0, jfro=0, kfro=0;
  Tileptr tpr; 

  const int halo = 3;

  auto& lhs = get_grids(); // target as a reference to update into
  const int Nx = lhs.Nx;
  const int Ny = lhs.Ny;
  const int Nz = lhs.Nz;

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      for(int kn=-1; kn <= 1; kn++) {

        if (in == 0 && jn == 0 && kn == 0) continue;

        tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn, kn) ));
        if (tpr) {
          auto& rhs = tpr->get_grids();

          /* diagonal rules are:
          if + then to   n
          if + then from 0

          if - then to   -1
          if - then from n-1
          */

          //shifted neg's from vals with -1
          if (in == +1) { ito = Nx-1; ifro = -1; }
          if (jn == +1) { jto = Ny-1; jfro = -1; }
          if (kn == +1) { kto = Nz-1; kfro = -1; }

          if (in == -1) { ito = 0;    ifro = Nx; }
          if (jn == -1) { jto = 0;    jfro = Ny; }
          if (kn == -1) { kto = 0;    kfro = Nz; }

          // generalized halo >= 1 loops

          // 2D case; kn == 0; no tiles behind or infront
          if (kn == 0) {

            // vertical
            if (jn == 0) { 
              //for(int h=0; h<halo; h++) add_vert_grids(mesh, mpr, ito-in*h, ifro-in*h);   

              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){  
                atomic_add( lhs_in.jx(ito-in*h, j, k), rhs_in.jx(ifro-in*h, j, k) );
                atomic_add( lhs_in.jy(ito-in*h, j, k), rhs_in.jy(ifro-in*h, j, k) );
                atomic_add( lhs_in.jz(ito-in*h, j, k), rhs_in.jz(ifro-in*h, j, k) );
              }, Ny, Nz, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int k ,int h, Grids &lhs_in, Grids &rhs_in){  
                lhs_in.jx(ito-in*h, j, k) += rhs_in.jx(ifro-in*h, j, k);
                lhs_in.jy(ito-in*h, j, k) += rhs_in.jy(ifro-in*h, j, k);
                lhs_in.jz(ito-in*h, j, k) += rhs_in.jz(ifro-in*h, j, k);
              }, Ny, Nz, halo, lhs, rhs);
              #endif

            // horizontal
            } else if (in == 0) { 
              //for(int g=0; g<halo; g++) add_horz_grids(mesh, mpr, jto-jn*g, jfro-jn*g);   
                
              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                atomic_add( lhs_in.jx(i, jto-jn*g, k), rhs_in.jx(i, jfro-jn*g, k) );
                atomic_add( lhs_in.jy(i, jto-jn*g, k), rhs_in.jy(i, jfro-jn*g, k) );
                atomic_add( lhs_in.jz(i, jto-jn*g, k), rhs_in.jz(i, jfro-jn*g, k) );
              }, Nx, Nz, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int k ,int g, Grids &lhs_in, Grids &rhs_in){
                lhs_in.jx(i, jto-jn*g, k) += rhs_in.jx(i, jfro-jn*g, k);
                lhs_in.jy(i, jto-jn*g, k) += rhs_in.jy(i, jfro-jn*g, k);
                lhs_in.jz(i, jto-jn*g, k) += rhs_in.jz(i, jfro-jn*g, k);
              }, Nx, Nz, halo, lhs, rhs);
              #endif

            // diagonal
            } else { 
              //for(int h=0; h<halo; h++) {
              //  for(int g=0; g<halo; g++) {
              //    add_z_pencil_grids(mesh, mpr, ito-in*h, jto-jn*g, ifro-in*h, jfro-jn*g); 
              //} }

              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                atomic_add( lhs_in.jx(ito-in*h, jto-jn*g, k), rhs_in.jx(ifro-in*h, jfro-jn*g, k) );
                atomic_add( lhs_in.jy(ito-in*h, jto-jn*g, k), rhs_in.jy(ifro-in*h, jfro-jn*g, k) );
                atomic_add( lhs_in.jz(ito-in*h, jto-jn*g, k), rhs_in.jz(ifro-in*h, jfro-jn*g, k) );
              }, Nz, halo, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int k, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                lhs_in.jx(ito-in*h, jto-jn*g, k) += rhs_in.jx(ifro-in*h, jfro-jn*g, k);
                lhs_in.jy(ito-in*h, jto-jn*g, k) += rhs_in.jy(ifro-in*h, jfro-jn*g, k);
                lhs_in.jz(ito-in*h, jto-jn*g, k) += rhs_in.jz(ifro-in*h, jfro-jn*g, k);
              }, Nz, halo, halo, lhs, rhs);
              #endif
            } 

          // 3D case with kn != 0
          } else {
            
            // infront/behind directions
            if (in == 0 && jn == 0 && kn != 0) { 
              //for(int g=0; g<halo; g++) add_face_grids(mesh, mpr, kto-kn*g, kfro-kn*g);   
              //add_face_grids_halo(mesh, mpr, Nx, Ny, halo, kto, kfro, kn, ind);

              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int i, int j ,int g, Grids &lhs_in, Grids &rhs_in){ 
                  atomic_add( lhs_in.jx(i, j, kto-kn*g), rhs_in.jx(i, j, kfro-kn*g) );
                  atomic_add( lhs_in.jy(i, j, kto-kn*g), rhs_in.jy(i, j, kfro-kn*g) );
                  atomic_add( lhs_in.jz(i, j, kto-kn*g), rhs_in.jz(i, j, kfro-kn*g) );
              }, Nx, Ny, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int j ,int g, Grids &lhs_in, Grids &rhs_in){ 
                  lhs_in.jx(i, j, kto-kn*g) += rhs_in.jx(i, j, kfro-kn*g);
                  lhs_in.jy(i, j, kto-kn*g) += rhs_in.jy(i, j, kfro-kn*g);
                  lhs_in.jz(i, j, kto-kn*g) += rhs_in.jz(i, j, kfro-kn*g);
              }, Nx, Ny, halo, lhs, rhs);
              #endif


            //// 3D generalized diagonal locations
            //// If the finite-difference scheme is purely non-diagonal
            //// then these can be dropped off.
                
            //// vertical wedges
            } else if (in != 0 && jn == 0 && kn != 0) {

              // y pencils
              //for(int h=0; h<halo; h++) {
              //for(int g=0; g<halo; g++) {
              //  add_y_pencil_grids(mesh, mpr, ito-in*h, kto-kn*g, ifro-in*h, kfro-kn*g); 
              //}}
              //add_y_pencil_grids_halo(mesh, mpr, Ny, halo, ito, ifro, kto, kfro, in, kn, ind);

              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                atomic_add( lhs_in.jx(ito-in*h, j, kto-kn*g), rhs_in.jx(ifro-in*h, j, kfro-kn*g) );
                atomic_add( lhs_in.jy(ito-in*h, j, kto-kn*g), rhs_in.jy(ifro-in*h, j, kfro-kn*g) );
                atomic_add( lhs_in.jz(ito-in*h, j, kto-kn*g), rhs_in.jz(ifro-in*h, j, kfro-kn*g) );
              }, Ny, halo, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int j, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                lhs_in.jx(ito-in*h, j, kto-kn*g) += rhs_in.jx(ifro-in*h, j, kfro-kn*g);
                lhs_in.jy(ito-in*h, j, kto-kn*g) += rhs_in.jy(ifro-in*h, j, kfro-kn*g);
                lhs_in.jz(ito-in*h, j, kto-kn*g) += rhs_in.jz(ifro-in*h, j, kfro-kn*g);
              }, Ny, halo, halo, lhs, rhs);
              #endif

            } else if (in == 0 && jn != 0 && kn != 0) {
              // x pencils
              //for(int h=0; h<halo; h++) {
              //for(int g=0; g<halo; g++) {
              //  add_x_pencil_grids(mesh, mpr, jto-jn*h, kto-kn*g, jfro-jn*h, kfro-kn*g); 
              //}}
              //add_x_pencil_grids_halo(lhs, rhs, Nx, halo, jto, jfro, kto, kfro, jn, kn, ind);

              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                atomic_add( lhs_in.jx(i, jto-jn*h, kto-kn*g), rhs_in.jx(i, jfro-jn*h, kfro-kn*g) );
                atomic_add( lhs_in.jy(i, jto-jn*h, kto-kn*g), rhs_in.jy(i, jfro-jn*h, kfro-kn*g) );
                atomic_add( lhs_in.jz(i, jto-jn*h, kto-kn*g), rhs_in.jz(i, jfro-jn*h, kfro-kn*g) );
              }, Nx, halo, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int i, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                lhs_in.jx(i, jto-jn*h, kto-kn*g) += rhs_in.jx(i, jfro-jn*h, kfro-kn*g);
                lhs_in.jy(i, jto-jn*h, kto-kn*g) += rhs_in.jy(i, jfro-jn*h, kfro-kn*g);
                lhs_in.jz(i, jto-jn*h, kto-kn*g) += rhs_in.jz(i, jfro-jn*h, kfro-kn*g);
              }, Nx, halo, halo, lhs, rhs);
              #endif

            //// corners
            } else if (in != 0 && jn != 0 && kn != 0) {

              // pointwise
              //for(int h=0; h<halo; h++) {
              //for(int g=0; g<halo; g++) {
              //for(int f=0; f<halo; f++) {
              //  add_point_grids(mesh, mpr, 
              //          ito -in*h, jto -jn*g, kto -kn*f,
              //          ifro-in*h, jfro-jn*g, kfro-kn*f);
              //}}}
              //add_point_grids_halo(mesh, mpr, halo, ito, ifro, jto, jfro, kto, kfro, in, jn, kn, ind);
                
              #ifdef VEC_CUR3D
              UniIter::iterate3D([=] DEVCALLABLE (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                atomic_add( lhs_in.jx(ito -in*h, jto -jn*g, kto -kn*f), rhs_in.jx(ifro-in*h, jfro-jn*g, kfro-kn*f) );
                atomic_add( lhs_in.jy(ito -in*h, jto -jn*g, kto -kn*f), rhs_in.jy(ifro-in*h, jfro-jn*g, kfro-kn*f) );
                atomic_add( lhs_in.jz(ito -in*h, jto -jn*g, kto -kn*f), rhs_in.jz(ifro-in*h, jfro-jn*g, kfro-kn*f) );
              }, halo, halo, halo, lhs, rhs);
              #else
              UniIter::UniIterHost::iterate3D_nonvec([=] (int f, int g ,int h, Grids &lhs_in, Grids &rhs_in){
                lhs_in.jx(ito -in*h, jto -jn*g, kto -kn*f) += rhs_in.jx(ifro-in*h, jfro-jn*g, kfro-kn*f);
                lhs_in.jy(ito -in*h, jto -jn*g, kto -kn*f) += rhs_in.jy(ifro-in*h, jfro-jn*g, kfro-kn*f);
                lhs_in.jz(ito -in*h, jto -jn*g, kto -kn*f) += rhs_in.jz(ifro-in*h, jfro-jn*g, kfro-kn*f);
              }, halo, halo, halo, lhs, rhs);
              #endif
            } 

            //ind++;
          } // 3D cases with kn != 0


        } // if tpr
      }
    }
  }
  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}


template<std::size_t D>
void Tile<D>::clear_current() 
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& gs = this->get_grids();
  gs.jx.clear();
  gs.jy.clear();
  gs.jz.clear();


#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// MPI routines

// create MPI tag given tile id and extra layer of differentiation
int get_tag(int tag, int extra_param)
{
  assert(extra_param <= 9); // max 9 different modes
  assert(tag < (pow(2,16) - 1)); // cray-mpich supports maximum of 2^22-1 tag value

  return tag + (extra_param)*pow(2,16);
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::send_data( 
    mpi::communicator& comm, 
    int dest, 
    int mode,
    int tag)
{

#ifdef GPU
  nvtxRangePush(__FUNCTION__);
#endif

  auto& gs = get_grids(); 
  std::vector<mpi::request> reqs;

  UniIter::sync();

  if (mode == 0) {
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 0), gs.jx.data(), gs.jx.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 1), gs.jy.data(), gs.jy.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 2), gs.jz.data(), gs.jz.size()) );
  } else if (mode == 1) {
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 3), gs.ex.data(), gs.ex.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 4), gs.ey.data(), gs.ey.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 5), gs.ez.data(), gs.ez.size()) );
  } else if (mode == 2) {
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 6), gs.bx.data(), gs.bx.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 7), gs.by.data(), gs.by.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 8), gs.bz.data(), gs.bz.size()) );
  }

#ifdef GPU
  nvtxRangePop();
#endif

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_data( 
    mpi::communicator& comm, 
    int orig, 
    int mode,
    int tag)
{

#ifdef GPU
  nvtxRangePush(__FUNCTION__);
#endif

  auto& gs = get_grids(); 

  std::vector<mpi::request> reqs;

  UniIter::sync();

  if (mode == 0) {
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 0), gs.jx.data(), gs.jx.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 1), gs.jy.data(), gs.jy.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 2), gs.jz.data(), gs.jz.size()) );
  } else if (mode == 1) {
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 3), gs.ex.data(), gs.ex.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 4), gs.ey.data(), gs.ey.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 5), gs.ez.data(), gs.ez.size()) );
  } else if (mode == 2) {
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 6), gs.bx.data(), gs.bx.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 7), gs.by.data(), gs.by.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 8), gs.bz.data(), gs.bz.size()) );
  }

#ifdef GPU
  nvtxRangePop();
#endif

  return reqs;
}

//--------------------------------------------------
// explicit template instantiation

template class Tile<1>;
template class Tile<2>;
template class Tile<3>;

} // end of ns emf
