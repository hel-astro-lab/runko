#include <iostream>
#include <cmath>

#include "tile.h"
//#include <nvtx3/nvToolsExt.h> 


#include "../tools/iter/iter.h"
#include "../tools/iter/allocator.h"

namespace fields {
  using namespace mpi4cpp;


/// Deposit current into electric field
template<std::size_t D>
void Tile<D>::deposit_current() 
{
  YeeLattice& mesh = get_yee();

  mesh.ex -= mesh.jx;
  mesh.ey -= mesh.jy;
  mesh.ez -= mesh.jz;

}


/// Get current time snapshot of Yee lattice
template<std::size_t D>
YeeLattice& Tile<D>::get_yee(int /*i*/) 
{
  //return this->yee.at(0);
  return this->yee;
}
  
/// Get current time snapshot of Yee lattice
template<std::size_t D>
YeeLattice& Tile<D>::get_yee2() 
{
  //return this->yee[0];
  return this->yee;
}

/// Set current time snapshot of Yee lattice
template<std::size_t D>
void Tile<D>::set_yee(YeeLattice& val)
{
  //this->yee[0] = val;
  this->yee = val;
}

template<std::size_t D>
std::shared_ptr<YeeLattice> Tile<D>::get_yeeptr() 
{
  //return std::shared_ptr<YeeLattice>(&this->yee[0]);
  return std::shared_ptr<YeeLattice>(&yee);
}

template<std::size_t D>
const YeeLattice& Tile<D>::get_const_yee(int /*i*/) const 
{
  //return this->yee.at(0);
  return this->yee;
}




//--------------------------------------------------
// Specialize Yee Lattice insertion
template<>
void Tile<1>::add_yee_lattice() 
{
  //yee.emplace_back( mesh_lengths[0], 1, 1);
}

template<>
void Tile<2>::add_yee_lattice() 
{
  //yee.emplace_back( mesh_lengths[0], mesh_lengths[1], 1);
  //yee = YeeLattice( mesh_lengths[0], mesh_lengths[1], 1);
}

template<>
void Tile<3>::add_yee_lattice() 
{
  //yee.emplace_back( mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]);
  //yee = YeeLattice( mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]);
}

//--------------------------------------------------

/// Quick helper function to copy everything inside Yee lattice 
void copy_vert_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int rhsI) 
{

  lhs.ex.copy_vert(rhs.ex, lhsI, rhsI); 
  lhs.ey.copy_vert(rhs.ey, lhsI, rhsI); 
  lhs.ez.copy_vert(rhs.ez, lhsI, rhsI); 

  lhs.bx.copy_vert(rhs.bx, lhsI, rhsI); 
  lhs.by.copy_vert(rhs.by, lhsI, rhsI); 
  lhs.bz.copy_vert(rhs.bz, lhsI, rhsI); 

  //TODO: separate into own function
  lhs.jx.copy_vert(rhs.jx, lhsI, rhsI); 
  lhs.jy.copy_vert(rhs.jy, lhsI, rhsI); 
  lhs.jz.copy_vert(rhs.jz, lhsI, rhsI); 
}


/// Quick helper function to add everything inside Yee lattice 
void add_vert_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int rhsI) 
{
  lhs.jx.add_vert(rhs.jx, lhsI, rhsI); 
  lhs.jy.add_vert(rhs.jy, lhsI, rhsI); 
  lhs.jz.add_vert(rhs.jz, lhsI, rhsI); 

}



/// Quick helper function to copy everything inside Yee lattice 
void copy_horz_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsJ, int rhsJ) 
{
  lhs.ex.copy_horz(rhs.ex, lhsJ, rhsJ); 
  lhs.ey.copy_horz(rhs.ey, lhsJ, rhsJ); 
  lhs.ez.copy_horz(rhs.ez, lhsJ, rhsJ); 
                                    
  lhs.bx.copy_horz(rhs.bx, lhsJ, rhsJ); 
  lhs.by.copy_horz(rhs.by, lhsJ, rhsJ); 
  lhs.bz.copy_horz(rhs.bz, lhsJ, rhsJ); 

  //TODO: separate into own function
  lhs.jx.copy_horz(rhs.jx, lhsJ, rhsJ); 
  lhs.jy.copy_horz(rhs.jy, lhsJ, rhsJ); 
  lhs.jz.copy_horz(rhs.jz, lhsJ, rhsJ); 
}                       


/// Quick helper function to add everything inside Yee lattice 
void add_horz_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsJ, int rhsJ) 
{
  lhs.jx.add_horz(rhs.jx, lhsJ, rhsJ); 
  lhs.jy.add_horz(rhs.jy, lhsJ, rhsJ); 
  lhs.jz.add_horz(rhs.jz, lhsJ, rhsJ); 
}

/// Quick helper function to copy everything inside Yee lattice 
void copy_face_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsK, int rhsK) 
{
  lhs.ex.copy_face(rhs.ex, lhsK, rhsK); 
  lhs.ey.copy_face(rhs.ey, lhsK, rhsK); 
  lhs.ez.copy_face(rhs.ez, lhsK, rhsK); 
                                    
  lhs.bx.copy_face(rhs.bx, lhsK, rhsK); 
  lhs.by.copy_face(rhs.by, lhsK, rhsK); 
  lhs.bz.copy_face(rhs.bz, lhsK, rhsK); 

  //TODO: separate into own function
  lhs.jx.copy_face(rhs.jx, lhsK, rhsK); 
  lhs.jy.copy_face(rhs.jy, lhsK, rhsK); 
  lhs.jz.copy_face(rhs.jz, lhsK, rhsK); 
}


/// Quick helper function to add everything inside Yee lattice 
void add_face_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsK, int rhsK) 
{
  lhs.jx.add_face(rhs.jx, lhsK, rhsK); 
  lhs.jy.add_face(rhs.jy, lhsK, rhsK); 
  lhs.jz.add_face(rhs.jz, lhsK, rhsK); 
}


void copy_x_pencil_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsJ, int lhsK,
    int rhsJ, int rhsK) 
{
  lhs.ex.copy_x_pencil(rhs.ex, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.ey.copy_x_pencil(rhs.ey, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.ez.copy_x_pencil(rhs.ez, lhsJ, lhsK, rhsJ, rhsK); 
                                                     
  lhs.bx.copy_x_pencil(rhs.bx, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.by.copy_x_pencil(rhs.by, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.bz.copy_x_pencil(rhs.bz, lhsJ, lhsK, rhsJ, rhsK); 
                                                     
  //TODO: separate into own function                 
  lhs.jx.copy_x_pencil(rhs.jx, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.jy.copy_x_pencil(rhs.jy, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.jz.copy_x_pencil(rhs.jz, lhsJ, lhsK, rhsJ, rhsK); 
}

void add_x_pencil_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsJ, int lhsK,
    int rhsJ, int rhsK) 
{
  lhs.jx.add_x_pencil(rhs.jx, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.jy.add_x_pencil(rhs.jy, lhsJ, lhsK, rhsJ, rhsK); 
  lhs.jz.add_x_pencil(rhs.jz, lhsJ, lhsK, rhsJ, rhsK); 
}


void copy_y_pencil_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int lhsK,
    int rhsI, int rhsK) 
{           
  lhs.ex.copy_y_pencil(rhs.ex, lhsI, lhsK, rhsI, rhsK); 
  lhs.ey.copy_y_pencil(rhs.ey, lhsI, lhsK, rhsI, rhsK); 
  lhs.ez.copy_y_pencil(rhs.ez, lhsI, lhsK, rhsI, rhsK); 
                                                     
  lhs.bx.copy_y_pencil(rhs.bx, lhsI, lhsK, rhsI, rhsK); 
  lhs.by.copy_y_pencil(rhs.by, lhsI, lhsK, rhsI, rhsK); 
  lhs.bz.copy_y_pencil(rhs.bz, lhsI, lhsK, rhsI, rhsK); 
                                                     
  //TODO: separate into own function                 
  lhs.jx.copy_y_pencil(rhs.jx, lhsI, lhsK, rhsI, rhsK); 
  lhs.jy.copy_y_pencil(rhs.jy, lhsI, lhsK, rhsI, rhsK); 
  lhs.jz.copy_y_pencil(rhs.jz, lhsI, lhsK, rhsI, rhsK); 
}

void add_y_pencil_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int lhsK,
    int rhsI, int rhsK) 
{
  lhs.jx.add_y_pencil(rhs.jx, lhsI, lhsK, rhsI, rhsK); 
  lhs.jy.add_y_pencil(rhs.jy, lhsI, lhsK, rhsI, rhsK); 
  lhs.jz.add_y_pencil(rhs.jz, lhsI, lhsK, rhsI, rhsK); 
}


void copy_z_pencil_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int lhsJ,
    int rhsI, int rhsJ) 
{
  lhs.ex.copy_z_pencil(rhs.ex, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.ey.copy_z_pencil(rhs.ey, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.ez.copy_z_pencil(rhs.ez, lhsI, lhsJ, rhsI, rhsJ); 

  lhs.bx.copy_z_pencil(rhs.bx, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.by.copy_z_pencil(rhs.by, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.bz.copy_z_pencil(rhs.bz, lhsI, lhsJ, rhsI, rhsJ); 

  //TODO: separate into own function
  lhs.jx.copy_z_pencil(rhs.jx, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.jy.copy_z_pencil(rhs.jy, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.jz.copy_z_pencil(rhs.jz, lhsI, lhsJ, rhsI, rhsJ); 
}

void add_z_pencil_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int lhsJ,
    int rhsI, int rhsJ) 
{
  lhs.jx.add_z_pencil(rhs.jx, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.jy.add_z_pencil(rhs.jy, lhsI, lhsJ, rhsI, rhsJ); 
  lhs.jz.add_z_pencil(rhs.jz, lhsI, lhsJ, rhsI, rhsJ); 
}


template<class T, int H>
inline void copy_point(
    toolbox::Mesh<T,H>& lhs, 
    toolbox::Mesh<T,H>& rhs, 
    int lhsI, int lhsJ, int lhsK,
    int rhsI, int rhsJ, int rhsK) 
{
  lhs(lhsI, lhsJ, lhsK) = rhs(rhsI, rhsJ, rhsK);
}

inline void copy_point_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int lhsJ, int lhsK,
    int rhsI, int rhsJ, int rhsK) 
{
  copy_point(lhs.ex, rhs.ex, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  copy_point(lhs.ey, rhs.ey, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  copy_point(lhs.ez, rhs.ez, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);

  copy_point(lhs.bx, rhs.bx, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  copy_point(lhs.by, rhs.by, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  copy_point(lhs.bz, rhs.bz, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);

  copy_point(lhs.jx, rhs.jx, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  copy_point(lhs.jy, rhs.jy, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  copy_point(lhs.jz, rhs.jz, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);

}


template<class T, int H>
inline void add_point(
    toolbox::Mesh<T,H>& lhs, 
    toolbox::Mesh<T,H>& rhs, 
    int lhsI, int lhsJ, int lhsK,
    int rhsI, int rhsJ, int rhsK) 
{
  lhs(lhsI, lhsJ, lhsK) += rhs(rhsI, rhsJ, rhsK);
}

inline void add_point_yee(
    YeeLattice& lhs, 
    YeeLattice& rhs, 
    int lhsI, int lhsJ, int lhsK,
    int rhsI, int rhsJ, int rhsK) 
{
  add_point(lhs.jx, rhs.jx, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  add_point(lhs.jy, rhs.jy, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
  add_point(lhs.jz, rhs.jz, lhsI, lhsJ, lhsK, rhsI, rhsJ, rhsK);
}


// further collapsed helper functions 


void copy_vert_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int Ny, int Nz, int halo, int ito, int ifro, int in, int ind=0)
{
    //
    UniIter::iterate3D([=] DEVCALLABLE (int j, int k ,int h, YeeLattice &lhs_in, YeeLattice &rhs_in){
      //
      lhs_in.ex(ito+in*h, j, k) = rhs_in.ex(ifro+in*h, j, k);
      lhs_in.ey(ito+in*h, j, k) = rhs_in.ey(ifro+in*h, j, k);
      lhs_in.ez(ito+in*h, j, k) = rhs_in.ez(ifro+in*h, j, k);
    
      lhs_in.bx(ito+in*h, j, k) = rhs_in.bx(ifro+in*h, j, k);
      lhs_in.by(ito+in*h, j, k) = rhs_in.by(ifro+in*h, j, k);
      lhs_in.bz(ito+in*h, j, k) = rhs_in.bz(ifro+in*h, j, k);
    
      lhs_in.jx(ito+in*h, j, k) = rhs_in.jx(ifro+in*h, j, k);
      lhs_in.jy(ito+in*h, j, k) = rhs_in.jy(ifro+in*h, j, k);
      lhs_in.jz(ito+in*h, j, k) = rhs_in.jz(ifro+in*h, j, k);
    }, Ny, Nz, halo, lhs, rhs);
}


void copy_horz_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int Nx, int Nz, int halo, int jto, int jfro, int jn, int ind=0)
{
    UniIter::iterate3D([=] DEVCALLABLE (int i, int k ,int g, YeeLattice &lhs_in, YeeLattice &rhs_in){
      lhs_in.ex(i, jto+jn*g, k) = rhs_in.ex(i, jfro+jn*g, k);
      lhs_in.ey(i, jto+jn*g, k) = rhs_in.ey(i, jfro+jn*g, k);
      lhs_in.ez(i, jto+jn*g, k) = rhs_in.ez(i, jfro+jn*g, k);

      lhs_in.bx(i, jto+jn*g, k) = rhs_in.bx(i, jfro+jn*g, k);
      lhs_in.by(i, jto+jn*g, k) = rhs_in.by(i, jfro+jn*g, k);
      lhs_in.bz(i, jto+jn*g, k) = rhs_in.bz(i, jfro+jn*g, k);

      lhs_in.jx(i, jto+jn*g, k) = rhs_in.jx(i, jfro+jn*g, k);
      lhs_in.jy(i, jto+jn*g, k) = rhs_in.jy(i, jfro+jn*g, k);
      lhs_in.jz(i, jto+jn*g, k) = rhs_in.jz(i, jfro+jn*g, k);
    }, Nx, Nz, halo, lhs, rhs);
}


void copy_face_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int Nx, int Ny, int halo, int kto, int kfro, int kn, int ind=0)
{
    UniIter::iterate3D([=] DEVCALLABLE (int i, int j ,int g, YeeLattice &lhs_in, YeeLattice &rhs_in){ 
        lhs_in.ex(i, j, kto+kn*g) = rhs_in.ex(i, j, kfro+kn*g);
        lhs_in.ey(i, j, kto+kn*g) = rhs_in.ey(i, j, kfro+kn*g);
        lhs_in.ez(i, j, kto+kn*g) = rhs_in.ez(i, j, kfro+kn*g);

        lhs_in.bx(i, j, kto+kn*g) = rhs_in.bx(i, j, kfro+kn*g);
        lhs_in.by(i, j, kto+kn*g) = rhs_in.by(i, j, kfro+kn*g);
        lhs_in.bz(i, j, kto+kn*g) = rhs_in.bz(i, j, kfro+kn*g);

        lhs_in.jx(i, j, kto+kn*g) = rhs_in.jx(i, j, kfro+kn*g);
        lhs_in.jy(i, j, kto+kn*g) = rhs_in.jy(i, j, kfro+kn*g);
        lhs_in.jz(i, j, kto+kn*g) = rhs_in.jz(i, j, kfro+kn*g);
    }, Nx, Ny, halo, lhs, rhs);
}



void copy_x_pencil_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int Nx, int halo, int jto, int jfro, int kto, int kfro, int jn, int kn, int ind=0)
{   
    UniIter::iterate3D([=] DEVCALLABLE (int i, int g ,int h, YeeLattice &lhs_in, YeeLattice &rhs_in){
      lhs_in.ex(i, jto+jn*h, kto+kn*g) = rhs_in.ex(i, jfro+jn*h, kfro+kn*g);
      lhs_in.ey(i, jto+jn*h, kto+kn*g) = rhs_in.ey(i, jfro+jn*h, kfro+kn*g);
      lhs_in.ez(i, jto+jn*h, kto+kn*g) = rhs_in.ez(i, jfro+jn*h, kfro+kn*g);

      lhs_in.bx(i, jto+jn*h, kto+kn*g) = rhs_in.bx(i, jfro+jn*h, kfro+kn*g);
      lhs_in.by(i, jto+jn*h, kto+kn*g) = rhs_in.by(i, jfro+jn*h, kfro+kn*g);
      lhs_in.bz(i, jto+jn*h, kto+kn*g) = rhs_in.bz(i, jfro+jn*h, kfro+kn*g);

      lhs_in.jx(i, jto+jn*h, kto+kn*g) = rhs_in.jx(i, jfro+jn*h, kfro+kn*g);
      lhs_in.jy(i, jto+jn*h, kto+kn*g) = rhs_in.jy(i, jfro+jn*h, kfro+kn*g);
      lhs_in.jz(i, jto+jn*h, kto+kn*g) = rhs_in.jz(i, jfro+jn*h, kfro+kn*g);
    }, Nx, halo, halo, lhs, rhs);
  }

  void copy_y_pencil_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int Ny, int halo, int ito, int ifro, int kto, int kfro, int in, int kn, int ind=0)
{   
    UniIter::iterate3D([=] DEVCALLABLE (int j, int g ,int h, YeeLattice &lhs_in, YeeLattice &rhs_in){

      lhs_in.ex(ito+in*h, j, kto+kn*g) = rhs_in.ex(ifro+in*h, j, kfro+kn*g);
      lhs_in.ey(ito+in*h, j, kto+kn*g) = rhs_in.ey(ifro+in*h, j, kfro+kn*g);
      lhs_in.ez(ito+in*h, j, kto+kn*g) = rhs_in.ez(ifro+in*h, j, kfro+kn*g);

      lhs_in.bx(ito+in*h, j, kto+kn*g) = rhs_in.bx(ifro+in*h, j, kfro+kn*g);
      lhs_in.by(ito+in*h, j, kto+kn*g) = rhs_in.by(ifro+in*h, j, kfro+kn*g);
      lhs_in.bz(ito+in*h, j, kto+kn*g) = rhs_in.bz(ifro+in*h, j, kfro+kn*g);

      lhs_in.jx(ito+in*h, j, kto+kn*g) = rhs_in.jx(ifro+in*h, j, kfro+kn*g);
      lhs_in.jy(ito+in*h, j, kto+kn*g) = rhs_in.jy(ifro+in*h, j, kfro+kn*g);
      lhs_in.jz(ito+in*h, j, kto+kn*g) = rhs_in.jz(ifro+in*h, j, kfro+kn*g);
    }, Ny, halo, halo, lhs, rhs);
  }

  void copy_z_pencil_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int Nz, int halo, int ito, int ifro, int jto, int jfro, int in, int jn, int ind=0)
{
    UniIter::iterate3D([=] DEVCALLABLE (int k, int g ,int h, YeeLattice &lhs_in, YeeLattice &rhs_in){
      lhs_in.ex(ito+in*h, jto+jn*g, k) = rhs_in.ex(ifro+in*h, jfro+jn*g, k);
      lhs_in.ey(ito+in*h, jto+jn*g, k) = rhs_in.ey(ifro+in*h, jfro+jn*g, k);
      lhs_in.ez(ito+in*h, jto+jn*g, k) = rhs_in.ez(ifro+in*h, jfro+jn*g, k);

      lhs_in.bx(ito+in*h, jto+jn*g, k) = rhs_in.bx(ifro+in*h, jfro+jn*g, k);
      lhs_in.by(ito+in*h, jto+jn*g, k) = rhs_in.by(ifro+in*h, jfro+jn*g, k);
      lhs_in.bz(ito+in*h, jto+jn*g, k) = rhs_in.bz(ifro+in*h, jfro+jn*g, k);

      lhs_in.jx(ito+in*h, jto+jn*g, k) = rhs_in.jx(ifro+in*h, jfro+jn*g, k);
      lhs_in.jy(ito+in*h, jto+jn*g, k) = rhs_in.jy(ifro+in*h, jfro+jn*g, k);
      lhs_in.jz(ito+in*h, jto+jn*g, k) = rhs_in.jz(ifro+in*h, jfro+jn*g, k);
  
    }, Nz, halo, halo, lhs, rhs);
  }

  void copy_point_yee_halo(YeeLattice& lhs, YeeLattice& rhs, int halo, int ito, int ifro, int jto, int jfro, int kto, int kfro, int in, int jn, int kn, int ind=0)
  {
    UniIter::iterate3D([=] DEVCALLABLE (int f, int g ,int h, YeeLattice &lhs_in, YeeLattice &rhs_in){
      lhs_in.ex(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ex(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in.ey(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ey(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in.ez(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.ez(ifro+in*h, jfro+jn*g, kfro+kn*f);

      lhs_in.bx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.bx(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in.by(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.by(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in.bz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.bz(ifro+in*h, jfro+jn*g, kfro+kn*f);

      lhs_in.jx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jx(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in.jy(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jy(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in.jz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in.jz(ifro+in*h, jfro+jn*g, kfro+kn*f);
    }, halo, halo, halo, lhs, rhs);
  }


/// Update Yee grid boundaries
template<>
void Tile<1>::update_boundaries(corgi::Grid<1>& grid) 
{
  // target
  YeeLattice& mesh = get_yee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<Tile<1> >(
        grid.get_tileptr( neighs(-1) ));
  YeeLattice& mleft = cleft->get_yee();

  // copy from right side to left
  copy_vert_yee(mesh, mleft, -1, mleft.Nx-1); 

  // right
  auto cright = 
    std::dynamic_pointer_cast<Tile<1> >(
        grid.get_tileptr( neighs(+1) ));
  YeeLattice& mright = cright->get_yee();
    
  // copy from left side to right
  copy_vert_yee(mesh, mright, mesh.Nx, 0); 

}


template<>
void Tile<2>::update_boundaries(corgi::Grid<2>& grid) 
{
  using Tile_t  = Tile<2>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito, jto, ifro, jfro;
  Tileptr tpr;

  auto& mesh = get_yee(); // target as a reference to update into

  int halo = 3; // halo region size for fields

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      if (in == 0 && jn == 0) continue;

      tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn) ));
      if (tpr) {
        auto& mpr = tpr->get_yee();

        /* diagonal rules are:
        if + then to   n
        if + then from 0

        if - then to   -1
        if - then from n-1
        */

        if (in == +1) { ito = mesh.Nx; ifro = 0; }
        if (jn == +1) { jto = mesh.Ny; jfro = 0; }

        if (in == -1) { ito = -1;      ifro = mpr.Nx-1; }
        if (jn == -1) { jto = -1;      jfro = mpr.Ny-1; }

        // copy (halo = 1) assignment
        //if      (jn == 0) copy_vert_yee(    mesh, mpr, ito, ifro);   // vertical
        //else if (in == 0) copy_horz_yee(    mesh, mpr, jto, jfro);   // horizontal
        //else              copy_z_pencil_yee(mesh, mpr, ito, jto, ifro, jfro); // diagonal
        
        // generalized halo >= 1 loops
        if (jn == 0) { // vertical
          for(int h=0; h<halo; h++)
            copy_vert_yee(mesh, mpr, ito+in*h, ifro+in*h);   

        } else if (in == 0) { // horizontal
          for(int g=0; g<halo; g++)
            copy_horz_yee(mesh, mpr, jto+jn*g, jfro+jn*g);   

        } else { // diagonal
          for(int h=0; h<halo; h++) {
            for(int g=0; g<halo; g++) {
              copy_z_pencil_yee(mesh, mpr, ito+in*h, jto+jn*g, ifro+in*h, jfro+jn*g); 
            }
          }
        }

      } // end of if(tpr)
    }
  }
}

struct BlockCopyParam
{
  //
  bool active;
  int ito, ifro, jto, jfro, kto, kfro, in, jn, kn;
  YeeLattice *yeePtrFrom;
  YeeLattice *yeePtrTo;
};

template<>
void Tile<3>::update_boundaries(corgi::Grid<3>& grid) 
{
  //std::cout << "upB: updating boundaries\n";
  //nvtxRangePush(__FUNCTION__);

  using Tile_t  = Tile<3>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito, jto, kto, ifro, jfro, kfro;
  Tileptr tpr;

  auto& mesh = get_yee(); // target as a reference to update into

  int halo = 3; // halo region size for fields

  int ind = 0;

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      for(int kn=-1; kn <= 1; kn++) {

        if (in == 0 && jn == 0 && kn == 0) continue;

        // continue only if the tile exists (get_tileptr return nullptr otherwise)
        tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn, kn) ));
        if (tpr) {
          auto& mpr = tpr->get_yee();

          /* diagonal rules are:
          if + then to   n
          if + then from 0

          if - then to   -1
          if - then from n-1
          */

          if (in == +1) { ito = mesh.Nx; ifro = 0; }
          if (jn == +1) { jto = mesh.Ny; jfro = 0; }
          if (kn == +1) { kto = mesh.Nz; kfro = 0; }

          if (in == -1) { ito = -1;      ifro = mpr.Nx-1; }
          if (jn == -1) { jto = -1;      jfro = mpr.Ny-1; }
          if (kn == -1) { kto = -1;      kfro = mpr.Nz-1; }

          // generalized halo >= 1 loops

          if (kn == 0) {

            // vertical
            if (jn == 0) { 
            // horizontal
              copy_vert_yee_halo(mesh, mpr, mesh.ex.Ny, mesh.ex.Nz, halo, ito, ifro, in, ind);
            } else if (in == 0) { 
            // diagonal
              copy_horz_yee_halo(mesh, mpr, mesh.ex.Nx, mesh.ex.Nz, halo, jto, jfro, jn, ind);
            } else { 
             copy_z_pencil_yee_halo(mesh, mpr, mesh.ex.Nz, halo, ito, ifro, jto, jfro, in, jn, ind);
            } 
         
          // 3D case with kn != 0
          } else {
            
            // infront/behind directions
            if (in == 0 && jn == 0) { 
              copy_face_yee_halo(mesh, mpr, mesh.ex.Nx, mesh.ex.Ny, halo, kto, kfro, kn, ind);
            // 3D generalized diagonal locations
            // If the finite-difference scheme is purely non-diagonal
            // then these can be dropped off.
              
            // vertical wedges
            } else if (jn == 0) {
              // y pencils
              copy_y_pencil_yee_halo(mesh, mpr, mesh.ex.Ny, halo, ito, ifro, kto, kfro, in, kn, ind);
            // horizontal wedges
            } else if (in == 0) {
              // x pencils
              copy_x_pencil_yee_halo(mesh, mpr, mesh.ex.Nx, halo, jto, jfro, kto, kfro, jn, kn, ind);

            // corners
            } else {
              // pointwise
              copy_point_yee_halo(mesh, mpr, halo, ito, ifro, jto, jfro, kto, kfro, in, jn, kn, ind);

            } 
            ind++;
          } // 3D cases with kn != 0
        } // if tpr
      } // kn
    } // jn
  } // in
  
  UniIter::sync();

  //nvtxRangePop();
  }


template<>
void Tile<1>::exchange_currents(corgi::Grid<1>& grid) 
{

  // target
  YeeLattice& mesh = get_yee();

  int halo = 1; // halo region size for currents


  // left 
  auto cleft = 
    std::dynamic_pointer_cast<Tile<1> >(
        grid.get_tileptr( neighs(-1) ));
  YeeLattice& mleft = cleft->get_yee();

  // add from right side to left
  for(int h=1; h<= halo; h++)
  add_vert_yee(mesh, mleft, -h, mleft.Nx-h); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<Tile<1> >(
        grid.get_tileptr( neighs(+1) ));
  YeeLattice& mright = cright->get_yee();
    
  // add from left side to right
  for(int h=1; h<= halo; h++)
  add_vert_yee(mesh, mright, mesh.Nx+h-1, h-1); 

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

  int ito, jto, ifro, jfro;
  Tileptr tpr; 

  int halo = 3;

  auto& mesh = get_yee(); // target as a reference to update into

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      if (in == 0 && jn == 0) continue;

      tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn) ));
      if (tpr) {
        auto& mpr = tpr->get_yee();

        /* diagonal rules are:
        if + then to   n
        if + then from 0

        if - then to   -1
        if - then from n-1
        */

        if (in == +1) { ito = mesh.Nx; ifro = 0; }
        if (jn == +1) { jto = mesh.Ny; jfro = 0; }

        if (in == -1) { ito = -1;      ifro = mpr.Nx-1; }
        if (jn == -1) { jto = -1;      jfro = mpr.Ny-1; }

        // add
        if (jn == 0) { // vertical
          for(int h=1; h<=halo; h++)
            add_vert_yee(mesh, mpr, ito-in*h, ifro-in*h);   

        } else if (in == 0) { // horizontal
          for(int g=1; g<=halo; g++)
            add_horz_yee(mesh, mpr, jto-jn*g, jfro-jn*g);   

        } else { // diagonal
          for(int h=1; h<=halo; h++) {
            for(int g=1; g<=halo; g++) {
              add_z_pencil_yee(mesh, mpr, ito-in*h, jto-jn*g, ifro-in*h, jfro-jn*g); 
            }
          }
        }
      } // end of if(tpr)
    }
  }
}

template<>
void Tile<3>::exchange_currents(corgi::Grid<3>& grid) 
{
  //nvtxRangePush(__FUNCTION__);

  using Tile_t  = Tile<3>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito, jto, kto, ifro, jfro, kfro;
  Tileptr tpr; 

  int halo = 3;

  auto& mesh = get_yee(); // target as a reference to update into

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      for(int kn=-1; kn <= 1; kn++) {
        if (in == 0 && jn == 0) continue;

        if (in == 0 && jn == 0 && kn == 0) continue;

        tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn, kn) ));
        if (tpr) {
          auto& mpr = tpr->get_yee();

          /* diagonal rules are:
          if + then to   n
          if + then from 0

          if - then to   -1
          if - then from n-1
          */

          if (in == +1) { ito = mesh.Nx; ifro = 0; }
          if (jn == +1) { jto = mesh.Ny; jfro = 0; }
          if (kn == +1) { kto = mesh.Nz; kfro = 0; }

          if (in == -1) { ito = -1;      ifro = mpr.Nx-1; }
          if (jn == -1) { jto = -1;      jfro = mpr.Ny-1; }
          if (kn == -1) { kto = -1;      kfro = mpr.Nz-1; }


          // generalized halo >= 1 loops

          // 2D case; kn == 0; no tiles behind or infront
          if (kn == 0) {

            // vertical
            if (jn == 0) { 
              for(int h=1; h<=halo; h++)
                add_vert_yee(mesh, mpr, ito-in*h, ifro-in*h);   

            // horizontal
            } else if (in == 0) { 
              for(int g=1; g<=halo; g++)
                add_horz_yee(mesh, mpr, jto-jn*g, jfro-jn*g);   

            // diagonal
            } else { 
              for(int h=1; h<=halo; h++) {
                for(int g=1; g<=halo; g++) {
                  add_z_pencil_yee(mesh, mpr, ito-in*h, jto-jn*g, ifro-in*h, jfro-jn*g); 
                }
              }
            } 
         
          // 3D case with kn != 0
          } else {
            
            // infront/behind directions
            if (in == 0 && jn == 0) { 
              for(int g=1; g<=halo; g++)
                add_face_yee(mesh, mpr, kto-kn*g, kfro-kn*g);   

            // 3D generalized diagonal locations
            // If the finite-difference scheme is purely non-diagonal
            // then these can be dropped off.
              
            // vertical wedges
            } else if (jn == 0) {

              // y pencils
              for(int h=1; h<=halo; h++) {
                for(int g=1; g<=halo; g++) {
                  add_y_pencil_yee(mesh, mpr, ito-in*h, kto-kn*g, ifro-in*h, kfro-kn*g); 
                }
              }

            // horizontal wedges
            } else if (in == 0) {

              // x pencils
              for(int h=1; h<=halo; h++) {
                for(int g=1; g<=halo; g++) {
                  add_x_pencil_yee(mesh, mpr, jto-jn*h, kto-kn*g, jfro-jn*h, kfro-kn*g); 
                }
              }

            // corners
            } else {

              // pointwise
              for(int h=1; h<=halo; h++) {
                for(int g=1; g<=halo; g++) {
                  for(int f=1; f<=halo; f++) {
                    add_point_yee(mesh, mpr, 
                        ito -in*h, jto -jn*g, kto -kn*f,
                        ifro-in*h, jfro-jn*g, kfro-kn*f);
                  }
                }
              }
            } 

          } // 3D cases with kn != 0
        } // if tpr
      }
    }
  }
  //nvtxRangePop();
}



template<std::size_t D>
void Tile<D>::cycle_yee() 
{
  //yee.cycle();
  // do nothing since Yee's are not in a container atm
}



template<std::size_t D>
void Tile<D>::clear_current() 
{
  auto& yee = this->get_yee();
  yee.jx.clear();
  yee.jy.clear();
  yee.jz.clear();
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
  //nvtxRangePush(__FUNCTION__);

  auto& yee = get_yee(); 
  //std::cout << "SEND field to " << dest 
  //  << "nx " << yee.jx.size()
  //  << "ny " << yee.jy.size()
  //  << "nz " << yee.jz.size()
  //  << "\n";
  std::vector<mpi::request> reqs;
  UniIter::sync();
  if (mode == 0) {
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 0), yee.jx.data(), yee.jx.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 1), yee.jy.data(), yee.jy.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 2), yee.jz.data(), yee.jz.size()) );
  } else if (mode == 1) {
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 3), yee.ex.data(), yee.ex.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 4), yee.ey.data(), yee.ey.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 5), yee.ez.data(), yee.ez.size()) );
  } else if (mode == 2) {
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 6), yee.bx.data(), yee.bx.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 7), yee.by.data(), yee.by.size()) );
    reqs.emplace_back( comm.isend(dest, get_tag(tag, 8), yee.bz.data(), yee.bz.size()) );
  }

  //nvtxRangePop();

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_data( 
    mpi::communicator& comm, 
    int orig, 
    int mode,
    int tag)
{
  //nvtxRangePush(__FUNCTION__);

  //std::cout << "RECV from " << orig << "\n";
  auto& yee = get_yee(); 
  //std::cout << "RECV field to " << orig
  //  << "nx " << yee.jx.size()
  //  << "ny " << yee.jy.size()
  //  << "nz " << yee.jz.size()
  //  << "\n";

  std::vector<mpi::request> reqs;
  UniIter::sync();

  if (mode == 0) {
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 0), yee.jx.data(), yee.jx.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 1), yee.jy.data(), yee.jy.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 2), yee.jz.data(), yee.jz.size()) );
  } else if (mode == 1) {
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 3), yee.ex.data(), yee.ex.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 4), yee.ey.data(), yee.ey.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 5), yee.ez.data(), yee.ez.size()) );
  } else if (mode == 2) {
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 6), yee.bx.data(), yee.bx.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 7), yee.by.data(), yee.by.size()) );
    reqs.emplace_back( comm.irecv(orig, get_tag(tag, 8), yee.bz.data(), yee.bz.size()) );
  }

  //nvtxRangePop();

  return reqs;
}

//--------------------------------------------------
// explicit template instantiation

template class Tile<1>;
template class Tile<2>;
template class Tile<3>;

} // end of ns fields
