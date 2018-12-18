#include <iostream>

#include "tile.h"


namespace fields {
  using namespace mpi4cpp;


/* 
 * 1D version:
 *	ey(i,j,k)=ey(i,j,k) + c *(bz(im1,j,k)-bz(i,j,k))  
 *	ez(i,j,k)=ez(i,j,k) + c *(by(i,j,k)-by(im1,j,k)) 
 *
 * 2D version
 * ex(i,j,k)=ex(i,j,k)+const*(-bz(i,jm1,k)+bz(i,j,k))
 * ey(i,j,k)=ey(i,j,k)+const*(bz(im1,j,k)-bz(i,j,k))
 * ez(i,j,k)=ez(i,j,k)+const*(bx(i,jm1,k)-bx(i,j,k)-  by(im1,j,k)+by(i,j,k))
*/

/*! \brief Update E field with full step
 *
 * Contains a dimension switch for solvers depending on internal mesh dimensions
 */


/// 1D E pusher
template<>
void Tile<1>::push_e() 
{
  YeeLattice& mesh = get_yee();
  Realf C = 1.0 * cfl;

  int k = 0;
  int j = 0;
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Ex
    // NONE

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k));

  }

}

/// 2D E pusher
template<>
void Tile<2>::push_e() 
{
  YeeLattice& mesh = get_yee();

  Realf C = 1.0 * cfl;

  int k = 0;
  for(int j=0; j<static_cast<int>(mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Ex
    mesh.ex(i,j,k) += 
      + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k) 
          -mesh.by(i-1,j,   k) + mesh.by(i,j,k));

  }
}


/// 3D E pusher
template<>
void Tile<3>::push_e() 
{
  YeeLattice& mesh = get_yee();
  Realf C = 1.0 * cfl;

  for(int k=0; k<static_cast<int>(mesh_lengths[2]); k++)
  for(int j=0; j<static_cast<int>(mesh_lengths[1]); j++)
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Ex
    mesh.ex(i,j,k) += 
      + C*( mesh.by(i,j,  k-1) - mesh.by(i,j,k))
      + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

    // Ey
    mesh.ey(i,j,k) += 
      + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k))
      + C*(-mesh.bx(i,  j, k-1) + mesh.bx(i,j,k));

    // Ez
    mesh.ez(i,j,k) += 
      + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k))
      + C*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k));

  }

}


//--------------------------------------------------


/// Deposit current into electric field
template<std::size_t D>
void Tile<D>::deposit_current() 
{
  YeeLattice& mesh = get_yee();

  mesh.ex -= mesh.jx;
  mesh.ey -= mesh.jy;
  mesh.ez -= mesh.jz;

}

//--------------------------------------------------

/*
 * 1D version:
  by(i,j,k)=by(i,j,k)+const*(ez(ip1,j,k)-ez(i,j,k)) !-0*ex(i,j,k+1)+0*ex(i,j,k))
	bz(i,j,k)=bz(i,j,k)+const*(ey(i,j,k)-ey(ip1,j,k)) !+0*ex(i,j+1,k)-0*ex(i,j,k))
  
 * 2D version:
 * bx(i,j,k)=bx(i,j,k)+const*(-ez(i,jp1,k)+ez(i,j,k))
 * by(i,j,k)=by(i,j,k)+const*(ez(ip1,j,k)-ez(i,j,k))
 * bz(i,j,k)=bz(i,j,k)+const*(ex(i,jp1,k)-ex(i,j,k) -ey(ip1,j,k)+ey(i,j,k))
*/

/// Update B field with a half step


/// 1D B pusher
template<>
void Tile<1>::push_half_b() 
{
  YeeLattice& mesh = get_yee();
  Realf C = 0.5 * cfl;

  int k = 0;
  int j = 0;
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Bx
    // NONE

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));
  }

}

/// 2D B pusher
template<>
void Tile<2>::push_half_b() 
{
  YeeLattice& mesh = get_yee();

  Realf C = 0.5 * cfl;

  int k = 0;
  for(int j=0; j<static_cast<int>(mesh_lengths[1]); j++)
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Bx
    mesh.bx(i,j,k) += 
      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k)
          -mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

  }

}


/// 3D B pusher
template<>
void Tile<3>::push_half_b() 
{
  YeeLattice& mesh = get_yee();
  Realf C = 0.5 * cfl;

  for(int k=0; k<static_cast<int>(mesh_lengths[2]); k++) 
  for(int j=0; j<static_cast<int>(mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(mesh_lengths[0]); i++) {

    // Bx
    mesh.bx(i,j,k) += 
      + C*( mesh.ey(i,  j,  k+1) - mesh.ey(i,j,k))
      + C*(-mesh.ez(i,  j+1,k  ) + mesh.ez(i,j,k));

    // By
    mesh.by(i,j,k) += 
      + C*( mesh.ez(i+1,j, k  ) - mesh.ez(i,j,k))
      + C*(-mesh.ex(i,  j, k+1) + mesh.ex(i,j,k));

    // Bz
    mesh.bz(i,j,k) += 
      + C*( mesh.ex(i,  j+1, k) - mesh.ex(i,j,k))
      + C*(-mesh.ey(i+1,j,   k) + mesh.ey(i,j,k));

  }
}




/// Get current time snapshot of Yee lattice
template<std::size_t D>
YeeLattice& Tile<D>::get_yee(size_t i) 
{
  return this->yee.get(i);
}

template<std::size_t D>
const YeeLattice& Tile<D>::get_const_yee(size_t i) const 
{
  return this->yee.get(i);
}


/// Get analysis lattice of i:th species
template<std::size_t D>
PlasmaMomentLattice& Tile<D>::get_analysis(size_t i) 
{
  return this->analysis.at(i);
}

template<std::size_t D>
const PlasmaMomentLattice& Tile<D>::get_const_analysis(size_t i) const 
{
  return this->analysis.at(i);
}

//--------------------------------------------------
// Specialize analysis species grid extension
template<>
void Tile<1>::add_analysis_species() 
{
  analysis.emplace_back(mesh_lengths[0], 1, 1);
}

template<>
void Tile<2>::add_analysis_species() 
{
  analysis.emplace_back(mesh_lengths[0], mesh_lengths[1], 1);
}

template<>
void Tile<3>::add_analysis_species() 
{
  analysis.emplace_back(mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]);
}


//--------------------------------------------------
// Specialize Yee Lattice insertion
template<>
void Tile<1>::add_yee_lattice() 
{
  yee.push_back( YeeLattice( mesh_lengths[0], 1, 1) );
}

template<>
void Tile<2>::add_yee_lattice() 
{
  yee.push_back( YeeLattice( mesh_lengths[0], mesh_lengths[1], 1) );
}

template<>
void Tile<3>::add_yee_lattice() 
{
  yee.push_back( YeeLattice( mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]) );
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




/// Update Yee grid boundaries
template<>
void Tile<1>::update_boundaries(corgi::Node<1>& node) 
{
  // target
  YeeLattice& mesh = get_yee();

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<Tile<1> >(
        node.get_tileptr( neighs(-1) ));
  YeeLattice& mleft = cleft->get_yee();

  // copy from right side to left
  copy_vert_yee(mesh, mleft, -1, mleft.Nx-1); 

  // right
  auto cright = 
    std::dynamic_pointer_cast<Tile<1> >(
        node.get_tileptr( neighs(+1) ));
  YeeLattice& mright = cright->get_yee();
    
  // copy from left side to right
  copy_vert_yee(mesh, mright, mesh.Nx, 0); 

}


/// Update Yee grid boundaries
template<>
void Tile<2>::update_boundaries(corgi::Node<2>& node) 
{

  // target
  YeeLattice& mesh = get_yee();
  int halo = 1; // halo region size for current

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(-1, 0) ));
  YeeLattice& mleft = cleft->get_yee();

  // copy from right side to left
  for(int h=1; h<= halo; h++)
  copy_vert_yee(mesh, mleft, -h, mleft.Nx-h); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(+1, 0) ));
  YeeLattice& mright = cright->get_yee();
    
  // copy from left side to right
  for(int h=1; h<= halo; h++)
  copy_vert_yee(mesh, mright, mesh.Nx+h-1, h-1); 


  // TODO: fix these: they produce saw-like oscillations
  // top 
  auto ctop = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(0, +1) ));
  YeeLattice& mtop = ctop->get_yee();

  //copy from bottom side to top
  for(int h=1; h<= halo; h++)
  copy_horz_yee(mesh, mtop, mesh.Ny+h-1, h-1); 


  // bottom
  auto cbot = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(0, -1) ));
  YeeLattice& mbot = cbot->get_yee();
    
  // copy from top side to bottom
  for(int h=1; h<= halo; h++)
  copy_horz_yee(mesh, mbot, -h, mbot.Ny-h); 


  // --------------------------------------------------  
  // diagonals
  // TODO: loop over H

  auto ctopleft = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(-1, +1) ));
  YeeLattice& mtopleft = ctopleft->get_yee();

  for(int h=1; h<= halo; h++)
  for(int g=1; g<= halo; g++)
  copy_z_pencil_yee(mesh, mtopleft, -h,           mesh.Ny +g-1,
                                   mtopleft.Nx-h, +g-1);

  auto ctopright = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(+1, +1) ));
  YeeLattice& mtopright = ctopright->get_yee();

  for(int h=1; h<= halo; h++)
  for(int g=1; g<= halo; g++)
  copy_z_pencil_yee(mesh, mtopright, mesh.Nx +h-1, mesh.Ny +g-1,
                                     +h-1,         +g-1);

  auto cbotleft = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(-1, -1) ));
  YeeLattice& mbotleft = cbotleft->get_yee();

  for(int h=1; h<= halo; h++)
  for(int g=1; g<= halo; g++)
  copy_z_pencil_yee(mesh, mbotleft, -h, -g,
                          mbotleft.Nx-h, mbotleft.Ny-g);

  auto cbotright = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(+1, -1) ));
  YeeLattice& mbotright = cbotright->get_yee();

  for(int h=1; h<= halo; h++)
  for(int g=1; g<= halo; g++)
  copy_z_pencil_yee(mesh, mbotright, mesh.Nx +h-1, -g,
                                     +h-1,         mbotright.Ny-g);


  // --------------------------------------------------  
  // front
  // TODO: hack to deal with 2D corgi tiles
  //copy_face_yee(mesh, mesh, -1, mesh.Nz-1);

  // back
  //copy_face_yee(mesh, mesh, mesh.Nz, 0);
  // --------------------------------------------------  
  // TODO: x-pencils
  // TODO: y-pencils
  // TODO: corners


}

template<>
void Tile<1>::exchange_currents(corgi::Node<1>& node) 
{

  // target
  YeeLattice& mesh = get_yee();

  int halo = 1; // halo region size for currents


  // left 
  auto cleft = 
    std::dynamic_pointer_cast<Tile<1> >(
        node.get_tileptr( neighs(-1) ));
  YeeLattice& mleft = cleft->get_yee();

  // add from right side to left
  for(int h=1; h<= halo; h++)
  add_vert_yee(mesh, mleft, -h, mleft.Nx-h); 


  // right
  auto cright = 
    std::dynamic_pointer_cast<Tile<1> >(
        node.get_tileptr( neighs(+1) ));
  YeeLattice& mright = cright->get_yee();
    
  // add from left side to right
  for(int h=1; h<= halo; h++)
  add_vert_yee(mesh, mright, mesh.Nx+h-1, h-1); 

}


/// Update currents on Yee grid boundaries
// TODO: assumes implicitly 2D (x-y) arrays only by setting k=0 and then ignoring it
// TODO: write unit test for this
template<>
void Tile<2>::exchange_currents(corgi::Node<2>& node) 
{
  // target
  YeeLattice& mesh = get_yee();

  int halo = 3;

  // left 
  auto cleft = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(-1, 0) ));
  YeeLattice& mleft = cleft->get_yee();

  // add from left to right
  for(int h=0; h< halo; h++) add_vert_yee(mesh, mleft, h, mleft.Nx+h); 

  // right
  auto cright = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(+1, 0) ));
  YeeLattice& mright = cright->get_yee();
    
  // add from right to left
  for(int h=1; h<= halo; h++) add_vert_yee(mesh, mright, mesh.Nx-h, -h); 


  // top 
  auto ctop = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(0, +1) ));
  YeeLattice& mtop = ctop->get_yee();

  //add from top to bottom
  for(int h=0; h< halo; h++) add_horz_yee(mesh, mtop, h, mtop.Ny+h); 


  // bottom
  auto cbot = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(0, -1) ));
  YeeLattice& mbot = cbot->get_yee();
    
  // add from bottom to top
  for(int h=1; h<=halo; h++) add_horz_yee(mesh, mbot, mesh.Ny-h, -h); 



  // --------------------------------------------------  
  // diagonals
  auto ctopleft = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(-1, +1) ));
  YeeLattice& mtopleft = ctopleft->get_yee();

  for(int h=0; h<  halo; h++)
  for(int g=1; g<= halo; g++)
  add_z_pencil_yee(mesh, mtopleft, h, mesh.Ny-g,
                                   mtopleft.Nx+h, -g);


  auto ctopright = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(+1, +1) ));
  YeeLattice& mtopright = ctopright->get_yee();

  for(int h=1; h<= halo; h++)
  for(int g=1; g<= halo; g++)
  add_z_pencil_yee(mesh, mtopright, mesh.Nx-h, mesh.Ny-g,
                                    -h,         -g);

  auto cbotleft = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(-1, -1) ));
  YeeLattice& mbotleft = cbotleft->get_yee();
  
  for(int h=0; h<  halo; h++)
  for(int g=0; g<  halo; g++)
  add_z_pencil_yee(mesh, mbotleft, h,            g,
                                   mbotleft.Nx+h, mbotleft.Ny+g);

  auto cbotright = 
    std::dynamic_pointer_cast<Tile<2> >(
        node.get_tileptr( neighs(+1, -1) ));
  YeeLattice& mbotright = cbotright->get_yee();

  for(int h=1; h<= halo; h++)
  for(int g=0; g<  halo; g++)
  add_z_pencil_yee(mesh, mbotright, mesh.Nx-h, g,
                                    -h,      mbotright.Ny+g);


  // front
  // TODO: hack to deal with 2D corgi tiles
  //add_face_yee(mesh, mesh, -1, mesh.Nz-1);

  // back
  //add_face_yee(mesh, mesh, mesh.Nz, 0);

}

template<std::size_t D>
void Tile<D>::cycle_yee() 
{
  yee.cycle();
}

/// cycle temporary and true current arrays
template<std::size_t D>
void Tile<D>::cycle_current() 
{
  auto& yee = this->get_yee();

  std::swap( yee.jx.mat, yee.jx1.mat );
  std::swap( yee.jy.mat, yee.jy1.mat );
  std::swap( yee.jz.mat, yee.jz1.mat );

}



//--------------------------------------------------
// MPI routines

// create MPI tag given tile id and extra layer of differentiation
int get_tag(int cid, int extra_param)
{
  assert(extra_param < 100);
  return cid + extra_param*1e6;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::send_data( mpi::communicator& comm, int dest, int /*tag*/)
{
  auto& yee = get_yee(); 
  std::cout << "SEND field to " << dest 
    << "nx " << yee.jx.size()
    << "ny " << yee.jy.size()
    << "nz " << yee.jz.size()
    << "\n";

  std::vector<mpi::request> reqs;
  reqs.push_back( comm.isend(dest, get_tag(corgi::Tile<D>::cid, 1), yee.jx.data(), yee.jx.size()) );
  reqs.push_back( comm.isend(dest, get_tag(corgi::Tile<D>::cid, 2), yee.jy.data(), yee.jy.size()) );
  reqs.push_back( comm.isend(dest, get_tag(corgi::Tile<D>::cid, 3), yee.jz.data(), yee.jz.size()) );

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_data( mpi::communicator& comm, int orig, int /*tag*/)
{
  //std::cout << "RECV from " << orig << "\n";
  auto& yee = get_yee(); 
  std::cout << "RECV field to " << orig
    << "nx " << yee.jx.size()
    << "ny " << yee.jy.size()
    << "nz " << yee.jz.size()
    << "\n";

  std::vector<mpi::request> reqs;
  reqs.push_back( comm.irecv(orig, get_tag(corgi::Tile<D>::cid, 1), yee.jx.data(), yee.jx.size()) );
  reqs.push_back( comm.irecv(orig, get_tag(corgi::Tile<D>::cid, 2), yee.jy.data(), yee.jy.size()) );
  reqs.push_back( comm.irecv(orig, get_tag(corgi::Tile<D>::cid, 3), yee.jz.data(), yee.jz.size()) );

  return reqs;
}







//--------------------------------------------------
// explicit template instantiation

template class Tile<1>;
template class Tile<2>;
//template class Tile<3>;

} // end of ns fields
