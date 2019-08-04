#include <iostream>
#include <cmath>

#include "tile.h"


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
YeeLattice& Tile<D>::get_yee(size_t /*i*/) 
{
  //return this->yee.get(i);
  return this->yee;
}

template<std::size_t D>
const YeeLattice& Tile<D>::get_const_yee(size_t /*i*/) const 
{
  //return this->yee.get(i);
  return this->yee;
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
  //std::cout << "add 1D Yee \n";
  //yee.push_back( YeeLattice( mesh_lengths[0], 1, 1) );
  yee = YeeLattice(mesh_lengths[0], 1, 1);
}

template<>
void Tile<2>::add_yee_lattice() 
{
  //std::cout << "add 2D Yee \n";
  //yee.push_back( YeeLattice( mesh_lengths[0], mesh_lengths[1], 1) );
  yee = YeeLattice( mesh_lengths[0], mesh_lengths[1], 1);
}

template<>
void Tile<3>::add_yee_lattice() 
{
  //std::cout << "add 3D Yee \n";
  //yee.push_back( YeeLattice( mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]) );
  yee = YeeLattice( mesh_lengths[0], mesh_lengths[1], mesh_lengths[2]);
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

        // copy
        if      (jn == 0) copy_vert_yee(    mesh, mpr, ito, ifro);   // vertical
        else if (in == 0) copy_horz_yee(    mesh, mpr, jto, jfro);   // horizontal
        else              copy_z_pencil_yee(mesh, mpr, ito, jto, ifro, jfro); // diagonal
        

        /*
        FIXME: this is how we should loop over H>1 halo boundaries
        for(int h=1; h<= halo; h++)
        copy_vert_yee(mesh, mleft, -h, mleft.Nx-h); 

        for(int h=1; h<= halo; h++)
        copy_horz_yee(mesh, mtop, mesh.Ny+h-1, h-1); 

        for(int h=1; h<= halo; h++)
        for(int g=1; g<= halo; g++)
        copy_z_pencil_yee(mesh, mtopleft, -h, mesh.Ny +g-1, mtopleft.Nx-h, +g-1);

        for(int h=1; h<= halo; h++)
        for(int g=1; g<= halo; g++)
        copy_z_pencil_yee(mesh, mtopright, mesh.Nx +h-1, mesh.Ny +g-1, +h-1,+g-1);
        */


      } // end of if(tpr)
    }
  }
}


template<>
void Tile<3>::update_boundaries(corgi::Grid<3>& grid) 
{
  using Tile_t  = Tile<3>;
  using Tileptr = std::shared_ptr<Tile_t>;

  int ito, jto, kto, ifro, jfro, kfro;
  Tileptr tpr;

  auto& mesh = get_yee(); // target as a reference to update into

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      for(int kn=-1; kn <= 1; kn++) {
        if (in == 0 && jn == 0 && kn == 0) continue;

        tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn, kn) ));
        if (tpr) {
          auto& mpr = tpr->get_yee();

          //TODO: implement 3D
          assert(false);

        }
      }
    }
  }

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
// TODO: assumes implicitly 2D (x-y) arrays only by setting k=0 and then ignoring it
// TODO: write unit test for this
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

        tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn, kn) ));
        if (tpr) {
          auto& mpr = tpr->get_yee();

          //TODO: implement 3D
          assert(false);

        }
      }
    }
  }

}



template<std::size_t D>
void Tile<D>::cycle_yee() 
{
  //yee.cycle();
  // do nothing since Yee's are not in a container atm
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
  auto& yee = get_yee(); 
  //std::cout << "SEND field to " << dest 
  //  << "nx " << yee.jx.size()
  //  << "ny " << yee.jy.size()
  //  << "nz " << yee.jz.size()
  //  << "\n";
  std::vector<mpi::request> reqs;

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

  return reqs;
}


template<std::size_t D>
std::vector<mpi::request> Tile<D>::recv_data( 
    mpi::communicator& comm, 
    int orig, 
    int mode,
    int tag)
{
  //std::cout << "RECV from " << orig << "\n";
  auto& yee = get_yee(); 
  //std::cout << "RECV field to " << orig
  //  << "nx " << yee.jx.size()
  //  << "ny " << yee.jy.size()
  //  << "nz " << yee.jz.size()
  //  << "\n";

  std::vector<mpi::request> reqs;


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


  return reqs;
}


//--------------------------------------------------
// explicit template instantiation

template class Tile<1>;
template class Tile<2>;
template class Tile<3>;

} // end of ns fields
