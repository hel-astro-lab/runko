#include "piston.h"
#include "../../definitions.h"

/// Perform phase-space reflection along x-axis
template<typename T>
inline void reflect_x(
    toolbox::AdaptiveMesh<T,3>& m)
{

  // find x=0 index for every rfl
  std::vector<T> flip_index;
  flip_index.resize( m.maximum_refinement_level );

  typename toolbox::AdaptiveMesh<T,3>::indices_t index = {{0,0,0}};

  for(int rfl=0; rfl<=m.maximum_refinement_level; rfl++) {
    for(size_t i=0; i<m.get_size(rfl)[0]; i++) {
      index[0] = i;
      auto uvel  = m.get_center(index, rfl);
      if ( uvel[0] >= T(0) ) {
        flip_index[rfl] = i;
        break;
      }
    }
  }

  /// Remove positive (top) of the distribution
  //for(auto cid : m.get_cells(true) ) {
  //  auto index = m.get_indices(cid);
  //  int rfl    = m.get_refinement_level(cid);
  //  auto uvel  = m.get_center(index, rfl);

  //  if( uvel[0] > T(0) ) m.data.erase(cid);
  //}


  // mirror & shave off bottom
  for(auto cid : m.get_cells(true) ) {
    auto index = m.get_indices(cid);
    int rfl    = m.get_refinement_level(cid);
    auto uvel  = m.get_center(index, rfl);

    if( uvel[0] < T(0) ) {

      // determine new location where x > 0
      size_t ii    = index[0];
      size_t flipi = flip_index[rfl] - ii;
      index[0] = flip_index[rfl] + flipi;
      uint64_t cid_top = m.get_cell_from_indices(index, rfl);
    
      // put there and remove 
      m.data[cid_top] = m.data.at(cid);
      //m.data.erase(cid); // finally erase
    }
  }

}


/// Tile reflection member that applies reflect_x into the Vlasov fluid 
// inside tile
template<size_t D>
void vlv::piston::Tile<D>::reflect( corgi::Grid<D>& grid )
{

  auto& step0 = vlv::Tile<D>::steps.get(0);

  //std::cout<<"BC reflecting...\n";

  // loop over different particle species (zips current [0] and new [1] solutions)
  int ispc = 0; // ith particle species
  for(auto&& block0 : step0 ) {
    auto& block0_right  = vlv::Tile<D>::get_external_data(grid, ispc, +1);

    int Nx = int(block0.Nx),
        Ny = int(block0.Ny),
        Nz = int(block0.Nz);


    /// middle of the tile in x-dir; reflection axis
    int xmiddle = Nx/2;

    // clear
    for (int s=0; s<Nz; s++) {
      for(int r=0; r<Ny; r++) {
        for(int q=0; q<Nx; q++) {
          auto& N   = block0.block(q,r,s);   // f_i
          N.data.clear();
        }
      }
    }

    // reflect in respect to rightmost wall
    for (int s=0; s<Nz; s++) {
      for(int r=0; r<Ny; r++) {

        //int q = xmiddle+1;
        int q = 0;
        for(int qinv=Nx-1; qinv>xmiddle; qinv--) {
          //for(int qinv=Nx-1; qinv>=0; qinv--) {
          const auto& M   = block0_right.block(q,r,s);   // f_i+1
          //const auto& M   = block0.block(q,r,s);    // f_i+1
          auto& N         = block0.block(qinv,r,s); // f_i
          N = M;

          reflect_x(N);
          q++;
        }

        //remove rest
        //for(int q=0; q<=1; q++) {
        //  auto& M = block0.block(0,r,s);
        //  M.data.clear();
        //}

      }
    }
  }
}



//--------------------------------------------------
// explicit template specialization
template class vlv::piston::Tile<1>;
