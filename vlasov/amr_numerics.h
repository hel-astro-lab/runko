#pragma once

#include <cmath> 
#include <algorithm>
#include <stdexcept>

#include "mesh.h"


namespace toolbox {



/*! \brief Finite difference partial derivative
 *
 * NOTE: Returns silent zero if derivative beyond mesh boundary is asked.
 */
template<typename T, int D>
T deriv(
    AdaptiveMesh<T,D>& mesh, 
    typename AdaptiveMesh<T,D>::indices_t& indices,
    int refinement_level,
    std::array<int, D>& directions
      ) 
{
  uint64_t cid = mesh.get_cell_from_indices(indices, refinement_level);

  for(auto dir : directions) if( std::abs(dir) > 1) throw std::invalid_argument("Direction is more than 1");


  // get relevant neighbor in my refinement level
  typename AdaptiveMesh<T,D>::indices_t indices_neighbor( indices );

  size_t i = 0;
  for(auto dir: directions) {

    // silently return zero if derivative beyond array range is asked
    if( (indices[i] <= 0) & (dir < 0))                 return T(0);
    if( (indices[i] >= mesh.length[i]-1) & (dir > 0) ) return T(0);

    indices_neighbor[i] += dir;
    i++;
  }
  uint64_t cid_neighbor = mesh.get_cell_from_indices(indices_neighbor, refinement_level);


  // TODO: implement usage of most refined neighbors; not current refinement level
  // then work on the most refined relevant neighbors 
  // std::vector<uint64_t> children = m.get_most_refined_children(cid_neighbor);
  //
  // instead we (for now) simplify and assume same refinement level
  // function values
  T 
    f0 = mesh.get(cid),
    f1 = mesh.get(cid_neighbor);


  // then compute distance between cells
  // NOTE: this operation can be done even if the cell at this ref. lvl. does not exist
  typename AdaptiveMesh<T,D>::value_array_t 
    x0 = mesh.get_center(indices,          refinement_level),
    x1 = mesh.get_center(indices_neighbor, refinement_level);

  T dx = T(0);
  for(size_t i=0; i<D; i++) dx += std::pow( x1[i] - x0[i], 2);
  dx = std::sqrt(dx);


  // finally the actual derivative
  T dfdx = (f1 - f0) / dx;


  return dfdx;
}



/*! \brief Gradient operator for amr mesh
 *
 * Uses deriv internally to get all the partial derivatives for each direction.
 * Centered derivative is obtained by computing average of +1 and -1 directions.
 */
template<typename T, int D>
typename AdaptiveMesh<T,D>::value_array_t grad(
    AdaptiveMesh<T,D>& mesh,
    typename AdaptiveMesh<T,D>::indices_t& indices,
    int refinement_level) 
{
  typename AdaptiveMesh<T,D>::value_array_t gradient;
  
  T fp, fm;
  std::array<int, D> dirsp, dirsm;

  for(size_t i=0; i<D; i++) {
    dirsp.fill( 0 ); 
    dirsp[i] = +1;
    fp = deriv<T, D>(mesh, indices, refinement_level, dirsp);

    dirsm.fill( 0 ); 
    dirsm[i] = -1;
    fm = deriv<T, D>(mesh, indices, refinement_level, dirsm);

    gradient[i] = 0.5*( fp + fm );
  }

  return gradient;
}















} // end of namespace toolbox
