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
    const AdaptiveMesh<T,D>& mesh, 
    const typename AdaptiveMesh<T,D>::indices_t& indices,
    int refinement_level,
    std::array<int, D>& directions
      ) 
{
  uint64_t cid = mesh.get_cell_from_indices(indices, refinement_level);

  for(auto dir : directions) if( std::abs(dir) > 1) throw std::invalid_argument("Direction is more than 1");


  // get relevant neighbor in my refinement level
  typename AdaptiveMesh<T,D>::indices_t 
    indices_neighbor( indices ),
    length;

  length = mesh.get_size(refinement_level);
  size_t i = 0;
  for(auto dir: directions) {

    // silently return zero if derivative beyond array range is asked
    if( (indices[i] <= 0) && (dir < 0))            return T(0);
    if( (indices[i] >= length[i]-1) && (dir > 0) ) return T(0);

    indices_neighbor[i] += dir;
    i++;
  }
  uint64_t cid_neighbor = mesh.get_cell_from_indices(indices_neighbor, refinement_level);


  // TODO: implement usage of most refined neighbors; not current refinement level
  // then work on the most refined relevant neighbors 
  // std::vector<uint64_t> children = m.get_most_refined_children(cid_neighbor); = leafs
  //
  // instead we (for now) simplify and assume same refinement level
  // function values; then look below for roots
  T 
    f0 = mesh.get_from_roots(cid),
    f1 = mesh.get_from_roots(cid_neighbor);


  // compute distance between cells
  // NOTE: this operation can be done even if the cell at this ref. lvl. does not exist
  typename AdaptiveMesh<T,D>::value_array_t 
    x0 = mesh.get_center(indices,          refinement_level),
    x1 = mesh.get_center(indices_neighbor, refinement_level);

  T dx = T(0);
  for(size_t i=0; i<D; i++) dx += std::pow( x1[i] - x0[i], 2);
  dx = std::sqrt(dx);


  // finally the actual derivative
  T dfdx = (f1 - f0) / dx;
  // T dfdx = (f1 - f0) / (f0 + 0.001);
  // T dfdx = (f1 - f0 + 1.0e-4);


  return dfdx;
}



/*! \brief Gradient operator for amr mesh
 *
 * Uses deriv internally to get all the partial derivatives for each direction.
 * Centered derivative is obtained by computing average of +1 and -1 directions.
 */
template<typename T, int D>
typename AdaptiveMesh<T,D>::value_array_t grad(
    const AdaptiveMesh<T,D>& mesh,
    const typename AdaptiveMesh<T,D>::indices_t& indices,
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


/// simple lerp function (= linear interpolation)
template<typename T>
static inline T lerp(T t, T v1, T v2) {return (T(1) - t)*v1 + t*v2;}


/// inlined value getter from the AMR mesh
template<typename T, int D>
static inline T _getv(
    const AdaptiveMesh<T,D>& mesh,
    const typename AdaptiveMesh<T,D>::indices_t& indices,
    int rfl)
{
  const uint64_t cid = mesh.get_cell_from_indices(indices, rfl);
  return mesh.get_from_roots(cid);
}



/// \brief Trilinear interpolation
// 
// Transforms from cell-centered values to vertex-centered,
// then uses the normal lerp routine.
//
// Partially based on:
// https://svn.blender.org/svnroot/bf-blender/
// branches/volume25/source/blender/blenlib/intern/voxel.c  
template<typename T>
T trilinear_interp(
    const AdaptiveMesh<T,3>& mesh,
    typename AdaptiveMesh<T,3>::indices_t& indices,
    typename AdaptiveMesh<T,3>::value_array_t coordinates,
    int rfl)
{
  uint64_t 
    i = indices[0],
    j = indices[1],
    k = indices[2];
	
	T dx = coordinates[0];// - T(0.5); 
  T dy = coordinates[1];// - T(0.5); 
  T dz = coordinates[2];// - T(0.5);
	
	T d00 = lerp(dx, _getv(mesh, {{i, j,   k  }}, rfl), _getv(mesh, {{i+1, j,   k  }}, rfl) );
	T d10 = lerp(dx, _getv(mesh, {{i, j+1, k  }}, rfl),	_getv(mesh, {{i+1, j+1, k  }}, rfl) );
	T d01 = lerp(dx, _getv(mesh, {{i, j,   k+1}}, rfl),	_getv(mesh, {{i+1, j,   k+1}}, rfl) );
	T d11 = lerp(dx, _getv(mesh, {{i, j+1, k+1}}, rfl),	_getv(mesh, {{i+1, j+1, k+1}}, rfl) );
	T d0  = lerp(dy, d00, d10);
	T d1  = lerp(dy, d01, d11);
	T d   = lerp(dz, d0, d1);
	
	return d;
}












} // end of namespace toolbox
