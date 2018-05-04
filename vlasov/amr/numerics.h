#pragma once

#include <cmath> 
#include <algorithm>
#include <stdexcept>

#include "../../definitions.h"
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


/// simple lerp function (= linear interpolation) for corner/vertex-centered values
template<typename T>
static inline T lerp_corner(T t, T v1, T v2) {return (T(1) - t)*v1 + t*v2;}

/// simple lerp function (= linear interpolation) for cell-centered values
template<typename T>
static inline T lerp_center(T t, T v0, T v1, T v2) {return T(0.5)*(v1 + (T(1) - t)*v0 + t*v2);}


/// inlined value getter from the AMR mesh
template<typename T>
static inline T _getv(
    const AdaptiveMesh<T,3>& mesh,
    const typename AdaptiveMesh<T,3>::indices_t& indices,
    int rfl)
{
  const uint64_t cid = mesh.get_cell_from_indices(indices, rfl);
  return mesh.get_from_roots(cid);
}


/// inlined value getter from the AMR mesh with specialized weighting
// TODO: does not work
//--------------------------------------------------
// template<typename T, int D>
// static inline T _getvD(
//     const AdaptiveMesh<T,D>& mesh,
//     const typename AdaptiveMesh<T,D>::indices_t& indices,
//     int rfl)
// {
//   const uint64_t cid = mesh.get_cell_from_indices(indices, rfl);
//   std::pair<T,int> val = mesh.get_value_and_level(cid);
// 
//   return val.first;
//   //return val.second == 0 ? val.first : val.first * (2.0/3.0)*val.second;
// }
// 
// template<typename T, int D>
// static inline T _getvU(
//     const AdaptiveMesh<T,D>& mesh,
//     const typename AdaptiveMesh<T,D>::indices_t& indices,
//     int rfl)
// {
//   const uint64_t cid = mesh.get_cell_from_indices(indices, rfl);
//   std::pair<T,int> val = mesh.get_value_and_level(cid);
// 
//   return val.first;
//   //return val.second == 0 ? val.first : val.first * (1.0/3.0)*val.second;
// }
//--------------------------------------------------
  

template<typename T, int D>
inline T interp_linear(
    const AdaptiveMesh<T,3>& mesh,
    typename AdaptiveMesh<T,3>::indices_t& indices,
    typename AdaptiveMesh<T,3>::value_array_t coordinates,
    int rfl) = delete;


/// \brief Linear interpolation (specialization to D=1)
//
// XXX WARNING: not tested
//
template<>
inline Realf interp_linear<Realf,1>(
    const AdaptiveMesh<Realf,3>& mesh,
    typename AdaptiveMesh<Realf,3>::indices_t& indices,
    typename AdaptiveMesh<Realf,3>::value_array_t coordinates,
    int rfl)
{
  uint64_t 
    i = indices[0], // + uint64_t(coordinates[0] - 1.5),
    j = indices[1], // + uint64_t(coordinates[1] - 1.5),
    k = indices[2]; // + uint64_t(coordinates[2] - 1.5);
	
	Realf dx = coordinates[0]; // - T(0.5); 

  Realf d00 = lerp_corner(dx, _getv(mesh, {{i, j,   k  }}, rfl), _getv(mesh, {{i+1, j,   k  }}, rfl) );


  // Mesh is actually corner-stored so does not work:
	//Realf d00 = lerp_center(dx, 
  //    _getv(mesh, {{i-1, j,   k  }}, rfl), 
  //    _getv(mesh, {{i,   j,   k  }}, rfl), 
  //    _getv(mesh, {{i+1, j,   k  }}, rfl) 
  //    );
	
	return d00;
}

/// \brief Trilinear interpolation (specialization to D=3)
// 
// Transforms from cell-centered values to vertex-centered,
// then uses the normal lerp routine.
// XXX WARNING: no it does not
//
// Partially based on:
// https://svn.blender.org/svnroot/bf-blender/
// branches/volume25/source/blender/blenlib/intern/voxel.c  
template<>
inline Realf interp_linear<Realf,3>(
    const AdaptiveMesh<Realf,3>& mesh,
    typename AdaptiveMesh<Realf,3>::indices_t& indices,
    typename AdaptiveMesh<Realf,3>::value_array_t coordinates,
    int rfl)
{
  uint64_t 
    i = indices[0], // + uint64_t(coordinates[0] - 1.5),
    j = indices[1], // + uint64_t(coordinates[1] - 1.5),
    k = indices[2]; // + uint64_t(coordinates[2] - 1.5);
	
	Realf dx = coordinates[0]; // - T(0.5); 
  Realf dy = coordinates[1]; // - T(0.5); 
  Realf dz = coordinates[2]; // - T(0.5);
	
	Realf d00 = lerp_corner(dx, _getv(mesh, {{i, j,   k  }}, rfl), _getv(mesh, {{i+1, j,   k  }}, rfl) );
	Realf d10 = lerp_corner(dx, _getv(mesh, {{i, j+1, k  }}, rfl), _getv(mesh, {{i+1, j+1, k  }}, rfl) );
	Realf d01 = lerp_corner(dx, _getv(mesh, {{i, j,   k+1}}, rfl), _getv(mesh, {{i+1, j,   k+1}}, rfl) );
	Realf d11 = lerp_corner(dx, _getv(mesh, {{i, j+1, k+1}}, rfl), _getv(mesh, {{i+1, j+1, k+1}}, rfl) );
	Realf d0  = lerp_corner(dy, d00, d10);
	Realf d1  = lerp_corner(dy, d01, d11);
	Realf d   = lerp_corner(dz, d0,  d1);
	
	return d;
}



/// \brief Tricubic interpolation
//
// Local B-spline implementation that estimates derivatives (and second derivs)
// from neighboring grid points using finite difference.
//
// Partially based on:
// https://svn.blender.org/svnroot/bf-blender/
// branches/volume25/source/blender/blenlib/intern/voxel.c  
//
// and on:
//
// https://github.com/igmhub/likely/blob/master/likely/TriCubicInterpolator.cc
  
  
const static int _C[64][64] = {
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
    {-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0},
    { 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
    {-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1},
    {18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1},
    {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
    {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
    {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},
    { 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
    {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
    {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
    {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},
    { 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
    {-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
    { 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}
};
  
/*
 * cube getter; currently inlined due to optimization
template<typename T>
inline void _get_cube(
    const AdaptiveMesh<T,3>& mesh,
    const int i,
    const int j,
    const int k,
    const int rfl,
    int& x) 
{
  x = { }
}
*/

template<typename T, int D>
inline T interp_cubic(
    const AdaptiveMesh<T,3>& mesh,
    typename AdaptiveMesh<T,3>::indices_t& indices,
    typename AdaptiveMesh<T,3>::value_array_t& coordinates,
    int rfl) = delete;


/// Cubic interpolation (Catmull-Rom) (specialization to D=1)
// https://en.wikipedia.org/wiki/Cubic_Hermite_spline
//   Interpolation on the unit interval
//
//   and https://arxiv.org/pdf/0905.3564.pdf
//
template<>
inline Realf interp_cubic<Realf,1>(
    const AdaptiveMesh<Realf,3>& mesh,
    typename AdaptiveMesh<Realf,3>::indices_t& indices,
    typename AdaptiveMesh<Realf,3>::value_array_t& coordinates,
    int rfl)
{
  uint64_t 
    i = indices[0], 
    j = indices[1], 
    k = indices[2]; 
	
	Realf dx = coordinates[0]; // - T(0.5); 
  //Realf dy = coordinates[1]; // - T(0.5); 
  //Realf dz = coordinates[2]; // - T(0.5);

	Realf pm1 = _getv(mesh, {{i-1, j  , k}}  , rfl),
	      p0  = _getv(mesh, {{i  , j  , k}}  , rfl),
	      pp1 = _getv(mesh, {{i+1, j  , k}}  , rfl),
	      pp2 = _getv(mesh, {{i+2, j  , k}}  , rfl);

  return 0.5*(
      (dx*dx*(2.0-dx)-dx)*     pm1 +
      (dx*dx*(3.0*dx-5.0)+2.0)*p0 +
      (dx*dx*(4.0-3.0*dx)+dx)* pp1 +
       dx*dx*(dx-1.0)*         pp2);
}

/*
 * tricubic interpolator with gradient; 
 * TODO Needs to be optimized
template<typename T>
T tricubic_interp(
    const AdaptiveMesh<T,3>& mesh,
    typename AdaptiveMesh<T,3>::indices_t& indices,
    typename AdaptiveMesh<T,3>::value_array_t& coordinates,
    typename AdaptiveMesh<T,3>::value_array_t& grad,
    int rfl)
*/


/// Tricubic interpolator (specialization to D=3)
template<>
inline Realf interp_cubic<Realf,3>(
    const AdaptiveMesh<Realf,3>& mesh,
    typename AdaptiveMesh<Realf,3>::indices_t& indices,
    typename AdaptiveMesh<Realf,3>::value_array_t& coordinates,
    int rfl)
{
  uint64_t 
    i = indices[0], // - uint64_t(coordinates[0] - 1.5),
    j = indices[1], // - uint64_t(coordinates[1] - 1.5),
    k = indices[2]; // - uint64_t(coordinates[2] - 1.5);
	
	Realf dx = coordinates[0]; // - T(0.5); 
  Realf dy = coordinates[1]; // - T(0.5); 
  Realf dz = coordinates[2]; // - T(0.5);

  /*
  T x[64];
  _get_cube(mesh, i, j, k, rfl, x);
  */

  // Extract the local vocal values and calculate partial derivatives.
	Realf x[64] = {

		    // values of f(x,y,z) at each corner.
		    _getv(mesh, {{i  , j  , k}}  , rfl),
        _getv(mesh, {{i+1, j  , k}}  , rfl),
        _getv(mesh, {{i  , j+1, k}}  , rfl),
		    _getv(mesh, {{i+1, j+1, k}}  , rfl),
        _getv(mesh, {{i  , j  , k+1}}, rfl),
        _getv(mesh, {{i+1, j  , k+1}}, rfl),
		    _getv(mesh, {{i  , j+1, k+1}}, rfl),
        _getv(mesh, {{i+1, j+1, k+1}}, rfl),

        // values of df/dx at each corner.
		    Realf(0.5)*(_getv(mesh, {{i+1, j  , k   }}, rfl) - _getv(mesh, {{i-1, j  , k   }}, rfl)),
		    Realf(0.5)*(_getv(mesh, {{i+2, j  , k   }}, rfl) - _getv(mesh, {{i  , j  , k   }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+1, j+1, k   }}, rfl) - _getv(mesh, {{i-1, j+1, k   }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+2, j+1, k   }}, rfl) - _getv(mesh, {{i  , j+1, k   }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+1, j  , k+1 }}, rfl) - _getv(mesh, {{i-1, j  , k+1 }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+2, j  , k+1 }}, rfl) - _getv(mesh, {{i  , j  , k+1 }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+1, j+1, k+1 }}, rfl) - _getv(mesh, {{i-1, j+1, k+1 }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+2, j+1, k+1 }}, rfl) - _getv(mesh, {{i  , j+1, k+1 }}, rfl)),
  
        // values of df/dy at each corner.
		    Realf(0.5)*(_getv(mesh, {{i  , j+1 , k  }}, rfl) - _getv(mesh, {{i  , j-1,  k  }}, rfl)),
		    Realf(0.5)*(_getv(mesh, {{i+1, j+1 , k  }}, rfl) - _getv(mesh, {{i+1, j-1,  k  }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i  , j+2 , k  }}, rfl) - _getv(mesh, {{i  , j,    k  }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+1, j+2 , k  }}, rfl) - _getv(mesh, {{i+1, j,    k  }}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i  , j+1 , k+1}}, rfl) - _getv(mesh, {{i  , j- 1, k+1}}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+1, j+1 , k+1}}, rfl) - _getv(mesh, {{i+1, j- 1, k+1}}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i  , j+2 , k+1}}, rfl) - _getv(mesh, {{i  , j,    k+1}}, rfl)),
			  Realf(0.5)*(_getv(mesh, {{i+1, j+2 , k+1}}, rfl) - _getv(mesh, {{i+1, j,    k+1}}, rfl)),

        // values of df/dz at each corner.
		    Realf(0.5)*(_getv(mesh , {{i   , j   , k+1}} , rfl) - _getv(mesh , {{i  , j  , k- 1}}, rfl)),
		    Realf(0.5)*(_getv(mesh , {{i+1 , j   , k+1}} , rfl) - _getv(mesh , {{i+1, j  , k- 1}}, rfl)),
			  Realf(0.5)*(_getv(mesh , {{i   , j+1 , k+1}} , rfl) - _getv(mesh , {{i  , j+1, k- 1}}, rfl)),
			  Realf(0.5)*(_getv(mesh , {{i+1 , j+1 , k+1}} , rfl) - _getv(mesh , {{i+1, j+1, k- 1}}, rfl)),
			  Realf(0.5)*(_getv(mesh , {{i   , j   , k+2}} , rfl) - _getv(mesh , {{i  , j  , k   }}, rfl)),
			  Realf(0.5)*(_getv(mesh , {{i+1 , j   , k+2}} , rfl) - _getv(mesh , {{i+1, j  , k   }}, rfl)),
			  Realf(0.5)*(_getv(mesh , {{i   , j+1 , k+2}} , rfl) - _getv(mesh , {{i  , j+1, k   }}, rfl)),
			  Realf(0.5)*(_getv(mesh , {{i+1 , j+1 , k+2}} , rfl) - _getv(mesh , {{i+1, j+1, k   }}, rfl)),

        // values of d2f/dxdy at each corner.
		    Realf(0.25)*(_getv(mesh, {{i+1,j+1,k}},   rfl) - _getv(mesh, {{i-1,j+1,k}}, rfl)   - _getv(mesh, {{i+1,j-1,k}}, rfl)   + _getv(mesh, {{i-1,j-1,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j+1,k}},   rfl) - _getv(mesh, {{i,j+1,k}}, rfl)     - _getv(mesh, {{i+2,j-1,k}}, rfl)   + _getv(mesh, {{i,j-1,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+2,k}},   rfl) - _getv(mesh, {{i-1,j+2,k}}, rfl)   - _getv(mesh, {{i+1,j,k}}, rfl)     + _getv(mesh, {{i-1,j,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j+2,k}},   rfl) - _getv(mesh, {{i,j+2,k}}, rfl)     - _getv(mesh, {{i+2,j,k}}, rfl)     + _getv(mesh, {{i,j,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+1,k+1}}, rfl) - _getv(mesh, {{i-1,j+1,k+1}}, rfl) - _getv(mesh, {{i+1,j-1,k+1}}, rfl) + _getv(mesh, {{i-1,j-1,k+1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j+1,k+1}}, rfl) - _getv(mesh, {{i,j+1,k+1}}, rfl)   - _getv(mesh, {{i+2,j-1,k+1}}, rfl) + _getv(mesh, {{i,j-1,k+1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+2,k+1}}, rfl) - _getv(mesh, {{i-1,j+2,k+1}}, rfl) - _getv(mesh, {{i+1,j,k+1}}, rfl)   + _getv(mesh, {{i-1,j,k+1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j+2,k+1}}, rfl) - _getv(mesh, {{i,j+2,k+1}}, rfl)   - _getv(mesh, {{i+2,j,k+1}}, rfl)   + _getv(mesh, {{i,j,k+1}}, rfl)),

        // values of d2f/dxdz at each corner.
		    Realf(0.25)*(_getv(mesh, {{i+1,j,k+1}}, rfl)   - _getv(mesh, {{i-1,j,k+1}}, rfl)  - _getv(mesh, {{i+1,j,k-1}}, rfl)   + _getv(mesh, {{i-1,j,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j,k+1}}, rfl)   - _getv(mesh, {{i,j,k+1}}, rfl)    - _getv(mesh, {{i+2,j,k-1}}, rfl)   + _getv(mesh, {{i,j,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+1,k+1}}, rfl) - _getv(mesh, {{i-1,j+1,k+1}}, rfl)- _getv(mesh, {{i+1,j+1,k-1}}, rfl) + _getv(mesh, {{i-1,j+1,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j+1,k+1}}, rfl) - _getv(mesh, {{i,j+1,k+1}}, rfl)  - _getv(mesh, {{i+2,j+1,k-1}}, rfl) + _getv(mesh, {{i,j+1,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j,k+2}}, rfl)   - _getv(mesh, {{i-1,j,k+2}}, rfl)  - _getv(mesh, {{i+1,j,k}}, rfl)     + _getv(mesh, {{i-1,j,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j,k+2}}, rfl)   - _getv(mesh, {{i,j,k+2}}, rfl)    - _getv(mesh, {{i+2,j,k}}, rfl)     + _getv(mesh, {{i,j,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+1,k+2}}, rfl) - _getv(mesh, {{i-1,j+1,k+2}}, rfl)- _getv(mesh, {{i+1,j+1,k}}, rfl)   + _getv(mesh, {{i-1,j+1,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+2,j+1,k+2}}, rfl) - _getv(mesh, {{i,j+1,k+2}}, rfl)  - _getv(mesh, {{i+2,j+1,k}}, rfl)   + _getv(mesh, {{i,j+1,k}}, rfl)),

        // values of d2f/dydz at each corner.
		    Realf(0.25)*(_getv(mesh, {{i,j+1,k+1}}, rfl)   - _getv(mesh, {{i,j-1,k+1}}, rfl)  - _getv(mesh, {{i,j+1,k-1}}, rfl)   + _getv(mesh, {{i,j-1,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+1,k+1}}, rfl) - _getv(mesh, {{i+1,j-1,k+1}}, rfl)- _getv(mesh, {{i+1,j+1,k-1}}, rfl) + _getv(mesh, {{i+1,j-1,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i,j+2,k+1}}, rfl)   - _getv(mesh, {{i,j,k+1}}, rfl)    - _getv(mesh, {{i,j+2,k-1}}, rfl)   + _getv(mesh, {{i,j,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+2,k+1}}, rfl) - _getv(mesh, {{i+1,j,k+1}}, rfl)  - _getv(mesh, {{i+1,j+2,k-1}}, rfl) + _getv(mesh, {{i+1,j,k-1}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i,j+1,k+2}}, rfl)   - _getv(mesh, {{i,j-1,k+2}}, rfl)  - _getv(mesh, {{i,j+1,k}}, rfl)     + _getv(mesh, {{i,j-1,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+1,k+2}}, rfl) - _getv(mesh, {{i+1,j-1,k+2}}, rfl)- _getv(mesh, {{i+1,j+1,k}}, rfl)   + _getv(mesh, {{i+1,j-1,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i,j+2,k+2}}, rfl)   - _getv(mesh, {{i,j,k+2}}, rfl)    - _getv(mesh, {{i,j+2,k}}, rfl)     + _getv(mesh, {{i,j,k}}, rfl)),
			  Realf(0.25)*(_getv(mesh, {{i+1,j+2,k+2}}, rfl) - _getv(mesh, {{i+1,j,k+2}}, rfl)  - _getv(mesh, {{i+1,j+2,k}}, rfl)   + _getv(mesh, {{i+1,j,k}}, rfl)),

			  // values of d3f/dxdydz at each corner.
		    Realf(0.125)*(_getv(mesh, {{i+1,j+1,k+1}}, rfl) - _getv(mesh, {{i-1,j+1,k+1}}, rfl) - _getv(mesh, {{i+1,j-1,k+1}}, rfl) + _getv(mesh, {{i-1,j-1,k+1}}, rfl) 
                 -_getv(mesh, {{i+1,j+1,k-1}}, rfl) + _getv(mesh, {{i-1,j+1,k-1}}, rfl) + _getv(mesh, {{i+1,j-1,k-1}}, rfl) - _getv(mesh, {{i-1,j-1,k-1}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+2,j+1,k+1}}, rfl) - _getv(mesh, {{i,j+1,k+1}}, rfl)   - _getv(mesh, {{i+2,j-1,k+1}}, rfl) + _getv(mesh, {{i,j-1,k+1}}, rfl)
                 -_getv(mesh, {{i+2,j+1,k-1}}, rfl) + _getv(mesh, {{i,j+1,k-1}}, rfl)   + _getv(mesh, {{i+2,j-1,k-1}}, rfl) - _getv(mesh, {{i,j-1,k-1}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+1,j+2,k+1}}, rfl) - _getv(mesh, {{i-1,j+2,k+1}}, rfl) - _getv(mesh, {{i+1,j,k+1}}, rfl)   + _getv(mesh, {{i-1,j,k+1}}, rfl)
                 -_getv(mesh, {{i+1,j+2,k-1}}, rfl) + _getv(mesh, {{i-1,j+2,k-1}}, rfl) + _getv(mesh, {{i+1,j,k-1}}, rfl)   - _getv(mesh, {{i-1,j,k-1}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+2,j+2,k+1}}, rfl) - _getv(mesh, {{i,j+2,k+1}}, rfl)   - _getv(mesh, {{i+2,j,k+1}}, rfl)   + _getv(mesh, {{i,j,k+1}}, rfl)
                 -_getv(mesh, {{i+2,j+2,k-1}}, rfl) + _getv(mesh, {{i,j+2,k-1}}, rfl)   + _getv(mesh, {{i+2,j,k-1}}, rfl)   - _getv(mesh, {{i,j,k-1}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+1,j+1,k+2}}, rfl) - _getv(mesh, {{i-1,j+1,k+2}}, rfl) - _getv(mesh, {{i+1,j-1,k+2}}, rfl) + _getv(mesh, {{i-1,j-1,k+2}}, rfl)
                 -_getv(mesh, {{i+1,j+1,k}}, rfl)   + _getv(mesh, {{i-1,j+1,k}}, rfl)   + _getv(mesh, {{i+1,j-1,k}}, rfl)   - _getv(mesh, {{i-1,j-1,k}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+2,j+1,k+2}}, rfl) - _getv(mesh, {{i,j+1,k+2}}, rfl)   - _getv(mesh, {{i+2,j-1,k+2}}, rfl) + _getv(mesh, {{i,j-1,k+2}}, rfl)
                 -_getv(mesh, {{i+2,j+1,k}}, rfl)   + _getv(mesh, {{i,j+1,k}}, rfl)     + _getv(mesh, {{i+2,j-1,k}}, rfl)   - _getv(mesh, {{i,j-1,k}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+1,j+2,k+2}}, rfl) - _getv(mesh, {{i-1,j+2,k+2}}, rfl) - _getv(mesh, {{i+1,j,k+2}}, rfl)   + _getv(mesh, {{i-1,j,k+2}}, rfl)
                 -_getv(mesh, {{i+1,j+2,k}}, rfl)   + _getv(mesh, {{i-1,j+2,k}}, rfl)   + _getv(mesh, {{i+1,j,k}}, rfl)     - _getv(mesh, {{i-1,j,k}}, rfl)),
			  Realf(0.125)*(_getv(mesh, {{i+2,j+2,k+2}}, rfl) - _getv(mesh, {{i,j+2,k+2}}, rfl)   - _getv(mesh, {{i+2,j,k+2}}, rfl)   + _getv(mesh, {{i,j,k+2}}, rfl)
              -_getv(mesh, {{i+2,j+2,k}}, rfl)   + _getv(mesh, {{i,j+2,k}}, rfl)     + _getv(mesh, {{i+2,j,k}}, rfl)     - _getv(mesh, {{i,j,k}}, rfl))
		};


  // Convert voxel values and partial derivatives to interpolation coefficients.
  Realf _coefs[64];
  for (int i=0;i<64;++i) {
    _coefs[i] = 0.0;
    for (int j=0;j<64;++j) {
      _coefs[i] += _C[i][j]*x[j];
    }
  }

  int ijkn(0);
  Realf dzpow(1);
  Realf f0(0);
  Realf alpha;

  /*
  grad[0] = T(0);
  grad[1] = T(0);
  grad[2] = T(0);
  */

  for(int k = 0; k < 4; ++k) {
    Realf dypow(1);
    for(int j = 0; j < 4; ++j) {
      alpha = dypow*dzpow*
        (_coefs[ijkn] + dx*(_coefs[ijkn+1] + dx*(_coefs[ijkn+2] + dx*_coefs[ijkn+3])));
      f0 += alpha;

      /*
      grad[0] += dypow*dzpow*
        (_coefs[ijkn+1] + dx*(2*_coefs[ijkn+2] + dx*3*_coefs[ijkn+3]));

      if(j > 0) {
        grad[1] += j*alpha/dy;
      }

      if(k > 0) {
        grad[2] += k*alpha/dz;
      }
      */

      ijkn += 4;
      dypow *= dy;
    }
    dzpow *= dz;
  }

  return f0;
}





} // end of namespace toolbox
