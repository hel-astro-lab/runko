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
    i = indices[0], // + uint64_t(coordinates[0] - 1.5),
    j = indices[1], // + uint64_t(coordinates[1] - 1.5),
    k = indices[2]; // + uint64_t(coordinates[2] - 1.5);
	
	T dx = coordinates[0] - T(0.5); 
  T dy = coordinates[1] - T(0.5); 
  T dz = coordinates[2] - T(0.5);
	
	T d00 = lerp(dx, _getv(mesh, {{i, j,   k  }}, rfl), _getv(mesh, {{i+1, j,   k  }}, rfl) );
	T d10 = lerp(dx, _getv(mesh, {{i, j+1, k  }}, rfl), _getv(mesh, {{i+1, j+1, k  }}, rfl) );
	T d01 = lerp(dx, _getv(mesh, {{i, j,   k+1}}, rfl), _getv(mesh, {{i+1, j,   k+1}}, rfl) );
	T d11 = lerp(dx, _getv(mesh, {{i, j+1, k+1}}, rfl), _getv(mesh, {{i+1, j+1, k+1}}, rfl) );
	T d0  = lerp(dy, d00, d10);
	T d1  = lerp(dy, d01, d11);
	T d   = lerp(dz, d0,  d1);
	
	return d;
}



/// \brief Tricubic interpolation
//
// This is heavily optimized with all the element positioning and looping.
// 
// Transforms from cell-centered values to vertex-centered,
// then uses the normal lerp routine.
//
// Partially based on:
// https://svn.blender.org/svnroot/bf-blender/
// branches/volume25/source/blender/blenlib/intern/voxel.c  
//
// and here:
//
// https://github.com/igmhub/likely/blob/master/likely/TriCubicInterpolator.cc
//
template<typename T>
T tricubic_interp(
    const AdaptiveMesh<T,3>& mesh,
    typename AdaptiveMesh<T,3>::indices_t& indices,
    typename AdaptiveMesh<T,3>::value_array_t coordinates,
    int rfl)
{
  uint64_t 
    i = indices[0], // + uint64_t(coordinates[0] - 1.5),
    j = indices[1], // + uint64_t(coordinates[1] - 1.5),
    k = indices[2]; // + uint64_t(coordinates[2] - 1.5);
	
	T dx = coordinates[0] - T(0.5); 
  T dy = coordinates[1] - T(0.5); 
  T dz = coordinates[2] - T(0.5);

  // Extract the local vocal values and calculate partial derivatives.
	double x[64] = {

		    // values of f(x,y,z) at each corner.
		    _data[_index(xi,yi,zi)],
        _data[_index(xi+1,yi,zi)],
        _data[_index(xi,yi+1,zi)],
		    _data[_index(xi+1,yi+1,zi)],
        _data[_index(xi,yi,zi+1)],
        _data[_index(xi+1,yi,zi+1)],
		    _data[_index(xi,yi+1,zi+1)],
        _data[_index(xi+1,yi+1,zi+1)],

        // values of df/dx at each corner.
		    0.5*(_data[_index(xi+1,yi,zi)]     - _data[_index(xi-1,yi,zi)]),
		    0.5*(_data[_index(xi+2,yi,zi)]     - _data[_index(xi,yi,zi)]),
			  0.5*(_data[_index(xi+1,yi+1,zi)]   - _data[_index(xi-1,yi+1,zi)]),
			  0.5*(_data[_index(xi+2,yi+1,zi)]   - _data[_index(xi,yi+1,zi)]),
			  0.5*(_data[_index(xi+1,yi,zi+1)]   - _data[_index(xi-1,yi,zi+1)]),
			  0.5*(_data[_index(xi+2,yi,zi+1)]   - _data[_index(xi,yi,zi+1)]),
			  0.5*(_data[_index(xi+1,yi+1,zi+1)] - _data[_index(xi-1,yi+1,zi+1)]),
			  0.5*(_data[_index(xi+2,yi+1,zi+1)] - _data[_index(xi,yi+1,zi+1)]),

        // values of df/dy at each corner.
		    0.5*(_data[_index(xi,yi+1,zi)]     - _data[_index(xi,yi-1,zi)]),
		    0.5*(_data[_index(xi+1,yi+1,zi)]   - _data[_index(xi+1,yi-1,zi)]),
			  0.5*(_data[_index(xi,yi+2,zi)]     - _data[_index(xi,yi,zi)]),
			  0.5*(_data[_index(xi+1,yi+2,zi)]   - _data[_index(xi+1,yi,zi)]),
			  0.5*(_data[_index(xi,yi+1,zi+1)]   - _data[_index(xi,yi-1,zi+1)]),
			  0.5*(_data[_index(xi+1,yi+1,zi+1)] - _data[_index(xi+1,yi-1,zi+1)]),
			  0.5*(_data[_index(xi,yi+2,zi+1)]   - _data[_index(xi,yi,zi+1)]),
			  0.5*(_data[_index(xi+1,yi+2,zi+1)] - _data[_index(xi+1,yi,zi+1)]),

        // values of df/dz at each corner.
		    0.5*(_data[_index(xi,yi,zi+1)]     -_data[_index(xi,yi,zi-1)]),
		    0.5*(_data[_index(xi+1,yi,zi+1)]   -_data[_index(xi+1,yi,zi-1)]),
			  0.5*(_data[_index(xi,yi+1,zi+1)]   -_data[_index(xi,yi+1,zi-1)]),
			  0.5*(_data[_index(xi+1,yi+1,zi+1)] -_data[_index(xi+1,yi+1,zi-1)]),
			  0.5*(_data[_index(xi,yi,zi+2)]     -_data[_index(xi,yi,zi)]),
			  0.5*(_data[_index(xi+1,yi,zi+2)]   -_data[_index(xi+1,yi,zi)]),
			  0.5*(_data[_index(xi,yi+1,zi+2)]   -_data[_index(xi,yi+1,zi)]),
			  0.5*(_data[_index(xi+1,yi+1,zi+2)] -_data[_index(xi+1,yi+1,zi)]),

        // values of d2f/dxdy at each corner.
		    0.25*(_data[_index(xi+1,yi+1,zi)]   - _data[_index(xi-1,yi+1,zi)]   - _data[_index(xi+1,yi-1,zi)]   + _data[_index(xi-1,yi-1,zi)]),
			  0.25*(_data[_index(xi+2,yi+1,zi)]   - _data[_index(xi,yi+1,zi)]     - _data[_index(xi+2,yi-1,zi)]   + _data[_index(xi,yi-1,zi)]),
			  0.25*(_data[_index(xi+1,yi+2,zi)]   - _data[_index(xi-1,yi+2,zi)]   - _data[_index(xi+1,yi,zi)]     + _data[_index(xi-1,yi,zi)]),
			  0.25*(_data[_index(xi+2,yi+2,zi)]   - _data[_index(xi,yi+2,zi)]     - _data[_index(xi+2,yi,zi)]     + _data[_index(xi,yi,zi)]),
			  0.25*(_data[_index(xi+1,yi+1,zi+1)] - _data[_index(xi-1,yi+1,zi+1)] - _data[_index(xi+1,yi-1,zi+1)] + _data[_index(xi-1,yi-1,zi+1)]),
			  0.25*(_data[_index(xi+2,yi+1,zi+1)] - _data[_index(xi,yi+1,zi+1)]   - _data[_index(xi+2,yi-1,zi+1)] + _data[_index(xi,yi-1,zi+1)]),
			  0.25*(_data[_index(xi+1,yi+2,zi+1)] - _data[_index(xi-1,yi+2,zi+1)] - _data[_index(xi+1,yi,zi+1)]   + _data[_index(xi-1,yi,zi+1)]),
			  0.25*(_data[_index(xi+2,yi+2,zi+1)] - _data[_index(xi,yi+2,zi+1)]   - _data[_index(xi+2,yi,zi+1)]   + _data[_index(xi,yi,zi+1)]),

        // values of d2f/dxdz at each corner.
		    0.25*(_data[_index(xi+1,yi,zi+1)]   - _data[_index(xi-1,yi,zi+1)]  - _data[_index(xi+1,yi,zi-1)]   + _data[_index(xi-1,yi,zi-1)]),
			  0.25*(_data[_index(xi+2,yi,zi+1)]   - _data[_index(xi,yi,zi+1)]    - _data[_index(xi+2,yi,zi-1)]   + _data[_index(xi,yi,zi-1)]),
			  0.25*(_data[_index(xi+1,yi+1,zi+1)] - _data[_index(xi-1,yi+1,zi+1)]- _data[_index(xi+1,yi+1,zi-1)] + _data[_index(xi-1,yi+1,zi-1)]),
			  0.25*(_data[_index(xi+2,yi+1,zi+1)] - _data[_index(xi,yi+1,zi+1)]  - _data[_index(xi+2,yi+1,zi-1)] + _data[_index(xi,yi+1,zi-1)]),
			  0.25*(_data[_index(xi+1,yi,zi+2)]   - _data[_index(xi-1,yi,zi+2)]  - _data[_index(xi+1,yi,zi)]     + _data[_index(xi-1,yi,zi)]),
			  0.25*(_data[_index(xi+2,yi,zi+2)]   - _data[_index(xi,yi,zi+2)]    - _data[_index(xi+2,yi,zi)]     + _data[_index(xi,yi,zi)]),
			  0.25*(_data[_index(xi+1,yi+1,zi+2)] - _data[_index(xi-1,yi+1,zi+2)]- _data[_index(xi+1,yi+1,zi)]   + _data[_index(xi-1,yi+1,zi)]),
			  0.25*(_data[_index(xi+2,yi+1,zi+2)] - _data[_index(xi,yi+1,zi+2)]  - _data[_index(xi+2,yi+1,zi)]   + _data[_index(xi,yi+1,zi)]),

        // values of d2f/dydz at each corner.
		    0.25*(_data[_index(xi,yi+1,zi+1)]   - _data[_index(xi,yi-1,zi+1)]  - _data[_index(xi,yi+1,zi-1)]   + _data[_index(xi,yi-1,zi-1)]),
			  0.25*(_data[_index(xi+1,yi+1,zi+1)] - _data[_index(xi+1,yi-1,zi+1)]- _data[_index(xi+1,yi+1,zi-1)] + _data[_index(xi+1,yi-1,zi-1)]),
			  0.25*(_data[_index(xi,yi+2,zi+1)]   - _data[_index(xi,yi,zi+1)]    - _data[_index(xi,yi+2,zi-1)]   + _data[_index(xi,yi,zi-1)]),
			  0.25*(_data[_index(xi+1,yi+2,zi+1)] - _data[_index(xi+1,yi,zi+1)]  - _data[_index(xi+1,yi+2,zi-1)] + _data[_index(xi+1,yi,zi-1)]),
			  0.25*(_data[_index(xi,yi+1,zi+2)]   - _data[_index(xi,yi-1,zi+2)]  - _data[_index(xi,yi+1,zi)]     + _data[_index(xi,yi-1,zi)]),
			  0.25*(_data[_index(xi+1,yi+1,zi+2)] - _data[_index(xi+1,yi-1,zi+2)]- _data[_index(xi+1,yi+1,zi)]   + _data[_index(xi+1,yi-1,zi)]),
			  0.25*(_data[_index(xi,yi+2,zi+2)]   - _data[_index(xi,yi,zi+2)]    - _data[_index(xi,yi+2,zi)]     + _data[_index(xi,yi,zi)]),
			  0.25*(_data[_index(xi+1,yi+2,zi+2)] - _data[_index(xi+1,yi,zi+2)]  - _data[_index(xi+1,yi+2,zi)]   + _data[_index(xi+1,yi,zi)]),

			  // values of d3f/dxdydz at each corner.
		    0.125*(_data[_index(xi+1,yi+1,zi+1)] - _data[_index(xi-1,yi+1,zi+1)] - _data[_index(xi+1,yi-1,zi+1)] + _data[_index(xi-1,yi-1,zi+1)] 
              -_data[_index(xi+1,yi+1,zi-1)] + _data[_index(xi-1,yi+1,zi-1)] + _data[_index(xi+1,yi-1,zi-1)] - _data[_index(xi-1,yi-1,zi-1)]),
			  0.125*(_data[_index(xi+2,yi+1,zi+1)] - _data[_index(xi,yi+1,zi+1)]   - _data[_index(xi+2,yi-1,zi+1)] + _data[_index(xi,yi-1,zi+1)]
              -_data[_index(xi+2,yi+1,zi-1)] + _data[_index(xi,yi+1,zi-1)]   + _data[_index(xi+2,yi-1,zi-1)] - _data[_index(xi,yi-1,zi-1)]),
			  0.125*(_data[_index(xi+1,yi+2,zi+1)] - _data[_index(xi-1,yi+2,zi+1)] - _data[_index(xi+1,yi,zi+1)]   + _data[_index(xi-1,yi,zi+1)]
              -_data[_index(xi+1,yi+2,zi-1)] + _data[_index(xi-1,yi+2,zi-1)] + _data[_index(xi+1,yi,zi-1)]   - _data[_index(xi-1,yi,zi-1)]),
			  0.125*(_data[_index(xi+2,yi+2,zi+1)] - _data[_index(xi,yi+2,zi+1)]   - _data[_index(xi+2,yi,zi+1)]   + _data[_index(xi,yi,zi+1)]
              -_data[_index(xi+2,yi+2,zi-1)] + _data[_index(xi,yi+2,zi-1)]   + _data[_index(xi+2,yi,zi-1)]   - _data[_index(xi,yi,zi-1)]),
			  0.125*(_data[_index(xi+1,yi+1,zi+2)] - _data[_index(xi-1,yi+1,zi+2)] - _data[_index(xi+1,yi-1,zi+2)] + _data[_index(xi-1,yi-1,zi+2)]
              -_data[_index(xi+1,yi+1,zi)]   + _data[_index(xi-1,yi+1,zi)]   + _data[_index(xi+1,yi-1,zi)]   - _data[_index(xi-1,yi-1,zi)]),
			  0.125*(_data[_index(xi+2,yi+1,zi+2)] - _data[_index(xi,yi+1,zi+2)]   - _data[_index(xi+2,yi-1,zi+2)] + _data[_index(xi,yi-1,zi+2)]
              -_data[_index(xi+2,yi+1,zi)]   + _data[_index(xi,yi+1,zi)]     + _data[_index(xi+2,yi-1,zi)]   - _data[_index(xi,yi-1,zi)]),
			  0.125*(_data[_index(xi+1,yi+2,zi+2)] - _data[_index(xi-1,yi+2,zi+2)] - _data[_index(xi+1,yi,zi+2)]   + _data[_index(xi-1,yi,zi+2)]
              -_data[_index(xi+1,yi+2,zi)]   + _data[_index(xi-1,yi+2,zi)]   + _data[_index(xi+1,yi,zi)]     - _data[_index(xi-1,yi,zi)]),
			  0.125*(_data[_index(xi+2,yi+2,zi+2)] - _data[_index(xi,yi+2,zi+2)]   - _data[_index(xi+2,yi,zi+2)]   + _data[_index(xi,yi,zi+2)]
              -_data[_index(xi+2,yi+2,zi)]   + _data[_index(xi,yi+2,zi)]     + _data[_index(xi+2,yi,zi)]     - _data[_index(xi,yi,zi)])
		};




}





} // end of namespace toolbox
