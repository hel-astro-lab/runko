#pragma once

#include <cmath> 
#include <algorithm>
#include <stdexcept>

#include "mesh.h"
#include "amr_numerics.h"


using std::max_element;
using std::abs;



namespace toolbox {


template<typename T, int D>
class Adapter {

  public:

  typedef typename AdaptiveMesh<T,D>::indices_t indices_t;
  typedef typename AdaptiveMesh<T,D>::value_array_t value_array_t;
  
  std::vector<uint64_t> cells_to_refine;
  std::vector<uint64_t> cells_to_unrefine;
  std::vector<uint64_t> cells_created;


  /// General tolerance for refinement
  T tolerance = 1.0e-4;
  T maximum_data_value = T(1.0);


  void set_maximum_data_value(T val) {
    maximum_data_value = val;
  }


  // T relative_difference(const AdaptiveMesh<T,D>& mesh, const uint64_t cid) const { }



  T maximum_value(const AdaptiveMesh<T,D>& mesh, const uint64_t cid) const 
  {
    int rfl = mesh.get_refinement_level(cid);

    T val = mesh.get_from_roots(cid);

    value_array_t lens = mesh.get_length(rfl);

    // compute volume 
    T vol = T(1);
    for(auto x : lens) vol *= x;

    return val * vol / maximum_data_value;
  }



  // error norm is max{ |grad(f)| }
  T maximum_gradient(const AdaptiveMesh<T,D>& mesh, const uint64_t cid) const
  {
    int rfl = mesh.get_refinement_level(cid);
    indices_t ind = mesh.get_indices(cid);

    value_array_t gradient = grad<T,D>(mesh, ind, rfl);

    return *std::max_element(std::begin( gradient) , std::end(gradient),
        [](T a, T b){return std::abs(a) < std::abs(b);});
  }



  void check( AdaptiveMesh<T,D>& mesh )
  {
    cells_to_refine.clear();

    T error_indicator = T(0);
    
    for(const auto& cid : mesh.get_cells(true)) {

      if (!mesh.is_leaf(cid)) continue;

      // error_indicator = maximum_gradient(mesh, cid);
      error_indicator = maximum_value(mesh, cid);

      // error_indicator = relative_difference(mesh, cid);


      if( error_indicator > tolerance ) cells_to_refine.push_back(cid);

      // fmt::print("Checking cell {} ({},{},{}): E: {} \n", cid, ind[0], ind[1], ind[2], error_indicator);
    }
  }



  void refine( AdaptiveMesh<T,D>& mesh )
  {
    cells_created.clear();

    for(const auto& cid : cells_to_refine) {

      T parent_value = mesh.get(cid);

      // creating empty cells 
      for(auto& cidc : mesh.get_children(cid)) {
        mesh.set(cidc, parent_value);

        cells_created.push_back(cidc);
      }

    }


  }






};







} // end of namespace toolbox






