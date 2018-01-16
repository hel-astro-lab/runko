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
  T tolerance = 0.001;


  void check( AdaptiveMesh<T,D>& mesh )
  {
    
    for(const auto& cid : mesh.get_cells(true)) {

      int rfl = mesh.get_refinement_level(cid);
      indices_t ind = mesh.get_indices(cid);

      // error norm is max{ |grad(f)| }
      value_array_t gradient = grad<T,D>(mesh, ind, rfl);
      T error_indicator = *std::max_element(std::begin( gradient) , std::end(gradient),
          [](T a, T b){return std::abs(a) < std::abs(b);});

      if( error_indicator < tolerance ) continue;

      // fmt::print("Checking cell {} ({},{},{}): E: {} \n", cid, ind[0], ind[1], ind[2], error_indicator);

      cells_to_refine.push_back(cid);
    }
  }


  void refine( AdaptiveMesh<T,D>& mesh )
  {

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






