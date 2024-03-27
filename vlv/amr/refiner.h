#pragma once

#include <cmath> 
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

#include "vlv/amr/mesh.h"
#include "vlv/amr/numerics.h"


using std::max_element;
using std::abs;



namespace toolbox {


template<typename T, int V>
class Adapter {

  public:

  using indices_t = typename AdaptiveMesh<T,V>::indices_t;
  using value_array_t = typename AdaptiveMesh<T,V>::value_array_t;
  
  // std::vector<uint64_t> cells_to_refine;
  // std::vector<uint64_t> cells_to_unrefine;
  // std::vector<uint64_t> cells_created;


  std::unordered_set<uint64_t> cells_to_refine;
  std::unordered_set<uint64_t> cells_to_unrefine;
  std::vector<uint64_t>        cells_created;
  std::vector<uint64_t>        cells_removed;


  /// General tolerance for refinement
  T tolerance = 1.0e-4;
  T maximum_data_value = T(1.0);


  void set_maximum_data_value(T val) {
    maximum_data_value = val;
  }


  // T relative_difference(const AdaptiveMesh<T,V>& mesh, const uint64_t cid) const { }



  T maximum_value(const AdaptiveMesh<T,V>& mesh, const uint64_t cid) const 
  {
    int rfl = mesh.get_refinement_level(cid);

    T val = mesh.get_from_roots(cid);

    value_array_t lens = mesh.get_length(rfl);

    // compute volume 
    auto vol = T(1);
    for(auto x : lens) vol *= x;

    return val * vol / maximum_data_value;
  }



  // error norm is max{ |grad(f)| }
  T maximum_gradient(const AdaptiveMesh<T,V>& mesh, const uint64_t cid) const
  {
    int rfl = mesh.get_refinement_level(cid);
    indices_t ind = mesh.get_indices(cid);

    value_array_t gradient = grad<T,V>(mesh, ind, rfl);

    return *std::max_element(std::begin( gradient) , std::end(gradient),
        [](T a, T b){return std::abs(a) < std::abs(b);});
  }


  /// Check given cell for possible refinement; if positive then it is appended to internal lists
  void check_cell(
      AdaptiveMesh<T,V>& mesh,
      uint64_t cid, 
      T refine_indicator, 
      T unrefine_indicator)
  {

    // to be refined
    if( refine_indicator > tolerance ) {
      cells_to_refine.insert(cid);
      return;
    }

    // if we are at the bottom, no need to unrefine
    if (mesh.get_refinement_level(cid) < 1) return;


    // to be possibly unrefined
    if ( unrefine_indicator < 0.1*tolerance ) {
      bool to_be_unrefined = true;

      auto siblings = mesh.get_siblings(cid);

      // check if any siblings are marked for refining
      for(auto cids : siblings) {
        if( cells_to_refine.count(cids) > 0) {
          to_be_unrefined = false;
          break;
        }
      }

      if(to_be_unrefined) {
        cells_to_unrefine.insert( mesh.get_parent(cid) );
      }
    }

  }


  /// Check full mesh for refinement
  void check( AdaptiveMesh<T,V>& mesh )
  {
    cells_to_refine.clear();
    cells_to_unrefine.clear();


    auto refine_indicator   = T(0),
      unrefine_indicator = T(0);
    

    for(const auto& it : mesh.data) {

      if (!mesh.is_leaf(it.first)) continue;

      // error indicator for refinement/unrefinement
        
      // error_indicator = maximum_gradient(mesh, cid);
      refine_indicator = maximum_value(mesh, it.first);
      // error_indicator = relative_difference(mesh, cid);

      unrefine_indicator = refine_indicator;

      check_cell(mesh, it.first, refine_indicator, unrefine_indicator);
    }
  }


  // create empty new leafs 
  void refine( AdaptiveMesh<T,V>& mesh )
  {
    cells_created.clear();

    for(const auto& cid : cells_to_refine) {

      // T parent_value = mesh.get(cid);
        
      // creating empty cells 
      for(const auto& cidc : mesh.get_children(cid)) {
        //mesh.set(cidc, parent_value);
        cells_created.push_back(cidc);
      }
    }
  }



  void unrefine( AdaptiveMesh<T,V>& mesh )
  {
    cells_removed.clear();

    // NOTE: these are actually parents of the children to be removed
    for(const auto& cid : cells_to_unrefine) {

        auto children = mesh.get_children(cid);

        // collect cell values and remove
        auto avg_value = T(0);
        for(const auto& cidc : children) {
          avg_value += mesh.get(cidc);

          mesh.data.erase(cidc);
          cells_removed.push_back(cidc);
        }
        avg_value /= static_cast<T>( children.size() ); // this is the avg of children vals

        mesh.set(cid, avg_value);
    }
  }


};







} // end of namespace toolbox






