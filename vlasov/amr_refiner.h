#pragma once

#include <cmath> 
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

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
    cells_to_unrefine.clear();


    T refine_indicator   = T(0),
      unrefine_indicator = T(0);
    

    for(const auto& cid : mesh.get_cells(true)) {

      if (!mesh.is_leaf(cid)) continue;

      // error indicator for refinement/unrefinement
        
      // error_indicator = maximum_gradient(mesh, cid);
      refine_indicator = maximum_value(mesh, cid);
      // error_indicator = relative_difference(mesh, cid);

      unrefine_indicator = refine_indicator;


      // to be refined
      if( refine_indicator > tolerance ) {
        cells_to_refine.insert(cid);
        continue;
      }

      // if we are at the bottom, no need to unrefine
      if (mesh.get_refinement_level(cid) < 1) continue;


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
  }



  void refine( AdaptiveMesh<T,D>& mesh )
  {
    cells_created.clear();

    for(auto cid : cells_to_refine) {

      T parent_value = mesh.get(cid);

      // creating empty cells 
      for(auto cidc : mesh.get_children(cid)) {
        mesh.set(cidc, parent_value);

        cells_created.push_back(cidc);
      }

    }
  }



  void unrefine( AdaptiveMesh<T,D>& mesh )
  {
    cells_removed.clear();

    // NOTE: these are actually parents of the children to be removed
    for(const auto cid : cells_to_unrefine) {

        auto children = mesh.get_children(cid);

        // collect cell values and remove
        T avg_value = T(0);
        for(auto cidc : children) {
          avg_value += mesh.get(cidc);

          mesh.data.erase(cidc);
          cells_removed.push_back(cidc);
        }
        avg_value /= T( children.size() ); // this is the avg of children vals

        mesh.set(cid, avg_value);
    }
  }






};







} // end of namespace toolbox






