#pragma once

#include <vector>
#include <algorithm>
#include <numeric>

// for GPU supported sorting
#include "iter/dynArray.h"
#include "iter/devcall.h"

// sort vector but return sorted indices 
//
// ref: https://stackoverflow.com/questions/10580982/c-sort-keeping-track-of-indices
template <typename T>
std::vector<size_t> argsort(std::vector<T> const& values) {

  std::vector<size_t> indices( values.size() );

  // linspace from 0 to N 
  std::iota( std::begin(indices), std::end(indices), static_cast<size_t>(0) );
  std::sort(
    begin(indices), end(indices),
      [&](size_t a, size_t b) { return values[a] < values[b]; }
  );
  return indices;
}

// reverse version
template <typename T>
std::vector<size_t> argsort_rev(std::vector<T> const& values) {

  std::vector<size_t> indices( values.size() );

  // linspace from 0 to N 
  std::iota( std::begin(indices), std::end(indices), static_cast<size_t>(0) );
  std::sort(
    begin(indices), end(indices),
      [&](size_t a, size_t b) { return values[a] > values[b]; } // NOTE: only difference is here
  );
  return indices;
}


// GPU version of reverse argsort
// FIXME: should indices vector std::vector<size_t> be ManVec as well?
template <typename T>
DEVCALLABLE ManVec<size_t> argsort_rev(ManVec<T> const& values) {

  size_t N = values.size();

  ManVec<size_t> indices;
  indices.resize(N);

  // linspace from 0 to N 
  std::iota( indices.begin(), indices.end(), 0u );

  // FIXME sort is not devcallable
  std::sort( indices.begin(), indices.end(),
      [&](size_t a, size_t b) { return values[a] > values[b]; } // NOTE: only difference is here
  );

  return indices;
}




