#pragma once

#include <vector>
#include <algorithm>
  
// for GPU supported sorting
#include "iter/dynArray.h"
#include "iter/devcall.h"

namespace toolbox {


template <typename T>
DEVCALLABLE int find_sorted_nearest( ManVec<T> const& arr, T val) {

    if(arr.size() == 0) return 0;

    //first and smallest value
    if(val < arr[0]) return 0;

    // value is between array elements
    for(int i=1; i<arr.size(); i++){
      if( ( arr[i-1] <= val ) && ( val < arr[i] ) ) {
        return i;
      }
    }

    // last element
    return static_cast<int>( arr.size() ) - 1; // len(arr)-1 
}

// find element from sorted array
// based on https://en.algorithmica.org/hpc/data-structures/binary-search/
template <typename T>
inline int find_sorted_nearest_algo2(T* t, T x, const int n) 
{
  // limits
  if(x <= t[0]  ) return 0;
  if(x >= t[n-1]) return n;

  float *base = t;
  int len = n;
  while (len > 1) {
        int half = len / 2;
        base += (base[half - 1] < x) * half; // will be replaced with a "cmov"
        len -= half;
    }
  //return *base; // value
  return base - &t[0]; // index
}


template <typename T>
inline size_t find_rev_sorted_nearest( ManVec<T> const& arr, const T val) {

  const size_t arr_size = arr.size();
  if(arr_size == 0) return 0;

  //first and largest value
  //if(arr[0] < val) return 0;
  if(val > arr[0]) return 0;

  // value is between array elements
  size_t i;
  for(i=1; i<arr_size; i++){
    //if( ( arr[i-1] >= val ) && ( val > arr[i] ) ) { // TODO equal to below?
    if( val > arr[i] ) { // TODO not strictly equivalent to below; more aggressive
      return i;
    }
  }

  // last element
  return arr_size - 1; // len(arr)-1 
}


// reverse find of reverse sorted array
template <typename T>
inline size_t revfind_rev_sorted_nearest( ManVec<T> const& arr, const T val) {

  if(arr.size() == 0) return 0;

  // special case where we would have to wrap all around the array
  if(val > arr[0u]) return 0u;

  const size_t N = arr.size()-1;
  if(val <= arr[N]) return N;

  for(size_t i=N; i>=0u; i--){
    //fmt::print("    m3: iterating {} comparing {} < {}\n",i, val, arr[i]);
    if( val <= arr[i] ) return i+1;
  }

  return 0u;
}

// reverse find of reverse sorted array
template <typename T>
inline size_t find_rev_sorted_nearest_algo( ManVec<T> & arr, const T val) {

  // NOTE profiling shows that these do not improve performance
  //if(val > arr[0u]) return 0u;
  //if(val < arr[arr.size()-1u]) return arr.size();

  auto it = std::lower_bound( arr.begin(), arr.end(), val, std::greater_equal<T>() );
  return it - arr.begin();
}


// reverse find of reverse sorted array
// based on https://en.algorithmica.org/hpc/data-structures/binary-search/
template <typename T>
inline int find_rev_sorted_nearest_algo2( 
    ManVec<T> & t, 
    const T x) 
{
  const size_t n = t.size();
  if(x < t[n-1]) return n;

  float *base = t.data();
  int len = n;
  while (len > 1) {
        int half = len / 2;
        len -= half;

        // prefetch hinting is about 5% slower
        //__builtin_prefetch(&base[len / 2 - 1]);
        //__builtin_prefetch(&base[half + len / 2 - 1]);

        base += (base[half - 1] >= x) * half; // will be replaced with a "cmov"
    }
  //return *base; // value
  return base - &t[0]; // index
}


// Sample between [imin, imax[
template <typename T>
inline int sample_prob( ManVec<T> const& arr, const T val) {

  int ind = 0;;
  for(size_t i=0; i<arr.size(); i++)
  {
    if( val >= arr[i] ) { 
      ind++;
    } else {
      break;
    }
  }

  return ind;
}

// Sample between [imin, imax[
template <typename T>
int inline sample_prob_between( ManVec<T> const& ws, const T val, const int imin, const int imax) {

    T wmin, wmax;
    if(imin == 0) {
        wmin = 0.0; // # equal to ws[-1]
    } else {
        wmin = ws[imin-1];
    }
    wmax = ws[imax-1];
    //const wminmax = wmax - wmin;
    const T vminmax = val*(wmax - wmin) + wmin;

    for(int j=imin; j<imax; j++){
        //if( val < (ws[j] - wmin)/(wmax - wmin) ) return j;
        //if( val < (ws[j] - wmin)/wminmax ) return j;
        if( vminmax < ws[j] ) return j;
    }
    return imax;
}


// Sample between [imin, imax[
template <typename T>
size_t inline sample_prob_between_algo( 
    const ManVec<T>& ws, 
    const T val, 
    const size_t imin, 
    const size_t imax) 
{
    T wmin, wmax;
    if(imin == 0u) {
        wmin = static_cast<T>(0); // # equal to ws[-1]
    } else {
        wmin = ws[imin-1u];
    }
    wmax = ws[imax-1u];
    const T vminmax = val*(wmax - wmin) + wmin;

    // search in the range [imin, imax] using binary search
    auto beg = std::next(ws.cbegin(), imin);
    auto end = std::next(ws.cbegin(), imax);
    auto it = std::upper_bound(beg, end, vminmax);

    return it - ws.cbegin();
}



} // enf of ns toolbox
