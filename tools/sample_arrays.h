#pragma once

#include <vector>
  
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

template <typename T>
inline int find_rev_sorted_nearest( ManVec<T> const& arr, const T val) {
  
  const int arr_size = static_cast<int>( arr.size() );
  if(arr_size == 0) return 0;

  //first and largest value
  if(arr[0] < val) return 0;

  // value is between array elements
  for(int i=1; i<arr_size; i++){
    //if( ( arr[i-1] >= val ) && ( val > arr[i] ) ) { // TODO equal to below?
    if( val > arr[i] ) { // TODO not strictly equivalent to below; more aggressive
      return i;
    }
  }

  // last element
  return arr_size - 1; // len(arr)-1 
}

// Sample between [imin, imax[
template <typename T>
DEVCALLABLE int sample_prob( ManVec<T> const& arr, T val) {

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
    return imax-1;
}


} // enf of ns toolbox
