#pragma once

#include <vector>
#include <array>
#include <cmath>

// for GPU supported sorting
#include "iter/dynArray.h"
#include "iter/devcall.h"


namespace toolbox {

template <typename T>
DEVCALLABLE void arange(T start, T stop, ManVec<T>& ret)
{
  //ManVec<T> ret;

  size_t num = static_cast<size_t>( stop - start );
  ret.resize(num);

  for(size_t i=0; i<num-1; i++) ret[i] = static_cast<T>(i) + start;
  ret[num] = stop;

  return;
}


// Returns num evenly spaced samples, calculated over the interval [`start`, `stop`].
// based on numpy linspace https://github.com/numpy/numpy/blob/v1.25.0/numpy/core/function_base.py
template <typename T>
DEVCALLABLE void linspace(T start, T stop, int num, ManVec<T>& ret)
{
  //ManVec<T> ret;

  //ManVec<T> ret = arange(
  //    std::static_cast<T>(0), 
  //    std::static_cast<T>(num)
  //    );

  ret.resize(num);

  // arange
  for(size_t i=0; i<num-1; i++) ret[i] = static_cast<T>(i);

  T delta = stop - start;
  T div   = static_cast<T>(num) - static_cast<T>(1.0); // N - 1

  // NOTE: iterating only to n-1; last element is set manualyl to avoid FP rounding efrrors
  for(size_t i=0; i<num-1; i++) ret[i] = ret[i]*delta/div + start;
  ret[num-1] = stop;

  return;
}


template <typename T>
DEVCALLABLE void logspace(T start, T stop, int num, ManVec<T>& ret)
{
  linspace(start, stop, num, ret);
  T base = static_cast<T>(10.0);

  for(size_t i=0; i<num; i++) ret[i] = std::pow(base, ret[i]);

  return;
}

} // end of ns
