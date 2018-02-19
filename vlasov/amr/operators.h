#pragma once

#include <cmath> 
#include <algorithm>

#include "mesh.h"


namespace toolbox {


// Memory optimized (move semantic) binary operators


// addition +
//--------------------------------------------------

template<typename T,int D>
inline AdaptiveMesh<T,D> operator +(
    AdaptiveMesh<T,D>& lhs,
    AdaptiveMesh<T,D>& rhs)
{
  lhs += rhs;

  return lhs;
}

  
template<typename T,int D>
inline AdaptiveMesh<T,D> operator +(
    const AdaptiveMesh<T,D>& lhs,
    const AdaptiveMesh<T,D>& rhs)
{
  AdaptiveMesh<T,D> ret(lhs);
  ret += rhs;

  return ret;
}


template<typename T, int D>
inline AdaptiveMesh<T,D>& operator +(
          AdaptiveMesh<T,D>& lhs,
    const AdaptiveMesh<T,D>& rhs)
{
  lhs += rhs;
  return lhs;
}

template<typename T, int D>
inline AdaptiveMesh<T,D>& operator +(
    const AdaptiveMesh<T,D>& lhs,
          AdaptiveMesh<T,D>& rhs)
{
  rhs += lhs;
  return std::move(rhs);
}



// subtraction -
//--------------------------------------------------
template<typename T,int D>
inline AdaptiveMesh<T,D> operator -(
    AdaptiveMesh<T,D>& lhs,
    AdaptiveMesh<T,D>& rhs)
{
  lhs -= rhs;

  return lhs;
}

template<typename T,int D>
inline AdaptiveMesh<T,D> operator -(
    const AdaptiveMesh<T,D>& lhs,
    const AdaptiveMesh<T,D>& rhs)
{
  AdaptiveMesh<T,D> ret(lhs);
  ret -= rhs;

  return ret;
}

template<typename T, int D>
inline AdaptiveMesh<T,D>& operator -(
          AdaptiveMesh<T,D>& lhs,
    const AdaptiveMesh<T,D>& rhs)
{
  lhs -= rhs;
  return lhs;
}

template<typename T, int D>
inline AdaptiveMesh<T,D>& operator -(
    const AdaptiveMesh<T,D>& lhs,
          AdaptiveMesh<T,D>& rhs)
{
  rhs -= lhs;
  return std::move(rhs);
}





}
