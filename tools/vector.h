#pragma once

#include <functional>
// tyvi has missing include of <functional> :(

#include "tyvi/mdspan.h"

#include <array>
#include <cmath>
#include <concepts>
#include <iostream>
#include <type_traits>
#include <utility>

namespace toolbox {

template<std::floating_point T, std::size_t D>
struct VecD {
  std::array<T, D> data {};

  template<typename Self>
  constexpr auto& operator()(this Self& self, const std::size_t ind)
  {
    return self.data[ind];
  }

  template<typename Self>
  constexpr auto& operator[](this Self& self, const std::size_t ind)
  {
    return self.data[ind];
  }

  constexpr VecD() = default;

  template<std::convertible_to<T>... U>
    requires(sizeof...(U) == D)
  constexpr VecD(U&&... args)
  {
    data = std::array<T, D> { static_cast<T>(std::forward<U>(args))... };
  }

  constexpr VecD(const std::array<T, D>& arr) { data = arr; }

  template<typename MDS>
    requires std::convertible_to<T, typename MDS::element_type> and
             (MDS::rank() == 1) and (MDS::static_extent(0) == D)
  constexpr VecD(const MDS& mds)
  {
    for(auto i = 0uz; i < D; ++i) { data[i] = mds[i]; }
  }

  template<std::convertible_to<T>... U>
    requires(sizeof...(U) == D)
  constexpr void set(U&&... args)
  {
    *this = VecD(std::forward<U>(args)...);
  }

  friend std::ostream& operator<<(std::ostream& os, const VecD<T, D>& v)
  {
    os << "[";
    for(size_t i = 0; i < D; i++) os << v(i) << " ";
    os << "]";
    return os;
  }
};


template<std::floating_point T>
using Vec3 = VecD<T, 3>;

template<std::floating_point T>
using Vec4 = VecD<T, 4>;

/// Deduction guide to make `Vec3 v(1, 2, 3);` work.
template<typename... U>
VecD(U&&...) -> VecD<std::common_type_t<U...>, sizeof...(U)>;

template<std::floating_point T, std::size_t D>
  requires(D >= 2)
struct MatD {
  std::array<T, D * D> data {};

  /// Column major order.
  constexpr auto mds(this auto& self)
  {
    static constexpr bool is_const =
      std::is_const_v<std::remove_reference_t<decltype(self)>>;
    using TT = std::conditional_t<is_const, const T, T>;
    return tyvi::sstd::geometric_mdspan<TT, 2, D, std::layout_left>(self.data.data());
  }

  /// Column major order.
  constexpr auto& operator[](this auto& self, const size_t i, const size_t j)
  {
    return self.mds()[i, j];
  }

  /// Column major order.
  constexpr auto& operator()(this auto& self, const size_t i, const size_t j)
  {
    return self.mds()[i, j];
  }

  constexpr MatD() = default;

  template<std::same_as<VecD<T, D>>... U>
    requires(sizeof...(U) == D)
  constexpr MatD(const U&... x)
  {
    const auto vecs = std::array { x... };
    const auto m    = this->mds();
    for(const auto [i, j]: tyvi::sstd::geometric_index_space<2, D>()) {
      m[i, j] = vecs[i](j);
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const MatD<T, D>& M)
  {
    os << std::endl;
    os << "[";
    for(auto i = 0uz; i < D; ++i) {
      os << "[ ";
      for(auto j = 0uz; i < D - 1; ++i) { os << M(i, j) << " , "; }
      os << M(i, D - 1) << "]" << "\n";
    }
    os << "]" << std::endl;

    return os;
  }
};

template<std::floating_point T>
using Mat3 = MatD<T, 3>;

template<std::floating_point T>
using Mat4 = MatD<T, 4>;

//--------------------------------------------------
// arithmetics
// elementwise arithmetics

template<std::floating_point T, std::size_t D>
constexpr VecD<T, D>
  operator/(const VecD<T, D>& v, const T s)
{
  VecD<T, D> ret;
  for(size_t i = 0; i < D; i++) ret(i) = v(i) / s;
  return ret;
}

template<std::floating_point T, std::size_t D>
constexpr VecD<T, D>
  operator*(const VecD<T, D>& v, const T s)
{
  VecD<T, D> ret;
  for(size_t i = 0; i < D; i++) ret(i) = v(i) * s;
  return ret;
}

template<std::floating_point T, std::size_t D>
constexpr VecD<T, D>
  operator*(const T s, const VecD<T, D>& v)
{
  return v * s;
}


template<std::floating_point T, std::size_t D>
constexpr VecD<T, D>
  operator+(const VecD<T, D>& v1, const VecD<T, D>& v2)
{
  VecD<T, D> ret;
  for(size_t i = 0; i < D; i++) ret(i) = v1(i) + v2(i);
  return ret;
}


template<std::floating_point T, std::size_t D>
constexpr VecD<T, D>
  operator-(const VecD<T, D>& v1, const VecD<T, D>& v2)
{
  VecD<T, D> ret;
  for(size_t i = 0; i < D; i++) ret(i) = v1(i) - v2(i);
  return ret;
}

template<std::floating_point T, std::size_t D>
constexpr T
  dot(const VecD<T, D>& v1, const VecD<T, D>& v2)
{
  T ret { 0 };
  for(size_t i = 0; i < D; i++) ret += v1(i) * v2(i);
  return ret;
}

template<std::floating_point T, std::size_t D>
constexpr T
  sum(const VecD<T, D>& v)
{
  // Ugly syntax until C++26 :(
  return [&]<std::size_t... I>(std::index_sequence<I...>) {
    return (v(I) + ...);
  }(std::make_index_sequence<D>());
}

// matrix-vector product
template<std::floating_point T, std::size_t D>
constexpr VecD<T, D>
  dot(const MatD<T, D>& M, const VecD<T, D>& v)
{
  VecD<T, D> ret;

  for(const auto [i, j]: tyvi::sstd::geometric_index_space<2, D>()) {
    ret[i] += M[i, j] * v[j];
  }
  return ret;
}


template<std::floating_point T>
constexpr Vec3<T>
  cross(const Vec3<T>& v1, const Vec3<T>& v2)
{
  Vec3<T> ret;
  ret(0) = v1(1) * v2(2) - v1(2) * v2(1);
  ret(1) = -v1(0) * v2(2) + v1(2) * v2(0);
  ret(2) = v1(0) * v2(1) - v1(1) * v2(0);
  return ret;
}

// unit vector into direction of v1 x v2
template<class T>
inline Vec3<T>
  unit_cross(Vec3<T>& v1, Vec3<T>& v2)
{
  const Vec3<T> ret = cross(v1, v2);
  // TODO prevent dividing by zero
  return ret / norm(ret);
}


// matrix inversion
template<std::floating_point T>
constexpr Mat3<T>
  inv(const Mat3<T>& M)
{
  Mat3<T> ret;

  T det = 0.0;
  for(int i = 0; i < 3; i++)
    det +=
      (M(0, i) *
       (M(1, (i + 1) % 3) * M(2, (i + 2) % 3) - M(1, (i + 2) % 3) * M(2, (i + 1) % 3)));

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ret(i, j) = ((M((j + 1) % 3, (i + 1) % 3) * M((j + 2) % 3, (i + 2) % 3)) -
                   (M((j + 1) % 3, (i + 2) % 3) * M((j + 2) % 3, (i + 1) % 3))) /
                  det;
    }
  }

  return ret;
}

//--------------------------------------------------
// vector norm
template<std::floating_point... T>
constexpr std::common_type_t<T...>
  norm(const T... args)
{
  return std::sqrt(((args * args) + ...));
}

template<std::floating_point T, std::size_t D>
constexpr T
  norm(const VecD<T, D>& v)
{
  return std::sqrt(dot(v, v));
}

template<std::floating_point T, std::size_t D>
constexpr T
  norm2d(const VecD<T, D>& v)  // contracted norm in x y direction
{
  return std::sqrt(v(0) * v(0) + v(1) * v(1));
}

template<std::floating_point T, std::size_t D>
constexpr T
  norm1d(const VecD<T, D>& v)  // contracted norm in x direction
{
  return std::abs(v(0));
}

template<std::floating_point T, std::size_t D>
constexpr T
  norm_minkowski(VecD<T, D>& v)
{
  // NOTE: minkowski norm with -0 + 1 + 2 + 3 + ... signature. However, we take absolute
  // value of the output to be agnostic about the metric signature.

  // Ugly syntax until C++26 :(
  const auto spatial_part = [&]<std::size_t... I>(std::index_sequence<I...>) {
    return ((v(I + 1) * v(I + 1)) + ...);
  }(std::make_index_sequence<D - 1>());

  return std::sqrt(std::abs(-v(0) * v(0) + spatial_part));
}


}  // namespace toolbox
