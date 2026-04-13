#pragma once

#include "thrust/device_vector.h"
#include "thrust/host_vector.h"
#include "thrust/memory.h"
#include "tyvi/mdgrid.h"

#include "tools/math.h"

#include <array>
#include <cstddef>
#include <print>
#include <span>
#include <stdexcept>
#include <concepts>

namespace toolbox {

/// Hollow 3D grid of vectors of dimension D on the device.
///
/// Stores only the values of type T which are within h cells
/// from the edge in a contiguous manner on a device.
///
/// Data access is through the accessor object from hollow_grid::get_accessor.
/// It contains the grid extents, pointer to the data and member function
/// hollow_grid::accessor::offset(i, j, k, c) which returns the offset
/// to the data corresponding to the grid position (i, j, k)
/// and c:th component of the D dimensional (mathematical) vector.
///
/// Example:
///
///    auto hgrid = hollow_grid<std::vector<int>, 3, 2>(10, 11, 12);
///    const auto acc = hgrid.get_accessor();
///    acc.data[acc.offset(9, 0, 0, 2)] = 42;
///
/// The offset function will return bogus values if it is called with
/// indices that are outside of the grid, inside the hollow region
/// or with component c s.t. c >= D.
template<typename T, std::size_t D, std::size_t h>
class [[nodiscard]] hollow_grid {
  std::array<std::size_t, 3> extents_;
  thrust::device_vector<T> data_ {};

public:
  template<typename U>
  struct accessor {
    U* data;
    std::ptrdiff_t Nx, Ny, Nz;

    template<std::integral I, std::integral J>
    [[nodiscard]] constexpr std::size_t
      offset(const I iu, const I ju, const I ku, const J component) const
    {
      const auto i             = static_cast<std::ptrdiff_t>(iu);
      const auto j             = static_cast<std::ptrdiff_t>(ju);
      const auto k             = static_cast<std::ptrdiff_t>(ku);
      static constexpr auto hh = static_cast<std::ptrdiff_t>(h);

      const auto whole_i_slices = sstd::min(sstd::max(0z, i - hh), Nx - 2z * hh);
      auto skips                = whole_i_slices * (Ny - 2z * hh) * (Nz - 2z * hh);

      if(i >= hh and i < Nx - hh) {
        const auto whole_j_slices = sstd::min(sstd::max(0z, j - hh), Ny - 2z * hh);
        skips += whole_j_slices * (Nz - 2 * hh);

        if(k >= Nz - hh and j >= hh and j < Ny - hh) { skips += Nz - 2z * hh; }
      }

      const auto size = Nx * Ny * Nz - (Nx - 2 * h) * (Ny - 2 * h) * (Nz - 2 * h);

      return component * size +
             static_cast<std::size_t>(i * Ny * Nz + j * Nz + k - skips);
    }

    template<std::integral I, std::integral J>
    [[nodiscard]] constexpr std::size_t
      offset(const std::array<I, 3> idx, const J component) const
    {
      return this->offset(idx[0], idx[1], idx[2], component);
    }
  };

  constexpr hollow_grid() = default;

  explicit constexpr hollow_grid(
    const std::size_t Nx,
    const std::size_t Ny,
    const std::size_t Nz) :
    extents_ { Nx, Ny, Nz },
    data_(D * (Nx * Ny * Nz - (Nx - 2 * h) * (Ny - 2 * h) * (Nz - 2 * h)))
  {
  }

  explicit hollow_grid(const std::array<std::size_t, 3> e) :
    hollow_grid(e[0], e[1], e[2])
  {
  }

  [[nodiscard]] constexpr auto get_accessor()
  {
    return accessor<T> { thrust::raw_pointer_cast(data_.data()),
                         static_cast<std::ptrdiff_t>(extents_[0]),
                         static_cast<std::ptrdiff_t>(extents_[1]),
                         static_cast<std::ptrdiff_t>(extents_[2]) };
  }

  [[nodiscard]] constexpr auto get_accessor() const
  {
    return accessor<const T> { thrust::raw_pointer_cast(data_.data()),
                               static_cast<std::ptrdiff_t>(extents_[0]),
                               static_cast<std::ptrdiff_t>(extents_[1]),
                               static_cast<std::ptrdiff_t>(extents_[2]) };
  }

  /// Copy values from compatible tyvi::mdgrid::mds.
  ///
  /// Throws if the mds is not compatible.
  constexpr void set_from_mds(const tyvi::mdgrid_work& w, const auto& mds)
  {
    const auto e = mds.extents();
    if(
      e.extent(0) != extents_[0] or e.extent(1) != extents_[1] or
      e.extent(2) != extents_[2]) {
      throw std::runtime_error { std::format(
        "Can not set hollow grid from incompatible mds: {} {} {} != {} {} {}",
        e.extent(0),
        e.extent(1),
        e.extent(2),
        extents_[0],
        extents_[1],
        extents_[2]) };
    }

    auto acc = this->get_accessor();

    w.for_each_index(mds, [=](const auto idx, const auto tidx) {
      const auto i             = static_cast<std::ptrdiff_t>(idx[0]);
      const auto j             = static_cast<std::ptrdiff_t>(idx[1]);
      const auto k             = static_cast<std::ptrdiff_t>(idx[2]);
      static constexpr auto hh = static_cast<std::ptrdiff_t>(h);

      const auto A_i = not(i < hh);
      const auto B_i = not(i >= acc.Nx - hh);

      const auto A_j = not(j < hh);
      const auto B_j = not(j >= acc.Ny - hh);

      const auto A_k = not(k < hh);
      const auto B_k = not(k >= acc.Nz - hh);
      if(A_i and B_i and A_j and B_j and A_k and B_k) { return; }

      acc.data[acc.offset(idx, tidx[0])] = mds[idx][tidx];
    });
  }

  /// Dumps all elements in to stdout.
  void debug_print() const
  {
    const auto buff = thrust::host_vector<float>(data_.begin(), data_.end());

    for(const float x: buff) { std::print("{} ", x); }
    std::println();
  }

  [[nodiscard]] constexpr auto span() &
  {
    return std::span<T>(thrust::raw_pointer_cast(data_.data()), data_.size());
  }
  [[nodiscard]] constexpr auto span() const&
  {
    return std::span<const T>(thrust::raw_pointer_cast(data_.data()), data_.size());
  }
};


}  // namespace toolbox
