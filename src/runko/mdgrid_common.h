// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "thrust/copy.h"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"
#include "tyvi/device_allocator.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdsegments.h"
#include "tyvi/mdspan.h"

#include <cstddef>
#include <stdexcept>

namespace runko {

using index_t = std::uint32_t;

namespace detail {

using list_extents   = std::dextents<index_t, 1>;
using grid2d_extents = std::dextents<index_t, 2>;
using grid_extents   = std::dextents<index_t, 3>;

using scalar_element_extents   = std::extents<index_t>;
using vec_element_extents      = std::extents<index_t, 3>;
using io_field_element_extents = std::extents<index_t, 14>;
using prtcl_io_element_extents = std::extents<index_t, 12>;
using spectra_element_extents = std::extents<index_t, 12>;  // 4 spectra * max 3 species


using element_layout = std::layout_right;
using grid_layout    = std::layout_right;

template<typename T, typename ElementExtents, typename GridExtents>
using device_grid = tyvi::mdgrid_buffer<
  thrust::device_vector<T>,
  ElementExtents,
  element_layout,
  GridExtents,
  grid_layout>;

template<typename DeviceGrid>
struct host_mirror;

template<typename T, typename EE, typename EL, typename GE, typename GL>
struct host_mirror<tyvi::mdgrid_buffer<thrust::device_vector<T>, EE, EL, GE, GL>> {
  using type = tyvi::mdgrid_buffer<thrust::host_vector<T>, EE, EL, GE, GL>;
};

}  // namespace detail

template<typename T>
using ScalarGrid =
  detail::device_grid<T, detail::scalar_element_extents, detail::grid_extents>;

template<typename T>
using VecGrid =
  detail::device_grid<T, detail::vec_element_extents, detail::grid_extents>;

template<typename T>
using ScalarList =
  detail::device_grid<T, detail::scalar_element_extents, detail::list_extents>;

template<typename T>
using VecList =
  detail::device_grid<T, detail::vec_element_extents, detail::list_extents>;

template<typename T>
using IOFieldGrid =
  detail::device_grid<T, detail::io_field_element_extents, detail::grid_extents>;

template<typename T>
using PrtclFieldList =
  detail::device_grid<T, detail::prtcl_io_element_extents, detail::list_extents>;

template<typename T>
using SpectraGrid =
  detail::device_grid<T, detail::spectra_element_extents, detail::grid2d_extents>;

/// Host memory grid with the same shape and layout as the given device grid.
template<typename DeviceGrid>
using HostMirror = typename detail::host_mirror<DeviceGrid>::type;

/// Copies data from one grid to another using thrust::copy.
///
/// Using thrust::copy means that it will automaticlly dispatch to
/// correct copy mechanism (h2d, d2h, etc.) as tyvi::mdgrid_buffer::begin
/// and *::end gives iterators to underlying vector which thrust can detect.
template<typename GridA, typename GridB>
void
  thrust_copy(const tyvi::mdgrid_work& w, const GridA& a, GridB& b)
{
  if(a.grid_extents() != b.grid_extents()) {
    throw std::invalid_argument { "thrust_copy: extents of given grids differ." };
  }
  thrust::copy(w.on_this(), a.begin(), a.end(), b.begin());
}

static constexpr auto segment_size = 1uz << 20uz;

template<typename T>
using device_scalar_segments = tyvi::mdsegments<
  T,
  segment_size,
  std::extents<index_t>,
  std::layout_right,
  tyvi::device_allocator<T>>;

template<typename T>
using host_scalar_segments = tyvi::mdsegments<T, segment_size, std::extents<index_t>>;

template<typename T>
using device_vec_segments = tyvi::mdsegments<
  T,
  segment_size,
  std::extents<index_t, 3uz>,
  std::layout_right,
  tyvi::device_allocator<T>>;

template<typename T>
using host_vec_segments =
  tyvi::mdsegments<T, segment_size, std::extents<index_t, 3uz>, std::layout_right>;


}  // namespace runko
