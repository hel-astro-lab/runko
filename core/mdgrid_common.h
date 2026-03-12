#pragma once

#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <cstddef>

namespace runko {

using index_t = std::uint32_t;

namespace detail {

template<typename T>
static constexpr auto scalar_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 0, .dim = 3 };

template<typename T>
static constexpr auto vec_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 1, .dim = 3 };

using list_extents = std::dextents<index_t, 1>;
using grid_extents = std::dextents<index_t, 3>;

template<typename T>
static constexpr auto io_field_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 1, .dim = 14 };

template<typename T>
static constexpr auto prtcl_io_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 1, .dim = 12 };

template<typename T>
static constexpr auto spectra_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 1,
                                       .dim  = 12 };  // 4 spectra * max 3 species

using list_extents   = std::dextents<runko::index_t, 1>;
using grid2d_extents = std::dextents<runko::index_t, 2>;
using grid_extents   = std::dextents<runko::index_t, 3>;

}  // namespace detail

template<typename T>
using ScalarGrid = tyvi::mdgrid<detail::scalar_element<T>, detail::grid_extents>;

template<typename T>
using VecGrid = tyvi::mdgrid<detail::vec_element<T>, detail::grid_extents>;

template<typename T>
using ScalarList = tyvi::mdgrid<detail::scalar_element<T>, detail::list_extents>;

template<typename T>
using VecList = tyvi::mdgrid<detail::vec_element<T>, detail::list_extents>;

template<typename T>
using IOFieldGrid = tyvi::mdgrid<detail::io_field_element<T>, detail::grid_extents>;

template<typename T>
using PrtclFieldList = tyvi::mdgrid<detail::prtcl_io_element<T>, detail::list_extents>;

template<typename T>
using SpectraGrid = tyvi::mdgrid<detail::spectra_element<T>, detail::grid2d_extents>;

}  // namespace runko
