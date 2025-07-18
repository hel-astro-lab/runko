#pragma once

#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <cstddef>

namespace runko {

namespace detail {

template<typename T>
static constexpr auto scalar_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 0, .dim = 3 };

template<typename T>
static constexpr auto vec_element =
  tyvi::mdgrid_element_descriptor<T> { .rank = 1, .dim = 3 };

using list_extents = std::dextents<std::size_t, 1>;
using grid_extents = std::dextents<std::size_t, 3>;

}  // namespace detail

template<typename T>
using ScalarGrid = tyvi::mdgrid<detail::scalar_element<T>, detail::grid_extents>;

template<typename T>
using VecGrid = tyvi::mdgrid<detail::vec_element<T>, detail::grid_extents>;

template<typename T>
using ScalarList = tyvi::mdgrid<detail::scalar_element<T>, detail::list_extents>;

template<typename T>
using VecList = tyvi::mdgrid<detail::vec_element<T>, detail::list_extents>;

}  // namespace runko
