#pragma once

#include "core/emf/yee_lattice.h"
#include "corgi/corgi.h"
#include "corgi/tile.h"
#include "mpi4cpp/mpi.h"
#include "pybind11/numpy.h"
#include "tools/config_parser.h"

#include <array>
#include <cstddef>
#include <experimental/mdspan>
#include <functional>
#include <optional>
#include <tuple>

namespace emf {

enum class FieldPropagator { FDTD2 };
enum class CurrentFilter { binomial2 };

/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in
 * staggered Yee lattice.
 */
template<std::size_t D>
class Tile : virtual public corgi::Tile<D> {

  static_assert(D == 3);

protected:
  YeeLattice yee_lattice_;
  double cfl_;
  FieldPropagator field_propagator_;
  std::optional<CurrentFilter> current_filter_ {};

public:
  static constexpr auto halo_size = 3;

  /// Construct Tile by deducing extents from given tile grid index and config.
  ///
  /// Config has to contain values for:
  ///
  /// `cfl`: (FIXME: add description)
  /// `N{x,y,z}`: tile grid extents
  /// `N{x,y,z}Mesh`: extents of mesh/grid in each tile
  /// `{x,y,z}min: minimum coordinate values
  /// `d{x,y,z}`: coordinate distance between mesh/grid points
  /// `field_propagator`: scheme to propagate E and B fields.
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_indices,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;

  using vector_field_function =
    std::function<std::tuple<double, double, double>(double, double, double)>;

  /// Set E, B and J fields in non-halo regions.
  void
    set_EBJ(vector_field_function E, vector_field_function B, vector_field_function J);

  using batch_array = pybind11::array_t<double>;
  using batch_vector_field_function =
    std::function<batch_array(batch_array, batch_array, batch_array)>;

  /// Set E, B and J fields in non-halo regions.
  void batch_set_EBJ(
    batch_vector_field_function Ex,
    batch_vector_field_function Ey,
    batch_vector_field_function Ez,
    batch_vector_field_function Bx,
    batch_vector_field_function By,
    batch_vector_field_function Bz,
    batch_vector_field_function Jx,
    batch_vector_field_function Jy,
    batch_vector_field_function Jz);

  /// Get E, B and J fields in non-halo regions.
  YeeLattice::YeeLatticeHostCopy get_EBJ();

  /// Get E, B and J fields including halo regions.
  YeeLattice::YeeLatticeHostCopy get_EBJ_with_halo();

  /// Returns mdspans to host accessible E, B and J in non-halo region.
  auto view_EBJ_on_host() { return yee_lattice_.view_EBJ_on_host(); }

  /// Size of the non-halo yee lattice.
  std::array<std::size_t, 3> extents_wout_halo() const;

  /// Advance B by half time step using scheme from configuration in non-halo region.
  void push_half_b();

  /// Advance E by full time step using scheme from configuration in non-halo region.
  ///
  /// Does not add contributions from the current.
  void push_e();

  /// E -= J
  void subtract_J_from_E();

  /// Applies potentially configured filter to J.
  ///
  /// Throws if filter is not configured.
  void filter_current();

  /// Send field specified with mode (see comm_mode).
  std::vector<mpi4cpp::mpi::request>
    send_data(mpi4cpp::mpi::communicator& /*comm*/, int dest, int mode, int tag)
      override;

  /// Receive field specified with mode (see comm_mode).
  std::vector<mpi4cpp::mpi::request>
    recv_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int mode, int tag)
      override;

  /// Get halo region of field specified with mode (see comm_mode) from other.
  ///
  /// Assumes that the other tile is emf::Tile or its descendant.
  void pairwise_moore_communication(
    const corgi::Tile<D>& /* other */,
    const std::array<int, D> dir_to_other,
    const int /* mode */
    ) override;
};

}  // namespace emf
