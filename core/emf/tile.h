#pragma once

#include "core/emf/antenna.h"
#include "core/emf/edge_bc.h"
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
#include <map>
#include <optional>
#include <tuple>
#include <vector>

#include <mpi.h>

namespace emf {

enum class FieldPropagator { FDTD2 };
enum class CurrentFilter { binomial2, binomial2_unrolled };

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

  /// lap_coeffs are used in reverse order (see Tile::register_antenna implementation).
  std::vector<emf::antenna_mode> antenna_modes_ {};

  /// Registered edge boundary conditions, applied in order.
  std::vector<emf::edge_bc> edge_bcs_ {};

  /// Compute the number of cells this tile is covered by the given edge BC.
  ///
  /// Returns nullopt if the edge does not overlap this tile.
  std::optional<std::size_t> edge_bc_width(const edge_bc& bc) const;

  std::array<double, 3> global_coordinate_mins_;
  std::array<double, 3> global_coordinate_maxs_;

  struct PersistentRequestKey {
    int dest;
    int mode;
    int tag;
    bool operator<(const PersistentRequestKey& o) const {
      if (dest != o.dest) return dest < o.dest;
      if (mode != o.mode) return mode < o.mode;
      return tag < o.tag;
    }
  };

  std::map<PersistentRequestKey, MPI_Request> persistent_send_requests_;
  std::map<PersistentRequestKey, MPI_Request> persistent_recv_requests_;
  bool persistent_requests_initialized_ = false;

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
  ~Tile();

  void initialize_persistent_requests(
    mpi4cpp::mpi::communicator& comm,
    const std::vector<int>& dest_ranks,
    const std::vector<int>& orig_ranks,
    int mode) override;
  void cleanup_persistent_requests() override;

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

  /// Returns device mdspans to non-halo E, B and J without device→host copy.
  auto view_EBJ_on_device() { return yee_lattice_.view_EBJ_on_device(); }

  /// Size of the non-halo yee lattice.
  std::array<std::size_t, 3> extents_wout_halo() const;

  /// Advance B by half time step using scheme from configuration in non-halo region.
  void push_half_b();

  /// Advance E by full time step using scheme from configuration in non-halo region.
  ///
  /// Does not add contributions from the current.
  void push_e();

  /// E -= J
  void add_current();

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
  void local_communication(
    const corgi::Tile<D>& /* other */,
    const std::array<int, D> dir_to_other,
    const int /* mode */
    ) override;


  /// Returns a lambda that maps indices (i, j, k) of this tile to the global
  /// coordinates.
  ///
  /// global_coordinates()(0, 0, 0) corresponds to point this->mins.
  ///
  /// The lambda will return the coordinates as std::array<double, 3>
  /// and it can be invoked even with fractional indices,
  /// i.e. global_coordinates()(0.5, 0.5, 0.5) maps to middle of the cell (0, 0, 0).
  auto global_coordinate_map() const;

  template<typename T>
  std::array<T, 3> get_global_coordinate_mins() const;

  template<typename T>
  std::array<T, 3> get_global_coordinate_maxs() const;

  void register_antenna(emf::antenna_mode);
  void deposit_antenna_current();

  /// Register an edge boundary condition on this tile.
  void register_edge_bc(edge_bc bc);

  /// Apply all registered edge BCs for the given field mode.
  ///
  /// Mode values: comm_mode::emf_E, comm_mode::emf_B, comm_mode::emf_J.
  void apply_edge_bcs(int mode);

  /// Apply a single edge BC for the given field mode.
  ///
  /// Mode values: comm_mode::emf_E, comm_mode::emf_B, comm_mode::emf_J.
  void apply_edge_bc(const edge_bc& bc, int mode);

  /// Returns the total energy in the magnetic field.
  ///
  /// Energy is given in code units: sum_{ijk} B_{ijk}^2 / (8pi)
  ///
  /// Note that (Delta x)^3 = 1, so it does not appear in above expression.
  double total_energy_B() const;

  /// Returns the total energy in the electric field.
  ///
  /// Energy is given in code units: sum_{ijk} E_{ijk}^2 / (8pi)
  ///
  /// Note that (Delta x)^3 = 1, so it does not appear in above expression.
  double total_energy_E() const;
};

template<std::size_t D>
auto
  Tile<D>::global_coordinate_map() const
{

  const auto mins = this->mins;
  const auto maxs = this->maxs;

  const auto Lx = static_cast<double>(maxs[0] - mins[0]);
  const auto Ly = static_cast<double>(maxs[1] - mins[1]);
  const auto Lz = static_cast<double>(maxs[2] - mins[2]);

  const auto e = yee_lattice_.extents_wout_halo();

  return [=](const auto i, const auto j, const auto k) {
    const auto x_coeff = static_cast<double>(i) / static_cast<double>(e[0]);
    const auto y_coeff = static_cast<double>(j) / static_cast<double>(e[1]);
    const auto z_coeff = static_cast<double>(k) / static_cast<double>(e[2]);

    return std::array<double, 3> { static_cast<double>(mins[0]) + x_coeff * Lx,
                                   static_cast<double>(mins[1]) + y_coeff * Ly,
                                   static_cast<double>(mins[2]) + z_coeff * Lz };
  };
}

template<std::size_t D>
template<typename T>
std::array<T, 3>
  Tile<D>::get_global_coordinate_mins() const
{

  return { static_cast<T>(this->global_coordinate_mins_[0]),
           static_cast<T>(this->global_coordinate_mins_[1]),
           static_cast<T>(this->global_coordinate_mins_[2]) };
}

template<std::size_t D>
template<typename T>
std::array<T, 3>
  Tile<D>::get_global_coordinate_maxs() const
{

  return { static_cast<T>(this->global_coordinate_maxs_[0]),
           static_cast<T>(this->global_coordinate_maxs_[1]),
           static_cast<T>(this->global_coordinate_maxs_[2]) };
}


}  // namespace emf
