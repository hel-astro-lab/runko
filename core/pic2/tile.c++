#include "core/pic2/tile.h"

#include "core/emf2/yee_lattice.h"
#include "core/particles_common.h"
#include "core/pic2/particle.h"

#include <cmath>
#include <format>
#include <functional>
#include <numeric>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

void
  construct_particle_buffs(auto& where_to, const toolbox::ConfigParser& conf)
{
  using opt_args_t = std::optional<pic2::ParticleContainerArgs>;

  // Return optional containing the arguments if q and m are found in conf.
  const auto make_opt_args =
    [&](const std::string& q_label, const std::string& m_label) -> opt_args_t {
    const auto q_to_qm_tuple = [&](const auto q) {
      return conf.get<double>(m_label).transform(
        [&](const auto m) { return std::tuple { q, m }; });
    };

    const auto qm_tuple_to_args = [&](const auto qm) {
      const auto [q, m] = qm;

      // Due to effects of density function to particle number,
      // we do not yet know how many particles there are going to be.
      return pic2::ParticleContainerArgs { .N = 0, .charge = q, .mass = m };
    };

    return conf.get<double>(q_label)
      .and_then(q_to_qm_tuple)
      .transform(qm_tuple_to_args);
  };

  for(auto i = 0uz; true; ++i) {
    const auto q_label = std::format("q{}", i);
    const auto m_label = std::format("m{}", i);

    if(const auto pcontainer_args = make_opt_args(q_label, m_label)) {
      // insert_or_assign and not operator[], because if element is missing,
      // then operator[] will default construct it.
      // pic2::ParticleContainer is not default constructible.
      std::ignore = where_to.insert_or_assign(
        i,
        pic2::ParticleContainer { pcontainer_args.value() });

    } else {
      break;
    }
  }
}

pic2::ParticlePusher
  parse_particle_pusher(const std::string_view p)
{
  if(p == "boris") {
    return pic2::ParticlePusher::boris;
  } else {
    const auto msg = std::format("{} is not supported particle pusher.", p);
    throw std::runtime_error { msg };
  }
}

pic2::FieldInterpolator
  parse_field_interpolator(const std::string_view p)
{
  if(p == "linear_1st") {
    return pic2::FieldInterpolator::linear_1st;
  } else {
    const auto msg = std::format("{} is not supported field_interpolator.", p);
    throw std::runtime_error { msg };
  }
}

pic2::CurrentDepositer
  parse_current_depositer(const std::string_view p)
{
  if(p == "zigzag" or p == "zigzag_1st") {
    return pic2::CurrentDepositer::zigzag_1st;
  } else if(p == "zigzag_1st_atomic") {
    return pic2::CurrentDepositer::zigzag_1st_atomic;
  } else {
    const auto msg = std::format("{} is not supported current depositer.", p);
    throw std::runtime_error { msg };
  }
}

}  // namespace

namespace pic2 {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_grid_idx,
  const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf2::Tile<D>(tile_grid_idx, conf),
  particle_pusher_ { parse_particle_pusher(
    conf.get_or_throw<std::string>("particle_pusher")) },
  field_interpolator_ { parse_field_interpolator(
    conf.get_or_throw<std::string>("field_interpolator")) },
  current_depositer_ { parse_current_depositer(
    conf.get_or_throw<std::string>("current_depositer")) }
{
  construct_particle_buffs(particle_buffs_, conf);

  const auto xmin   = conf.get_or_throw<double>("xmin");
  const auto ymin   = conf.get_or_throw<double>("ymin");
  const auto zmin   = conf.get_or_throw<double>("zmin");
  const auto Nx     = conf.get_or_throw<double>("Nx");
  const auto Ny     = conf.get_or_throw<double>("Ny");
  const auto Nz     = conf.get_or_throw<double>("Nz");
  const auto NxMesh = conf.get_or_throw<double>("NxMesh");
  const auto NyMesh = conf.get_or_throw<double>("NyMesh");
  const auto NzMesh = conf.get_or_throw<double>("NzMesh");

  using PVT                     = ParticleContainer::value_type;
  this->global_coordinate_mins_ = { static_cast<PVT>(xmin),
                                    static_cast<PVT>(ymin),
                                    static_cast<PVT>(zmin) };

  this->global_coordinate_maxs_ = {
    global_coordinate_mins_[0] + static_cast<PVT>(Nx * NxMesh),
    global_coordinate_mins_[1] + static_cast<PVT>(Ny * NyMesh),
    global_coordinate_mins_[2] + static_cast<PVT>(Nz * NzMesh)
  };
}

template<std::size_t D>
std::array<std::vector<typename Tile<D>::value_type>, 3>
  Tile<D>::get_positions(const std::size_t p)
{
  return particle_buffs_.at(p).get_positions();
}

template<std::size_t D>
std::array<std::vector<typename Tile<D>::value_type>, 3>
  Tile<D>::get_velocities(const std::size_t p)
{
  return particle_buffs_.at(p).get_velocities();
}

template<std::size_t D>
void
  Tile<D>::inject_to_each_cell(const std::size_t particle_type, particle_generator pgen)
{
  if(this->mins == this->maxs) {
    throw std::logic_error {
      "failed pic2::Tile::inject_to_each_cell precondition: corgi::tile::{mins,maxs} "
      "aren't set"
    };
  }

  /// FIXME: unify global coordinates from here and emf2::Tile::set_EBJ to mixin.
  const auto Lx = static_cast<double>(this->maxs[0] - this->mins[0]);
  const auto Ly = static_cast<double>(this->maxs[1] - this->mins[1]);
  const auto Lz = static_cast<double>(this->maxs[2] - this->mins[2]);

  const auto e = this->extents_wout_halo();

  auto global_coordinates = [&](const auto i, const auto j, const auto k) {
    const auto x_coeff = static_cast<double>(i) / static_cast<double>(e[0]);
    const auto y_coeff = static_cast<double>(j) / static_cast<double>(e[1]);
    const auto z_coeff = static_cast<double>(k) / static_cast<double>(e[2]);

    return std::tuple<double, double, double> {
      static_cast<double>(this->mins[0]) + x_coeff * Lx,
      static_cast<double>(this->mins[1]) + y_coeff * Ly,
      static_cast<double>(this->mins[2]) + z_coeff * Lz
    };
  };

  std::vector<runko::ParticleState> new_particles {};

  namespace rv     = std::views;
  auto index_space = rv::cartesian_product(
    rv::iota(0uz, e[0]),
    rv::iota(0uz, e[1]),
    rv::iota(0uz, e[2]));
  for(const auto [i, j, k]: index_space) {
    const auto [x, y, z] = global_coordinates(i, j, k);
    for(const auto p: pgen(x, y, z)) { new_particles.push_back(p); }
  }

  particle_buffs_.at(particle_type).add_particles(new_particles);
}

template<std::size_t D>
void
  Tile<D>::inject(
    const std::size_t particle_type,
    const std::vector<runko::ParticleState>& new_particles)
{
  particle_buffs_.at(particle_type).add_particles(new_particles);
}

template<std::size_t D>
void
  Tile<D>::batch_inject_to_cells(
    const std::size_t particle_type,
    batch_particle_generator pgen)
{

  /// FIXME: unify global coordinates
  const auto Lx = static_cast<double>(this->maxs[0] - this->mins[0]);
  const auto Ly = static_cast<double>(this->maxs[1] - this->mins[1]);
  const auto Lz = static_cast<double>(this->maxs[2] - this->mins[2]);

  const auto e = this->yee_lattice_.extents_wout_halo();

  auto global_coordinates = [&](const auto i, const auto j, const auto k) {
    const auto x_coeff = static_cast<double>(i) / static_cast<double>(e[0]);
    const auto y_coeff = static_cast<double>(j) / static_cast<double>(e[1]);
    const auto z_coeff = static_cast<double>(k) / static_cast<double>(e[2]);

    return std::array { static_cast<double>(this->mins[0]) + x_coeff * Lx,
                        static_cast<double>(this->mins[1]) + y_coeff * Ly,
                        static_cast<double>(this->mins[2]) + z_coeff * Lz };
  };

  const auto cells_in_total = std::array { e[0] * e[1] * e[2] };
  auto x                    = pybind11::array_t<double>(cells_in_total);
  auto y                    = pybind11::array_t<double>(cells_in_total);
  auto z                    = pybind11::array_t<double>(cells_in_total);

  {
    auto xv = x.template mutable_unchecked<1>();
    auto yv = y.template mutable_unchecked<1>();
    auto zv = z.template mutable_unchecked<1>();

    // Use this and tyvi index_space to easily iterate over 3D indices.
    const auto helper_mds = std::mdspan((int*)nullptr, e[0], e[1], e[2]);
    // Use the mapping directly to get 1D offset from the 3D indices.
    const auto m = helper_mds.mapping();

    for(const auto [i, j, k]: tyvi::sstd::index_space(helper_mds)) {
      const auto n                  = m(i, j, k);
      const auto gc = global_coordinates(i, j, k);
      xv(n) = gc[0];
      yv(n) = gc[1];
      zv(n) = gc[2];
    }
  }

  const auto state_batch = pgen(x, y, z);
  const auto batch_size  = state_batch.pos[0].shape(0);

  {
    const auto assert_shapes = [&](const auto& a) {
      if(a.ndim() != 1) {
        throw std::runtime_error {
          "pic2::Tile::batch_inject_to_cells: given batch must be one dimensional."
        };
      }

      if(a.shape(0) != batch_size) {
        throw std::runtime_error {
          "pic2::Tile::batch_inject_to_cells: batches must have same length."
        };
      }
    };

    assert_shapes(state_batch.pos[0]);
    assert_shapes(state_batch.pos[1]);
    assert_shapes(state_batch.pos[2]);
    assert_shapes(state_batch.vel[0]);
    assert_shapes(state_batch.vel[1]);
    assert_shapes(state_batch.vel[2]);
  }

  const auto posx_view = state_batch.pos[0].template unchecked<1>();
  const auto posy_view = state_batch.pos[1].template unchecked<1>();
  const auto posz_view = state_batch.pos[2].template unchecked<1>();

  const auto velx_view = state_batch.vel[0].template unchecked<1>();
  const auto vely_view = state_batch.vel[1].template unchecked<1>();
  const auto velz_view = state_batch.vel[2].template unchecked<1>();

  // This is bit of a waste to go for SOA to AOS and then in inject back to SOA.
  // However, I don't think this matters.
  auto states   = std::vector<runko::ParticleState>(batch_size);
  using index_t = std::remove_cvref_t<decltype(batch_size)>;
  for(const auto n: std::views::iota(index_t { 0 }, batch_size)) {
    const auto posx = posx_view(n);
    const auto posy = posy_view(n);
    const auto posz = posz_view(n);

    const auto velx = velx_view(n);
    const auto vely = vely_view(n);
    const auto velz = velz_view(n);

    states[n] =
      runko::ParticleState { .pos = { posx, posy, posz }, .vel = { velx, vely, velz } };
  }

  this->inject(particle_type, states);
}

template<std::size_t D>
void
  Tile<D>::push_particles()
{
  using yee_value_type = emf2::YeeLattice::value_type;
  const auto origo_pos =
    std::array { static_cast<yee_value_type>(this->mins[0]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[1]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[2]) - this->halo_size };

  pic2::ParticleContainer::InterpolatedEB_function ipol_func {};

  switch(field_interpolator_) {
    case FieldInterpolator::linear_1st: {
      ipol_func = std::bind_front(
        &emf2::YeeLattice::interpolate_EB_linear_1st,
        std::cref(this->yee_lattice_),
        origo_pos);
      break;
    }
    default:
      throw std::logic_error { "pic2::Tile::push_particles: unkown interpolator" };
  }


  switch(particle_pusher_) {
    case ParticlePusher::boris:
      for(auto& [_, pbuff]: particle_buffs_) {
        pbuff.push_particles_boris(this->cfl_, ipol_func);
      }
      break;
    default:
      throw std::logic_error { "pic2::Tile::push_particles: unkown particle pusher" };
  }
}

template<std::size_t D>
void
  Tile<D>::deposit_current()
{
  this->yee_lattice_.clear_current();

  using yee_value_type = emf2::YeeLattice::value_type;
  const auto origo_pos =
    std::array { static_cast<yee_value_type>(this->mins[0]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[1]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[2]) - this->halo_size };

  switch(current_depositer_) {
    case CurrentDepositer::zigzag_1st:
      for(const auto& [_, pcontainer]: this->particle_buffs_) {
        this->yee_lattice_.deposit_current(
          pcontainer.current_zigzag_1st(origo_pos, this->cfl_));
      }
      break;
    case CurrentDepositer::zigzag_1st_atomic: {
      auto generated_J = runko::VecGrid<emf2::YeeLattice::value_type>(
        this->yee_lattice_.extents_with_halo());

      // This might be unneccesseary but I don't trust
      // underlying thrust::device vector to zero initialize the data.
      const auto genJmds = generated_J.mds();
      tyvi::mdgrid_work {}
        .for_each_index(
          genJmds,
          [=](const auto idx, const auto tidx) { genJmds[idx][tidx] = 0; })
        .wait();

      for(const auto& [_, pcontainer]: this->particle_buffs_) {
        pcontainer.current_zigzag_1st(generated_J, origo_pos, this->cfl_);
      }

      this->yee_lattice_.deposit_current(generated_J);
      break;
    }
    default:
      throw std::logic_error {
        "pic2::Tile::deposit_current: unkown current depositer"
      };
  }
}

template<std::size_t D>
void
  Tile<D>::sort_particles()
{
  const auto m = this->yee_lattice_.grid_mapping_with_halo();
  using M      = decltype(m);

  using F              = pic2::ParticleContainer::value_type;
  const auto origo_pos = std::array { static_cast<F>(this->mins[0]) - this->halo_size,
                                      static_cast<F>(this->mins[1]) - this->halo_size,
                                      static_cast<F>(this->mins[2]) - this->halo_size };
  using Vec3F          = toolbox::Vec3<F>;

  auto score = [=](const F x, const F y, const F z) {
    const auto dx  = Vec3F(x, y, z) - Vec3F(origo_pos);
    const auto idx = dx.template as<typename M::index_type>();

    return m(idx[0], idx[1], idx[2]);
  };

  for(auto& [_, pbuff]: this->particle_buffs_) { pbuff.sort(std::move(score)); }
}

}  // namespace pic2

template class pic2::Tile<3>;
