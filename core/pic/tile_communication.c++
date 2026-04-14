#include "core/communication_common.h"
#include "core/pic/tile.h"
#include "core/pic/virtual_tile.h"
#include "thrust/memory.h"
#include "tools/system.h"

#include <cstddef>
#include <iterator>
#include <ranges>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace pic {

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  Tile<D>::send_data(
    mpi4cpp::mpi::communicator& comm,
    const int dest,
    const int mode,
    const int tag)
{

  // GPU backend requires GPU-aware MPI to pass device pointers directly;
  // CPU backend uses host memory where standard MPI works.
#ifndef TYVI_BACKEND_CPU
  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error { "GPU backend requires GPU-aware MPI." };
  }
#endif

  /*
    All tiles should have the particle containers
    for all particles that are configure in config.
    This means that we can assume that all tiles should be
    in agreement on what particles there are.
   */

  using runko::comm_mode;

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::number_of_particles:
      return std::vector<mpi4cpp::mpi::request> { comm.isend(
        dest,
        tag,
        thrust::raw_pointer_cast(subregion_particle_ends_.data()),
        subregion_particle_ends_.size()) };
    case comm_mode::pic_particle: {
      const auto p = thrust::raw_pointer_cast(this->subregion_particle_buff_.data());
      const auto s = std::span(p, subregion_particle_buff_.size());
      return std::vector<mpi4cpp::mpi::request> { comm.isend(
        dest,
        tag,
        reinterpret_cast<char const*>(std::as_bytes(s).data()),
        s.size_bytes()) };
    }
    default: return emf::Tile<D>::send_data(comm, dest, mode, tag);
  }
}

template<std::size_t D>
void
  Tile<D>::pack_outgoing_particles()
{
  const auto x_div =
    std::array { static_cast<ParticleContainer::value_type>(this->mins[0]),
                 static_cast<ParticleContainer::value_type>(this->maxs[0]) };
  const auto y_div =
    std::array { static_cast<ParticleContainer::value_type>(this->mins[1]),
                 static_cast<ParticleContainer::value_type>(this->maxs[1]) };
  const auto z_div =
    std::array { static_cast<ParticleContainer::value_type>(this->mins[2]),
                 static_cast<ParticleContainer::value_type>(this->maxs[2]) };

  this->subregion_particle_ends_.resize(27 * this->particle_buffs_.size());
  this->subregion_particle_buff_.resize(0);
  // ptypes are assumed always to be contiguous 0, 1, ..., N
  for(auto& [ptype, pcontainer]: this->particle_buffs_) {
    auto spans = pcontainer.divide_to_subregions(
      this->subregion_particle_buff_,
      x_div,
      y_div,
      z_div);


    for(const auto [dir, span]: spans) {
      this->subregion_particle_ends_.at(27 * ptype + dir.neighbor_index()) =
        std::get<1>(span);
    }
  }
}

template<std::size_t D>
void
  Tile<D>::local_communication_postlude(const int mode)
{
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) { return; }

  for(const auto& [ptype, spans]: this->incoming_subregion_particles_) {
    particle_buffs_.at(ptype).append(spans);
  }

  for(auto& [_, pbuff]: this->particle_buffs_) {
    using T = pic::ParticleContainer::value_type;
    pbuff.wrap_positions(
      this->template get_global_coordinate_mins<T>(),
      this->template get_global_coordinate_maxs<T>());
  }

  this->subregion_particle_ends_.clear();
  this->incoming_subregion_particles_.clear();
  // Note that this->subregion_particle_buffs_ can not be cleared here.
  // Other tiles might be using it.
}

template<std::size_t D>
void
  Tile<D>::local_communication(
    const corgi::Tile<D>& other_base,
    const std::array<int, D> dir_to_other,
    const int mode)
{
  auto const* const other_base_ptr = &other_base;
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) {
    emf::Tile<D>::local_communication(other_base, dir_to_other, mode);
    return;
  }

  const auto dir          = runko::grid_neighbor<3>(dir_to_other);
  const auto inverted_dir = dir.inverted();

  namespace rv = std::views;

  auto get_span =
    [&](
      const std::size_t ptype,
      const runko::grid_neighbor<3> dir,
      const std::vector<std::size_t>& ends,
      const thrust::device_vector<runko::ParticleState<value_type>>& buff) {
      const auto index = 27 * ptype + dir.neighbor_index();
      const auto end   = ends[index];
      const auto begin = index == 0 ? 0uz : ends[index - 1];
      const auto p     = thrust::raw_pointer_cast(buff.data());
      return std::span<const runko::ParticleState<value_type>>(
        std::ranges::next(p, begin),
        std::ranges::next(p, end));
    };

  if(const auto* nonvirtual_other = dynamic_cast<const Tile<D>*>(other_base_ptr)) {
    switch(static_cast<comm_mode>(mode)) {
      case comm_mode::pic_particle:
        for(const auto ptype: rv::iota(0uz, this->particle_buffs_.size())) {
          this->incoming_subregion_particles_[ptype].push_back(get_span(
            ptype,
            inverted_dir,
            nonvirtual_other->subregion_particle_ends_,
            nonvirtual_other->subregion_particle_buff_));
        }
        break;
      default:
        throw std::logic_error { std::format(
          "pic::Tile::local_communication does not support given "
          "communication mode: {}",
          mode) };
    }
  } else if(
    const auto* virtual_other = dynamic_cast<const VirtualTile<D>*>(other_base_ptr)) {
    switch(static_cast<comm_mode>(mode)) {
      case comm_mode::pic_particle:
        for(const auto ptype: rv::iota(0uz, this->particle_buffs_.size())) {
          this->incoming_subregion_particles_[ptype].push_back(get_span(
            ptype,
            inverted_dir,
            virtual_other->subregion_ends_,
            virtual_other->buffer_));
        }
        break;
      default:
        throw std::logic_error { std::format(
          "pic::Tile::local_communication does not support given "
          "communication mode: {}",
          mode) };
    }
  } else {
    throw std::runtime_error {
      "pic::Tile::local_communication assumes that the other tile is "
      "pic::Tile or pic::VirtualTile."
    };
  }
}

}  // namespace pic

template class pic::Tile<3>;
