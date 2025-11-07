#include "core/communication_common.h"
#include "core/pic/tile.h"
#include "tools/system.h"

#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace {

constexpr int
  particle_tag(
    const int tag,
    const auto particle_ordinal,
    const std::size_t particles_in_total)
{
  return static_cast<int>(particle_ordinal + particles_in_total * tag);
}

enum class particle_data_type { pos, vel };

template<particle_data_type data_type>
constexpr int
  particle_data_tag(
    const int tag,
    const auto particle_ordinal,
    const std::size_t particles_in_total)
{
  static constexpr auto x = static_cast<int>(data_type == particle_data_type::pos);
  return static_cast<int>(x + 2 * particle_ordinal + 2 * particles_in_total * tag);
}

}  // namespace

namespace pic {

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  Tile<D>::send_data(
    mpi4cpp::mpi::communicator& comm,
    const int dest,
    const int mode,
    const int tag)
{

  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error {
      "Non gpu aware MPI communication is not yet implemented."
    };
  }

  /*
    All tiles should have the particle containers
    for all particles that are configure in config.
    This means that we can assume that all tiles should be
    in agreement on what particles there are.
   */


  auto make_isend = [&](const auto s, const int t) {
    return comm.isend(dest, t, s.data(), s.size());
  };

  using runko::comm_mode;

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::number_of_particles: {
      auto comms = std::vector<mpi4cpp::mpi::request>();
      for(const auto& [key, value]: this->particle_buffs_) {
        this->amount_of_particles_to_be_send_[key] = value.size();
        comms.push_back(comm.isend(
          dest,
          particle_tag(tag, key, this->particle_buffs_.size()),
          this->amount_of_particles_to_be_send_[key]));
      }
      return comms;
    }
    case comm_mode::pic_particle: {
      auto comms = std::vector<mpi4cpp::mpi::request>();
      for(const auto& [key, value]: this->particle_buffs_) {
        comms.push_back(make_isend(
          value.span_pos(),
          particle_data_tag<particle_data_type::pos>(
            tag,
            key,
            this->particle_buffs_.size())));
        comms.push_back(make_isend(
          value.span_vel(),
          particle_data_tag<particle_data_type::vel>(
            tag,
            key,
            this->particle_buffs_.size())));
      }
      return comms;
    }
    default: return emf::Tile<D>::send_data(comm, dest, mode, tag);
  }
}

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  Tile<D>::recv_data(
    mpi4cpp::mpi::communicator& comm,
    const int orig,
    const int mode,
    const int tag)
{
  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error {
      "Non gpu aware MPI communication is not yet implemented."
    };
  }

  using runko::comm_mode;

  const auto recv_particle_sizes = [&, this] {
    auto comms = std::vector<mpi4cpp::mpi::request>();
    for(const auto& [key, value]: this->particle_buffs_) {
      comms.push_back(comm.irecv(
        orig,
        particle_tag(tag, key, this->particle_buffs_.size()),
        this->amount_of_particles_to_be_received_[key]));
    }
    return comms;
  };


  auto make_irecv = [&](const auto s, const int t) {
    return comm.irecv(orig, t, s.data(), s.size());
  };

  const auto recv_particles = [&, this] {
    auto comms = std::vector<mpi4cpp::mpi::request>();
    for(auto& [key, value]: this->particle_buffs_) {
      auto new_args = value.args();
      new_args.N    = amount_of_particles_to_be_received_[key];
      value         = pic::ParticleContainer(new_args);

      comms.push_back(make_irecv(
        value.span_pos(),
        particle_data_tag<particle_data_type::pos>(
          tag,
          key,
          this->particle_buffs_.size())));
      comms.push_back(make_irecv(
        value.span_vel(),
        particle_data_tag<particle_data_type::vel>(
          tag,
          key,
          this->particle_buffs_.size())));
    }
    return comms;
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::number_of_particles: return recv_particle_sizes();
    case comm_mode::pic_particle: return recv_particles();
    default: return emf::Tile<D>::recv_data(comm, orig, mode, tag);
  }
}

template<std::size_t D>
void
  Tile<D>::divide_particles_to_subregions()
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

  for(auto& [ptype, pcontainer]: this->particle_buffs_) {
    auto [s, c] = pcontainer.divide_to_subregions(x_div, y_div, z_div);
    std::ignore = this->subregion_particle_spans_.insert_or_assign(ptype, std::move(s));
    std::ignore = this->subregion_particle_buffs_.insert_or_assign(ptype, std::move(c));
  }
}

template<std::size_t D>
void
  Tile<D>::pairwise_moore_communication_prelude(const int mode)
{
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) { return; }

  this->divide_particles_to_subregions();

  // Make particles staying in this tile "incoming" particles.
  for(const auto& [ptype, spans]: this->subregion_particle_spans_) {
    incoming_subregion_particles_[ptype] =
      std::vector { ParticleContainer::specific_span {
        .span      = spans.at({ 0, 0, 0 }),
        .container = &this->subregion_particle_buffs_.at(ptype) } };
  }
}

template<std::size_t D>
void
  Tile<D>::pairwise_moore_communication_postlude(const int mode)
{
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) { return; }

  for(const auto& [ptype, specific_spans]: this->incoming_subregion_particles_) {
    particle_buffs_.at(ptype).set(specific_spans);
  }

  for(auto& [_, pbuff]: this->particle_buffs_) {
    pbuff.wrap_positions(this->global_coordinate_mins_, this->global_coordinate_maxs_);
  }

  this->incoming_subregion_particles_.clear();
  this->subregion_particle_spans_.clear();
  // Note that this->subregion_particle_buffs_ can not be cleared here.
  // Other tiles might be using it.
}

template<std::size_t D>
void
  Tile<D>::pairwise_moore_communication(
    const corgi::Tile<D>& other_base,
    const std::array<int, D> dir_to_other,
    const int mode)
{
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) {
    emf::Tile<D>::pairwise_moore_communication(other_base, dir_to_other, mode);
    return;
  }

  const Tile<D>& other = [&]() -> const Tile<D>& {
    try {
      return dynamic_cast<const Tile<D>&>(other_base);
    } catch(const std::bad_cast& ex) {
      throw std::runtime_error { std::format(
        "pic::Tile::pairwise_moore_communication assumes that the other tile is "
        "pic::Tile or its descendant. Orginal exception: {}",
        ex.what()) };
    }
  }();

  const auto inverted_dir = this->yee_lattice_.invert_dir(dir_to_other);

  for(const auto& [ptype, spans]: other.subregion_particle_spans_) {
    if(spans.contains(inverted_dir)) {
      this->incoming_subregion_particles_.at(ptype).push_back(
        ParticleContainer::specific_span {
          .span      = spans.at(inverted_dir),
          .container = &other.subregion_particle_buffs_.at(ptype) });
    }
  }
};

}  // namespace pic

template class pic::Tile<3>;
