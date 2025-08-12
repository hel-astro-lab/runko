#include "core/communication_common.h"
#include "core/pic2/tile.h"
#include "tools/system.h"

#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace {

constexpr int
  particle_tag(const int tag, const auto particle_ordinal)
{
  return static_cast<int>(particle_ordinal) + tag;
}

enum class particle_data_type { pos, vel };

template<particle_data_type data_type>
constexpr int
  particle_data_tag(const int tag, const auto particle_ordinal)
{
  static constexpr auto x = static_cast<int>(data_type == particle_data_type::pos);
  return 2 * static_cast<int>(particle_ordinal) + tag + x;
}

}  // namespace

namespace pic2 {

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
        comms.push_back(comm.isend(dest, particle_tag(tag, key), value.size()));
      }
      return comms;
    }
    case comm_mode::pic_particle: {
      auto comms = std::vector<mpi4cpp::mpi::request>();
      for(const auto& [key, value]: this->particle_buffs_) {
        comms.push_back(make_isend(
          value.span_pos(),
          particle_data_tag<particle_data_type::pos>(tag, key)));
        comms.push_back(make_isend(
          value.span_vel(),
          particle_data_tag<particle_data_type::vel>(tag, key)));
      }
      return comms;
    }
    default: return emf2::Tile<D>::send_data(comm, dest, mode, tag);
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
        particle_tag(tag, key),
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
      value         = pic2::ParticleContainer(new_args);

      comms.push_back(make_irecv(
        value.span_pos(),
        particle_data_tag<particle_data_type::pos>(tag, key)));
      comms.push_back(make_irecv(
        value.span_vel(),
        particle_data_tag<particle_data_type::vel>(tag, key)));
    }
    return comms;
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::number_of_particles: return recv_particle_sizes();
    case comm_mode::pic_particle: return recv_particles();
    default: return emf2::Tile<D>::recv_data(comm, orig, mode, tag);
  }
}

template<std::size_t D>
void
  Tile<D>::split_particles_to_subregions()
{
  auto buffs = subregion_particle_buff {};

  const auto x_div =
    std::array { static_cast<ParticleContainer::value_type>(this->mins[0]),
                 static_cast<ParticleContainer::value_type>(this->maxs[0]) };
  const auto y_div =
    std::array { static_cast<ParticleContainer::value_type>(this->mins[1]),
                 static_cast<ParticleContainer::value_type>(this->maxs[1]) };
  const auto z_div =
    std::array { static_cast<ParticleContainer::value_type>(this->mins[2]),
                 static_cast<ParticleContainer::value_type>(this->maxs[2]) };

  for(auto& [key, value]: this->particle_buffs_) {
    for(auto& [dir, pbuff]: value.split_to_subregions(
          x_div,
          y_div,
          z_div,
          global_coordinate_mins_,
          global_coordinate_maxs_)) {
      buffs[dir].insert_or_assign(key, std::move(pbuff));
    }
  }

  this->subregion_particles_ = std::move(buffs);
}

template<std::size_t D>
void
  Tile<D>::pairwise_moore_communication_prelude(const int mode)
{
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) { return; }

  this->split_particles_to_subregions();
}

template<std::size_t D>
void
  Tile<D>::pairwise_moore_communication_postlude(const int mode)
{
  using runko::comm_mode;
  if(static_cast<comm_mode>(mode) != comm_mode::pic_particle) { return; }

  for(auto& [ptype, incoming_buffs]: this->incoming_subregion_particles_) {

    // Due to definition of pairwise moore we know that there is
    // going to be 26 incoming buffs (in case of D == 3 which is assumed).
    static constexpr auto Nincoming = 26uz;
    if(Nincoming != incoming_buffs.size()) {
      throw std::logic_error {
        "pic2::Tile::pairwise_moore_communicaion_postlude: unexpect amount of incoming "
        "particle buffs."
      };
    }

    // Ugly syntax until C++26 :/
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      this->particle_buffs_.at(ptype).add_particles(incoming_buffs[I]...);
    }(std::make_index_sequence<Nincoming>());

    incoming_buffs.clear();
  }
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
    emf2::Tile<D>::pairwise_moore_communication(other_base, dir_to_other, mode);
    return;
  }

  const Tile<D>& other = [&]() -> const Tile<D>& {
    try {
      return dynamic_cast<const Tile<D>&>(other_base);
    } catch(const std::bad_cast& ex) {
      throw std::runtime_error { std::format(
        "pic2::Tile::pairwise_moore_communication assumes that the other tile is "
        "pic2::Tile or its descendant. Orginal exception: {}",
        ex.what()) };
    }
  }();

  const auto inverted_dir = this->yee_lattice_.invert_dir(dir_to_other);

  for(const auto& [ptype, pbuff]: other.subregion_particles_.at(inverted_dir)) {
    this->incoming_subregion_particles_[ptype].push_back(std::cref(pbuff));
  }
};

}  // namespace pic2

template class pic2::Tile<3>;
