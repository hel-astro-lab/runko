#include "core/pic/virtual_tile.h"

#include "thrust/memory.h"
#include "tools/system.h"

#include <format>
#include <print>

namespace pic {

template<std::size_t D>
VirtualTile<D>::VirtualTile(
  const std::array<std::size_t, 3> tile_indices,
  const toolbox::ConfigParser& p) :
  emf::VirtualTile<D>(tile_indices, p)
{

  for(auto i = 0uz; true; ++i) {
    const auto q_label = std::format("q{}", i);
    const auto m_label = std::format("m{}", i);

    const auto q = p.get<double>(q_label);
    const auto m = p.get<double>(m_label);

    if(q and m) {
      amount_of_ptypes_ = i + 1;
    } else {
      break;
    }
  }
}

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  VirtualTile<D>::recv_data(
    mpi4cpp::mpi::communicator& comm,
    const int orig,
    const int mode,
    const int tag)
{
#ifndef TYVI_BACKEND_CPU
  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error { "GPU backend requires GPU-aware MPI." };
  }
#endif

  using runko::comm_mode;

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::number_of_particles:
      this->subregion_ends_.resize(this->amount_of_ptypes_ * 27);
      return std::vector<mpi4cpp::mpi::request> { comm.irecv(
        orig,
        tag,
        thrust::raw_pointer_cast(subregion_ends_.data()),
        this->subregion_ends_.size()) };
    case comm_mode::pic_particle: {
      this->buffer_.resize(this->subregion_ends_.back());
      const auto p = thrust::raw_pointer_cast(this->buffer_.data());
      const auto s = std::span(p, this->buffer_.size());
      return std::vector<mpi4cpp::mpi::request> { comm.irecv(
        orig,
        tag,
        reinterpret_cast<char*>(std::as_writable_bytes(s).data()),
        s.size_bytes()) };
    }
    default: return emf::VirtualTile<D>::recv_data(comm, orig, mode, tag);
  }
}
}  // namespace pic

template class pic::VirtualTile<3>;
