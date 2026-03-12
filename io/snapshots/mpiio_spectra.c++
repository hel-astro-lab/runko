#include <mpi.h>

#include "io/snapshots/mpiio_spectra.h"
#include "io/snapshots/mpiio_header.h"
#include "core/emf/tile.h"
#include "core/pic/tile.h"
#include "tools/simd_math.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <format>
#include <string>


//--------------------------------------------------
// write_spectra_header

int mpiio::write_spectra_header(
  MPI_File fh,
  int32_t nx, int32_t ny, int32_t nz,
  int32_t stride,
  int32_t Nx, int32_t Ny, int32_t Nz,
  int32_t NxMesh, int32_t NyMesh, int32_t NzMesh,
  int32_t lap,
  int32_t num_fields,
  int32_t nbins,
  float umin,
  float umax)
{
  char buf[512];
  std::memset(buf, 0, sizeof(buf));

  auto put32 = [&](int offset, uint32_t val) {
    std::memcpy(buf + offset, &val, 4);
  };

  auto puti32 = [&](int offset, int32_t val) {
    std::memcpy(buf + offset, &val, 4);
  };

  auto putf32 = [&](int offset, float val) {
    std::memcpy(buf + offset, &val, 4);
  };

  put32(0, magic);
  put32(4, version);
  put32(8, static_cast<uint32_t>(header_size));
  put32(12, static_cast<uint32_t>(num_fields));

  puti32(16, nx);
  puti32(20, ny);
  puti32(24, nz);
  puti32(28, stride);
  puti32(32, Nx);
  puti32(36, Ny);
  puti32(40, Nz);
  puti32(44, NxMesh);
  puti32(48, NyMesh);
  puti32(52, NzMesh);
  puti32(56, lap);
  put32(60, 4);  // dtype_size = sizeof(float)

  // Field names: "s0_u", "s0_bx", "s0_by", "s0_bz", "s1_u", ...
  const int nspec = num_fields / num_spectra_per_species;
  for (int s = 0; s < nspec; s++) {
    for (int q = 0; q < num_spectra_per_species; q++) {
      auto name = std::format("s{}_{}", s, spectra_suffixes[q]);
      int f = s * num_spectra_per_species + q;
      std::strncpy(buf + 64 + f * 16, name.c_str(), 15);
    }
  }

  // Spectra-specific metadata at byte 256
  puti32(256, nbins);
  putf32(260, umin);
  putf32(264, umax);

  MPI_Status status;
  return MPI_File_write_at(fh, 0, buf, 512, MPI_BYTE, &status);
}


//--------------------------------------------------
// SpectraWriter constructor

template<>
mpiio::SpectraWriter<3>::SpectraWriter(
  const std::string& prefix,
  int Nx, int NxMesh,
  int Ny, int NyMesh,
  int Nz, int NzMesh,
  int stride,
  int nbins,
  float umin,
  float umax,
  int nspecies)
  : WriterBase<3>(prefix,
      std::min(nspecies, static_cast<int>(max_spectra_species)) * num_spectra_per_species),
    Nx_(Nx), Ny_(Ny), Nz_(Nz),
    NxMesh_(NxMesh), NyMesh_(NyMesh), NzMesh_(NzMesh),
    stride_(stride),
    nbins_(nbins),
    umin_(umin), umax_(umax),
    nspecies_(std::min(nspecies, static_cast<int>(max_spectra_species))),
    tile_buf_(0, 0)
{
  nxt_ = std::max(1, NxMesh_ / stride_);
  nx_global_ = Nx_ * nxt_;

  log_umin_ = std::log10(umin_);
  float log_umax = std::log10(umax_);
  inv_dlog_ = static_cast<float>(nbins_) / (log_umax - log_umin_);
  inv_dbeta_ = static_cast<float>(nbins_) / 2.0f;

  tile_buf_.invalidating_resize(
    static_cast<runko::index_t>(nxt_),
    static_cast<runko::index_t>(nbins_));
}


//--------------------------------------------------
// histogram_tile

template<>
void mpiio::SpectraWriter<3>::histogram_tile(emf::Tile<3>& tile)
{
  auto buf_mds = tile_buf_.mds();
  const int nf = num_fields_;

  // Zero the histogram buffer
  tyvi::mdgrid_work {}
    .for_each_index(tile_buf_, [=](const auto idx) {
      for (int f = 0; f < nf; f++)
        buf_mds[idx][f] = 0.0f;
    })
    .wait();

  auto* pic_tile = dynamic_cast<pic::Tile<3>*>(&tile);
  if (!pic_tile) return;

  using vt = pic::ParticleContainer::value_type;
  const auto mx         = static_cast<vt>(tile.mins[0]);
  const auto inv_stride = vt{1} / static_cast<vt>(stride_);
  const auto log_umin   = static_cast<vt>(log_umin_);
  const auto inv_dlog   = static_cast<vt>(inv_dlog_);
  const auto inv_dbeta  = static_cast<vt>(inv_dbeta_);
  const auto last_bin_f = static_cast<vt>(nbins_ - 1);

  const auto nspec = static_cast<int>(pic_tile->number_of_species());
  const int ndeposit = std::min(nspec, nspecies_);

  for (int s = 0; s < ndeposit; s++) {
    const auto base_field = static_cast<runko::index_t>(s * num_spectra_per_species);
    const auto pos_mds = pic_tile->particles(static_cast<runko::index_t>(s)).pos_mds();
    const auto vel_mds = pic_tile->particles(static_cast<runko::index_t>(s)).vel_mds();

    tyvi::mdgrid_work {}
      .for_each_index(
        pos_mds,
        [=](const auto idx) {
          // particle x position -> x bin
          const auto px = pos_mds[idx][0] - mx;
          const auto ix = static_cast<runko::index_t>(sstd::floor(px * inv_stride));

          // particle velocities
          const auto ux = vel_mds[idx][0];
          const auto uy = vel_mds[idx][1];
          const auto uz = vel_mds[idx][2];

          const auto u2    = ux*ux + uy*uy + uz*uz;
          const auto u_mag = sstd::sqrt(u2);
          const auto gamma = sstd::sqrt(vt{1} + u2);
          const auto inv_gamma = vt{1} / gamma;

          // branchless clamp to [0, last_bin]; out-of-range -> boundary bin
          auto bin = [=](vt raw) -> runko::index_t {
            return static_cast<runko::index_t>(sstd::clamp(raw, vt{0}, last_bin_f));
          };

          // u spectrum: log10 bins [umin, umax]
          const auto ib_u  = bin(sstd::floor((sstd::log10(u_mag) - log_umin) * inv_dlog));

          // beta spectra: linear bins [-1, +1]
          const auto ib_bx = bin(sstd::floor((ux * inv_gamma + vt{1}) * inv_dbeta));
          const auto ib_by = bin(sstd::floor((uy * inv_gamma + vt{1}) * inv_dbeta));
          const auto ib_bz = bin(sstd::floor((uz * inv_gamma + vt{1}) * inv_dbeta));

          // unconditional atomic deposits
          sstd::atomic_add(&thrust::raw_reference_cast(buf_mds[ix, ib_u ][base_field + 0]), vt{1});
          sstd::atomic_add(&thrust::raw_reference_cast(buf_mds[ix, ib_bx][base_field + 1]), vt{1});
          sstd::atomic_add(&thrust::raw_reference_cast(buf_mds[ix, ib_by][base_field + 2]), vt{1});
          sstd::atomic_add(&thrust::raw_reference_cast(buf_mds[ix, ib_bz][base_field + 3]), vt{1});
        })
      .wait();
  }
}


//--------------------------------------------------
// NVI overrides

template<>
std::string mpiio::SpectraWriter<3>::make_filename_(int lap) const
{
  return prefix_ + "/pspectra_" + std::to_string(lap) + ".bin";
}


template<>
int mpiio::SpectraWriter<3>::write_header_(MPI_File fh, int lap)
{
  return mpiio::write_spectra_header(
    fh,
    nx_global_, Ny_, Nz_,
    stride_,
    Nx_, Ny_, Nz_,
    NxMesh_, NyMesh_, NzMesh_,
    lap,
    num_fields_,
    nbins_,
    umin_, umax_);
}


template<>
bool mpiio::SpectraWriter<3>::write_payload_(MPI_File fh, corgi::Grid<3>& grid)
{
  const MPI_Offset hdr_size = mpiio::header_size;
  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(Nz_) * Ny_ * nx_global_ * nbins_ * sizeof(float);
  const int tile_row_elems = nxt_ * nbins_;

  // Pre-allocate the file to its full size
  MPI_Offset total_size = hdr_size + static_cast<MPI_Offset>(num_fields_) * field_bytes;
  MPI_File_set_size(fh, total_size);

  int rc;
  for (auto cid : grid.get_local_tiles()) {
    auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile(cid));
    histogram_tile(tile);

    // On CPU the device buffer is host-accessible; on GPU use staging.
#if defined(TYVI_BACKEND_CPU)
    const float* write_ptr = tile_buf_.span().data();
#elif defined(TYVI_BACKEND_HIP)
    tyvi::mdgrid_work {}.sync_to_staging(tile_buf_).wait();
    const float* write_ptr = tile_buf_.staging_span().data();
#endif

    const int ti = static_cast<int>(tile.index[0]);
    const int tj = static_cast<int>(tile.index[1]);
    const int tk = static_cast<int>(tile.index[2]);

    // Write each field: one contiguous block of nxt_*nbins_ floats
    for (int f = 0; f < num_fields_; f++) {
      const MPI_Offset field_base = hdr_size + f * field_bytes;
      const MPI_Offset file_offset = field_base +
        static_cast<MPI_Offset>(
          (tk * Ny_ + tj) * nx_global_ + ti * nxt_) * nbins_ * sizeof(float);

      const int buf_offset = f * tile_row_elems;

      MPI_Status status;
      rc = MPI_File_write_at(
        fh, file_offset,
        write_ptr + buf_offset,
        tile_row_elems, MPI_FLOAT, &status);
      if (rc != MPI_SUCCESS) return false;
    }
  }

  return true;
}


//--------------------------------------------------
// explicit template class instantiation
template class mpiio::SpectraWriter<3>;
