#include <mpi.h>

#include "io/snapshots/mpiio_fields.h"
#include "io/snapshots/mpiio_header.h"
#include "core/emf/tile.h"
#include "core/pic/tile.h"
#include "tools/simd_math.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <cstring>
#include <string>


//--------------------------------------------------
// write_header (moved from mpiio_header.h)

int mpiio::write_header(
  MPI_File fh,
  int32_t nx, int32_t ny, int32_t nz,
  int32_t stride,
  int32_t Nx, int32_t Ny, int32_t Nz,
  int32_t NxMesh, int32_t NyMesh, int32_t NzMesh,
  int32_t lap)
{
  char buf[256];
  std::memset(buf, 0, sizeof(buf));

  auto put32 = [&](int offset, uint32_t val) {
    std::memcpy(buf + offset, &val, 4);
  };

  auto puti32 = [&](int offset, int32_t val) {
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

  // Field names: 11 entries of 16 chars each, starting at offset 64
  for (int f = 0; f < num_fields; f++) {
    std::strncpy(buf + 64 + f * 16, field_names[f], 15);
  }

  // [240:256] already zeroed

  MPI_Status status;
  return MPI_File_write_at(fh, 0, buf, 256, MPI_BYTE, &status);
}


//--------------------------------------------------

template<>
mpiio::FieldsWriter<3>::FieldsWriter(
  const std::string& prefix,
  int Nx, int NxMesh,
  int Ny, int NyMesh,
  int Nz, int NzMesh,
  int stride)
  : prefix_(prefix),
    Nx_(Nx), Ny_(Ny), Nz_(Nz),
    NxMesh_(NxMesh), NyMesh_(NyMesh), NzMesh_(NzMesh),
    stride_(stride),
    iter_grid_(0, 0, 0),
    tile_buf_(0, 0, 0)
{
  nxt_ = std::max(1, NxMesh_ / stride_);
  nyt_ = std::max(1, NyMesh_ / stride_);
  nzt_ = std::max(1, NzMesh_ / stride_);

  nx_ = Nx_ * nxt_;
  ny_ = Ny_ * nyt_;
  nz_ = Nz_ * nzt_;

  iter_grid_.invalidating_resize(
    static_cast<runko::size_t>(nzt_),
    static_cast<runko::size_t>(nyt_),
    static_cast<runko::size_t>(nxt_));

  tile_buf_.invalidating_resize(
    static_cast<runko::size_t>(nzt_),
    static_cast<runko::size_t>(nyt_),
    static_cast<runko::size_t>(nxt_));
}


template<>
void mpiio::FieldsWriter<3>::pack_tile(emf::Tile<3>& tile)
{
  const auto [Emds, Bmds, Jmds] = tile.view_EBJ_on_device();

  const int stride = stride_;
  auto buf_mds = tile_buf_.mds();

  tyvi::mdgrid_work w {};

  // E and B fields: subsample by stride-hopping; zero n0, n1
  w.for_each_index(iter_grid_, [=](const auto idx) {
    const auto iz = idx[0], iy = idx[1], ix = idx[2];

    const auto si = static_cast<runko::size_t>(ix * stride);
    const auto sj = static_cast<runko::size_t>(iy * stride);
    const auto sk = static_cast<runko::size_t>(iz * stride);

    // ex, ey, ez
    buf_mds[iz, iy, ix][0] = Emds[si, sj, sk][0];
    buf_mds[iz, iy, ix][1] = Emds[si, sj, sk][1];
    buf_mds[iz, iy, ix][2] = Emds[si, sj, sk][2];

    // bx, by, bz
    buf_mds[iz, iy, ix][3] = Bmds[si, sj, sk][0];
    buf_mds[iz, iy, ix][4] = Bmds[si, sj, sk][1];
    buf_mds[iz, iy, ix][5] = Bmds[si, sj, sk][2];

    // n0, n1 (overwritten below if PIC tile)
    buf_mds[iz, iy, ix][9]  = 0.0f;
    buf_mds[iz, iy, ix][10] = 0.0f;
  });

  // J fields: volume-sum over stride^3 cells
  w.for_each_index(iter_grid_, [=](const auto idx) {
    const auto iz = idx[0], iy = idx[1], ix = idx[2];

    float sjx = 0.0f, sjy = 0.0f, sjz = 0.0f;
    for (int kk = 0; kk < stride; kk++)
      for (int jj = 0; jj < stride; jj++)
        for (int ii = 0; ii < stride; ii++) {
          const auto si = static_cast<runko::size_t>(ix*stride + ii);
          const auto sj = static_cast<runko::size_t>(iy*stride + jj);
          const auto sk = static_cast<runko::size_t>(iz*stride + kk);
          sjx += Jmds[si, sj, sk][0];
          sjy += Jmds[si, sj, sk][1];
          sjz += Jmds[si, sj, sk][2];
        }

    buf_mds[iz, iy, ix][6] = sjx;
    buf_mds[iz, iy, ix][7] = sjy;
    buf_mds[iz, iy, ix][8] = sjz;
  });

  w.wait();

  // Number density: deposit particles directly at coarse resolution.
  // n0, n1 slots are already zeroed in the EB kernel above.
  auto* pic_tile = dynamic_cast<pic::Tile<3>*>(&tile);
  if (pic_tile) {
    using vt = pic::ParticleContainer::value_type;
    const auto mx = static_cast<vt>(tile.mins[0]);
    const auto my = static_cast<vt>(tile.mins[1]);
    const auto mz = static_cast<vt>(tile.mins[2]);
    const auto inv_stride = vt{1} / static_cast<vt>(stride);

    auto deposit_species = [&](runko::size_t species, int f_idx) {
      const auto pos_mds = pic_tile->particles(species).pos_mds();
      const auto fi = static_cast<runko::size_t>(f_idx);

      tyvi::mdgrid_work {}
        .for_each_index(
          pos_mds,
          [=](const auto idx) {
            const auto px = pos_mds[idx][0] - mx;
            const auto py = pos_mds[idx][1] - my;
            const auto pz = pos_mds[idx][2] - mz;

            const auto ci = static_cast<runko::size_t>(sstd::floor(px * inv_stride));
            const auto cj = static_cast<runko::size_t>(sstd::floor(py * inv_stride));
            const auto ck = static_cast<runko::size_t>(sstd::floor(pz * inv_stride));

            auto* const n = &thrust::raw_reference_cast(buf_mds[ck, cj, ci][fi]);
            sstd::atomic_add(n, vt{1});
          })
        .wait();
    };

    const auto nspec = pic_tile->number_of_species();
    if (nspec > 0) deposit_species(0, 9);
    if (nspec > 1) deposit_species(1, 10);
  }

#if defined(TYVI_BACKEND_HIP)
  // GPU: copy packed data from device to host staging buffer
  tyvi::mdgrid_work w2 {};
  w2.sync_to_staging(tile_buf_);
  w2.wait();
#endif
  // CPU: device buffer is already host-accessible; no copy needed.
}


template<>
bool mpiio::FieldsWriter<3>::write(corgi::Grid<3>& grid, int lap)
{
  std::string filename =
    prefix_ + "/flds_" + std::to_string(lap) + ".bin";

  MPI_File fh;
  int rc = MPI_File_open(
    MPI_Comm(grid.comm), filename.c_str(),
    MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &fh);

  if (rc != MPI_SUCCESS) return false;

  // Rank 0 writes the 256-byte header.
  // MPI_File_open is collective and synchronizing; no barrier needed
  // since all data writes target non-overlapping offsets >= 256.
  if (grid.comm.rank() == 0) {
    rc = mpiio::write_header(
      fh,
      nx_, ny_, nz_,
      stride_,
      Nx_, Ny_, Nz_,
      NxMesh_, NyMesh_, NzMesh_,
      lap);
    if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
  }

  const MPI_Offset hdr_size = mpiio::header_size;
  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(nx_) * ny_ * nz_ * sizeof(float);
  const int tile_elems = nxt_ * nyt_ * nzt_;

  for (auto cid : grid.get_local_tiles()) {
    auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile(cid));
    pack_tile(tile);

    // On CPU the device buffer is host-accessible; on GPU use staging.
#if defined(TYVI_BACKEND_CPU)
    const float* write_ptr = tile_buf_.span().data();
#elif defined(TYVI_BACKEND_HIP)
    const float* write_ptr = tile_buf_.staging_span().data();
#endif

    const int ti = static_cast<int>(tile.index[0]);
    const int tj = static_cast<int>(tile.index[1]);
    const int tk = static_cast<int>(tile.index[2]);

    // Write each field for this tile
    for (int f = 0; f < num_fields; f++) {
      const MPI_Offset field_base = hdr_size + f * field_bytes;

      // Write row by row (each row of nxt_ floats is contiguous in file)
      for (int ks = 0; ks < nzt_; ks++) {
        for (int js = 0; js < nyt_; js++) {
          const MPI_Offset file_offset = field_base +
            static_cast<MPI_Offset>(
              (tk*nzt_ + ks) * ny_ * nx_ +
              (tj*nyt_ + js) * nx_ +
               ti*nxt_) * sizeof(float);

          const int buf_offset =
            f * tile_elems + ks*nyt_*nxt_ + js*nxt_;

          MPI_Status status;
          rc = MPI_File_write_at(
            fh, file_offset,
            write_ptr + buf_offset,
            nxt_, MPI_FLOAT, &status);
          if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
        }
      }
    }
  }

  rc = MPI_File_close(&fh);
  return rc == MPI_SUCCESS;
}


//--------------------------------------------------
// explicit template class instantiation
template class mpiio::FieldsWriter<3>;
