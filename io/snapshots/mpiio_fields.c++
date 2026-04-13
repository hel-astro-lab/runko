#include <mpi.h>

#include "io/snapshots/mpiio_fields.h"
#include "io/snapshots/mpiio_header.h"
#include "core/emf/tile.h"
#include "core/pic/tile.h"
#include "tools/math.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <cstring>
#include <format>
#include <string>
#include <vector>


//--------------------------------------------------
// write_header

int mpiio::write_header(
  MPI_File fh,
  int32_t nx, int32_t ny, int32_t nz,
  int32_t stride,
  int32_t Nx, int32_t Ny, int32_t Nz,
  int32_t NxMesh, int32_t NyMesh, int32_t NzMesh,
  int32_t lap,
  int32_t num_fields)
{
  char buf[512];
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

  // Field names: num_fields entries of 16 chars each, starting at offset 64
  for (int f = 0; f < num_emf_fields; f++) {
    std::strncpy(buf + 64 + f * 16, emf_field_names[f], 15);
  }
  for (int s = 0; s < num_fields - num_emf_fields; s++) {
    auto name = std::format("n{}", s);
    std::strncpy(buf + 64 + (num_emf_fields + s) * 16, name.c_str(), 15);
  }

  MPI_Status status;
  return MPI_File_write_at(fh, 0, buf, 512, MPI_BYTE, &status);
}


//--------------------------------------------------

template<>
mpiio::FieldsWriter<3>::FieldsWriter(
  const std::string& prefix,
  int Nx, int NxMesh,
  int Ny, int NyMesh,
  int Nz, int NzMesh,
  int stride,
  int nspecies)
  : WriterBase<3>(prefix, num_emf_fields + std::min(nspecies, static_cast<int>(max_species))),
    Nx_(Nx), Ny_(Ny), Nz_(Nz),
    NxMesh_(NxMesh), NyMesh_(NyMesh), NzMesh_(NzMesh),
    stride_(stride),
    nspecies_(std::min(nspecies, static_cast<int>(max_species))),
    iter_grid_(0, 0, 0),
    tile_buf_(0, 0, 0),
    cached_empty_type_(MPI_DATATYPE_NULL),
    cached_max_ntiles_(0)
{
  nxt_ = std::max(1, NxMesh_ / stride_);
  nyt_ = std::max(1, NyMesh_ / stride_);
  nzt_ = std::max(1, NzMesh_ / stride_);

  nx_ = Nx_ * nxt_;
  ny_ = Ny_ * nyt_;
  nz_ = Nz_ * nzt_;

  iter_grid_.invalidating_resize(
    static_cast<std::size_t>(nzt_),
    static_cast<std::size_t>(nyt_),
    static_cast<std::size_t>(nxt_));

  tile_buf_.invalidating_resize(
    static_cast<std::size_t>(nzt_),
    static_cast<std::size_t>(nyt_),
    static_cast<std::size_t>(nxt_));
}


template<>
void mpiio::FieldsWriter<3>::free_cached_filetypes_()
{
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized) return;

  for (auto& ft : cached_filetypes_)
    MPI_Type_free(&ft);
  cached_filetypes_.clear();

  if (cached_empty_type_ != MPI_DATATYPE_NULL)
    MPI_Type_free(&cached_empty_type_);

  cached_max_ntiles_ = 0;
  cached_cids_.clear();
}


template<>
mpiio::FieldsWriter<3>::~FieldsWriter<3>()
{
  free_cached_filetypes_();
}


template<>
void mpiio::FieldsWriter<3>::rebuild_filetypes_(corgi::Grid<3>& grid)
{
  free_cached_filetypes_();

  auto local_cids = grid.get_local_tiles();
  const int ntiles = static_cast<int>(local_cids.size());
  const int nf = num_fields_;
  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(nx_) * ny_ * nz_ * sizeof(float);

  // Collective: agree on max tiles across all ranks
  int max_ntiles = 0;
  MPI_Allreduce(&ntiles, &max_ntiles, 1, MPI_INT, MPI_MAX,
                MPI_Comm(grid.comm));

  int gsizes[3]   = {nz_, ny_, nx_};
  int subsizes[3] = {nzt_, nyt_, nxt_};

  // Field displacement table (same for every tile)
  std::vector<MPI_Aint>     field_disps(nf);
  std::vector<int>          field_blens(nf, 1);
  std::vector<MPI_Datatype> ftypes(nf);
  for (int f = 0; f < nf; f++)
    field_disps[f] = static_cast<MPI_Aint>(f) * field_bytes;

  // Build per-tile filetypes
  cached_filetypes_.resize(ntiles);
  for (int t = 0; t < ntiles; t++) {
    auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile(local_cids[t]));
    int starts[3] = {
      static_cast<int>(tile.index[2]) * nzt_,
      static_cast<int>(tile.index[1]) * nyt_,
      static_cast<int>(tile.index[0]) * nxt_
    };

    MPI_Datatype tile_sub;
    MPI_Type_create_subarray(
      3, gsizes, subsizes, starts,
      MPI_ORDER_C, MPI_FLOAT, &tile_sub);
    MPI_Type_commit(&tile_sub);

    std::fill(ftypes.begin(), ftypes.end(), tile_sub);
    MPI_Type_create_struct(
      nf, field_blens.data(), field_disps.data(),
      ftypes.data(), &cached_filetypes_[t]);
    MPI_Type_commit(&cached_filetypes_[t]);

    MPI_Type_free(&tile_sub);
  }

  // Single empty type for non-participating rounds
  MPI_Type_contiguous(0, MPI_FLOAT, &cached_empty_type_);
  MPI_Type_commit(&cached_empty_type_);

  cached_max_ntiles_ = max_ntiles;
  cached_cids_ = std::move(local_cids);
}


template<>
void mpiio::FieldsWriter<3>::pack_tile(emf::Tile<3>& tile)
{
  const auto [Emds, Bmds, Jmds] = tile.view_EBJ_on_device();

  const int stride = stride_;
  const int nf = num_fields_;
  auto buf_mds = tile_buf_.mds();

  tyvi::mdgrid_work w {};

  // E and B fields: subsample by stride-hopping; zero density slots
  w.for_each_index(iter_grid_, [=](const auto idx) {
    const auto iz = idx[0], iy = idx[1], ix = idx[2];

    const auto si = static_cast<std::size_t>(ix * stride);
    const auto sj = static_cast<std::size_t>(iy * stride);
    const auto sk = static_cast<std::size_t>(iz * stride);

    // ex, ey, ez
    buf_mds[iz, iy, ix][0] = Emds[si, sj, sk][0];
    buf_mds[iz, iy, ix][1] = Emds[si, sj, sk][1];
    buf_mds[iz, iy, ix][2] = Emds[si, sj, sk][2];

    // bx, by, bz
    buf_mds[iz, iy, ix][3] = Bmds[si, sj, sk][0];
    buf_mds[iz, iy, ix][4] = Bmds[si, sj, sk][1];
    buf_mds[iz, iy, ix][5] = Bmds[si, sj, sk][2];

    // density slots: zero n0..n{nspecies-1} (overwritten below if PIC tile)
    for (int s = 0; s < nf - num_emf_fields; s++)
      buf_mds[iz, iy, ix][num_emf_fields + s] = 0.0f;
  });

  // J fields: volume-sum over stride^3 cells
  w.for_each_index(iter_grid_, [=](const auto idx) {
    const auto iz = idx[0], iy = idx[1], ix = idx[2];

    float sjx = 0.0f, sjy = 0.0f, sjz = 0.0f;
    for (int kk = 0; kk < stride; kk++)
      for (int jj = 0; jj < stride; jj++)
        for (int ii = 0; ii < stride; ii++) {
          const auto si = static_cast<std::size_t>(ix*stride + ii);
          const auto sj = static_cast<std::size_t>(iy*stride + jj);
          const auto sk = static_cast<std::size_t>(iz*stride + kk);
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
  // Density slots are already zeroed in the EB kernel above.
  auto* pic_tile = dynamic_cast<pic::Tile<3>*>(&tile);
  if (pic_tile) {
    using vt = pic::ParticleContainer::value_type;
    const auto mx = static_cast<vt>(tile.mins[0]);
    const auto my = static_cast<vt>(tile.mins[1]);
    const auto mz = static_cast<vt>(tile.mins[2]);
    const auto inv_stride = vt{1} / static_cast<vt>(stride);

    auto deposit_species = [&](std::size_t species, int f_idx) {
      const auto pos_mds = pic_tile->particles(species).pos_mds();
      const auto ids_mds = pic_tile->particles(species).ids_mds();
      const auto fi = static_cast<std::size_t>(f_idx);

      tyvi::mdgrid_work {}
        .for_each_index(
          pos_mds,
          [=](const auto idx) {
            if (ids_mds[idx][] == runko::dead_prtc_id) { return; }
            const auto px = pos_mds[idx][0] - mx;
            const auto py = pos_mds[idx][1] - my;
            const auto pz = pos_mds[idx][2] - mz;

            const auto ci = static_cast<std::size_t>(sstd::floor(px * inv_stride));
            const auto cj = static_cast<std::size_t>(sstd::floor(py * inv_stride));
            const auto ck = static_cast<std::size_t>(sstd::floor(pz * inv_stride));

            auto* const n = &thrust::raw_reference_cast(buf_mds[ck, cj, ci][fi]);
            sstd::atomic_add(n, vt{1});
          })
        .wait();
    };

    const auto nspec = static_cast<int>(pic_tile->number_of_species());
    const int ndeposit = std::min(nspec, nspecies_);
    for (int s = 0; s < ndeposit; s++)
      deposit_species(static_cast<std::size_t>(s), 9 + s);
  }

  // Data remains on device. Callers sync to staging if needed.
}


//--------------------------------------------------
// NVI overrides

template<>
std::string mpiio::FieldsWriter<3>::make_filename_(int lap) const
{
  return prefix_ + "/flds_" + std::to_string(lap) + ".bin";
}


template<>
int mpiio::FieldsWriter<3>::write_header_(MPI_File fh, int lap)
{
  return mpiio::write_header(
    fh,
    nx_, ny_, nz_,
    stride_,
    Nx_, Ny_, Nz_,
    NxMesh_, NyMesh_, NzMesh_,
    lap,
    num_fields_);
}


template<>
bool mpiio::FieldsWriter<3>::write_payload_(MPI_File fh, corgi::Grid<3>& grid)
{
  const MPI_Offset hdr_size = mpiio::header_size;
  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(nx_) * ny_ * nz_ * sizeof(float);
  const int tile_elems = nxt_ * nyt_ * nzt_;

  // Pre-allocate the file to its full size
  MPI_Offset total_size = hdr_size + static_cast<MPI_Offset>(num_fields_) * field_bytes;
  MPI_File_set_size(fh, total_size);

  int rc;
  for (auto cid : grid.get_local_tiles()) {
    auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile(cid));
    pack_tile(tile);

    // On CPU the device buffer is host-accessible; on GPU use staging.
#if defined(TYVI_BACKEND_CPU)
    const float* write_ptr = tile_buf_.span().data();
#elif defined(TYVI_BACKEND_HIP)
    // GPU: copy packed data from device to host staging buffer
    tyvi::mdgrid_work {}.sync_to_staging(tile_buf_).wait();
    const float* write_ptr = tile_buf_.staging_span().data();
#endif

    const int ti = static_cast<int>(tile.index[0]);
    const int tj = static_cast<int>(tile.index[1]);
    const int tk = static_cast<int>(tile.index[2]);

    // Write each field for this tile
    for (int f = 0; f < num_fields_; f++) {
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
          if (rc != MPI_SUCCESS) return false;
        }
      }
    }
  }

  return true;
}


//--------------------------------------------------
// write_collective: bulk MPI-IO using derived file types

template<>
bool mpiio::FieldsWriter<3>::write_collective(corgi::Grid<3>& grid, int lap)
{
  std::string filename =
    prefix_ + "/flds_" + std::to_string(lap) + ".bin";

  MPI_File fh;
  int rc = MPI_File_open(
    MPI_Comm(grid.comm), filename.c_str(),
    MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &fh);

  if (rc != MPI_SUCCESS) return false;

  // Rank 0 writes the 512-byte header.
  if (grid.comm.rank() == 0) {
    rc = mpiio::write_header(
      fh,
      nx_, ny_, nz_,
      stride_,
      Nx_, Ny_, Nz_,
      NxMesh_, NyMesh_, NzMesh_,
      lap,
      num_fields_);
    if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
  }

  const MPI_Offset hdr_size = mpiio::header_size;
  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(nx_) * ny_ * nz_ * sizeof(float);
  const int tile_elems = nxt_ * nyt_ * nzt_;
  const int nf = num_fields_;

  // Pre-allocate the file to its full size so that derived-type writes
  // don't leave the file truncated (some MPI implementations do not
  // extend the file for sparse write patterns).
  MPI_Offset total_size = hdr_size + static_cast<MPI_Offset>(nf) * field_bytes;
  rc = MPI_File_set_size(fh, total_size);
  if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }

  // --- Rebuild cached MPI filetypes if the tile layout changed ---
  auto local_cids = grid.get_local_tiles();
  if (local_cids != cached_cids_)
    rebuild_filetypes_(grid);

  const int ntiles = static_cast<int>(local_cids.size());

  // --- Pack + write per tile using cached filetypes ---
  //
  // Each tile's filetype is a struct of nf subarray copies at
  // field_bytes spacing, matching tile_buf_'s SoA layout.
  // Types are cached across calls and rebuilt only when the
  // local tile set changes.
  //
  // Data is written directly from tile_buf_.span().data():
  //   CPU: host pointer (thrust::host_vector)
  //   HIP: device pointer (thrust::device_vector) — GPU-aware MPI
  //        (OpenMPI 5+) handles staging internally.
  //
  // Because MPI_File_set_view is collective, all ranks must iterate
  // the same number of times (cached_max_ntiles_ across all ranks).

  bool ok = true;
  for (int t = 0; t < cached_max_ntiles_ && ok; t++) {

    // --- pack this tile on device (no D→H copy) ---
    if (t < ntiles) {
      auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile(local_cids[t]));
      pack_tile(tile);
    }

    MPI_Datatype filetype = (t < ntiles)
      ? cached_filetypes_[t]
      : cached_empty_type_;
    int my_count = (t < ntiles) ? nf * tile_elems : 0;

    // --- collective set_view + write directly from device buffer ---
    rc = MPI_File_set_view(
      fh, hdr_size, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    if (rc != MPI_SUCCESS) { ok = false; }

    if (ok) {
      MPI_Status status;
      rc = MPI_File_write_all(
        fh,
        tile_buf_.span().data(),  // CPU: host ptr; HIP: device ptr
        my_count, MPI_FLOAT, &status);
      if (rc != MPI_SUCCESS) ok = false;
    }
  }

  MPI_File_close(&fh);
  return ok;
}


//--------------------------------------------------
// explicit template class instantiation
template class mpiio::FieldsWriter<3>;
