#include <mpi.h>

#include "io/snapshots/mpiio_fields.h"
#include "io/snapshots/mpiio_header.h"
#include "core/emf/tile.h"
#include "core/pic/tile.h"
#include "tools/simd_math.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <cstring>
#include <format>
#include <string>


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
  : prefix_(prefix),
    Nx_(Nx), Ny_(Ny), Nz_(Nz),
    NxMesh_(NxMesh), NyMesh_(NyMesh), NzMesh_(NzMesh),
    stride_(stride),
    nspecies_(std::min(nspecies, static_cast<int>(max_species))),
    num_fields_(num_emf_fields + nspecies_),
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
  const int nf = num_fields_;
  auto buf_mds = tile_buf_.mds();

  tyvi::mdgrid_work w {};

  // E and B fields: subsample by stride-hopping; zero density slots
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
  // Density slots are already zeroed in the EB kernel above.
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

    const auto nspec = static_cast<int>(pic_tile->number_of_species());
    const int ndeposit = std::min(nspec, nspecies_);
    for (int s = 0; s < ndeposit; s++)
      deposit_species(static_cast<runko::size_t>(s), 9 + s);
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

  // Rank 0 writes the 512-byte header.
  // MPI_File_open is collective and synchronizing; no barrier needed
  // since all data writes target non-overlapping offsets >= header_size.
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
          if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
        }
      }
    }
  }

  rc = MPI_File_close(&fh);
  return rc == MPI_SUCCESS;
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

  // --- Collect local tiles ---
  auto local_cids = grid.get_local_tiles();
  const int ntiles = static_cast<int>(local_cids.size());

  struct TileIdx { int ti, tj, tk; };
  std::vector<TileIdx> tile_indices(ntiles);

  // --- Phase 1: Pack all tiles into staging buffer ---
  // Layout: staging[t * nf * tile_elems + f * tile_elems .. +tile_elems]
  // Grouped by tile, then by field — matches the per-tile all-fields
  // struct type's etype slot ordering.
  const int tile_stride = nf * tile_elems;
  write_staging_.resize(static_cast<size_t>(ntiles) * tile_stride);

  for (int t = 0; t < ntiles; t++) {
    auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile(local_cids[t]));
    pack_tile(tile);

#if defined(TYVI_BACKEND_CPU)
    const float* write_ptr = tile_buf_.span().data();
#elif defined(TYVI_BACKEND_HIP)
    const float* write_ptr = tile_buf_.staging_span().data();
#endif

    tile_indices[t].ti = static_cast<int>(tile.index[0]);
    tile_indices[t].tj = static_cast<int>(tile.index[1]);
    tile_indices[t].tk = static_cast<int>(tile.index[2]);

    for (int f = 0; f < nf; f++) {
      std::memcpy(
        write_staging_.data()
          + static_cast<size_t>(t) * tile_stride
          + static_cast<size_t>(f) * tile_elems,
        write_ptr + f * tile_elems,
        static_cast<size_t>(tile_elems) * sizeof(float));
    }
  }

  // --- Phase 2: Build per-tile MPI filetypes and write ---
  //
  // For each tile we create a subarray type describing its footprint
  // in one (nz_,ny_,nx_) field, then wrap nf copies of that subarray
  // in a struct with field_bytes spacing.  This gives a single type
  // covering all fields for one tile, with etype slots ordered as
  // [field0 data][field1 data]...[fieldN data] — matching the staging
  // buffer layout.
  //
  // Because MPI_File_set_view is collective, all ranks must iterate
  // the same number of times (max_ntiles across all ranks).

  int max_ntiles = 0;
  MPI_Allreduce(&ntiles, &max_ntiles, 1, MPI_INT, MPI_MAX,
                MPI_Comm(grid.comm));

  int gsizes[3]   = {nz_, ny_, nx_};
  int subsizes[3] = {nzt_, nyt_, nxt_};

  // Field displacement table and type vector (same shape for every tile)
  std::vector<MPI_Aint>     field_disps(nf);
  std::vector<int>          field_blens(nf, 1);
  std::vector<MPI_Datatype> ftypes(nf);
  for (int f = 0; f < nf; f++)
    field_disps[f] = static_cast<MPI_Aint>(f) * field_bytes;

  bool ok = true;
  for (int t = 0; t < max_ntiles && ok; t++) {

    // --- build filetype for this tile (or empty for non-participating) ---
    MPI_Datatype tile_sub = MPI_DATATYPE_NULL;
    MPI_Datatype filetype = MPI_DATATYPE_NULL;
    int my_count = 0;

    if (t < ntiles) {
      int starts[3] = {
        tile_indices[t].tk * nzt_,
        tile_indices[t].tj * nyt_,
        tile_indices[t].ti * nxt_
      };
      MPI_Type_create_subarray(
        3, gsizes, subsizes, starts,
        MPI_ORDER_C, MPI_FLOAT, &tile_sub);
      MPI_Type_commit(&tile_sub);

      // Struct of nf copies of tile_sub at field_bytes spacing
      std::fill(ftypes.begin(), ftypes.end(), tile_sub);
      MPI_Type_create_struct(
        nf, field_blens.data(), field_disps.data(),
        ftypes.data(), &filetype);
      MPI_Type_commit(&filetype);

      my_count = nf * tile_elems;
    } else {
      MPI_Type_contiguous(0, MPI_FLOAT, &filetype);
      MPI_Type_commit(&filetype);
    }

    // --- collective set_view + write ---
    rc = MPI_File_set_view(
      fh, hdr_size, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    if (rc != MPI_SUCCESS) { ok = false; }

    if (ok) {
      MPI_Status status;
      rc = MPI_File_write_all(
        fh,
        (t < ntiles)
          ? write_staging_.data() + static_cast<size_t>(t) * tile_stride
          : write_staging_.data(),
        my_count, MPI_FLOAT, &status);
      if (rc != MPI_SUCCESS) ok = false;
    }

    // --- free per-tile types ---
    MPI_Type_free(&filetype);
    if (tile_sub != MPI_DATATYPE_NULL) MPI_Type_free(&tile_sub);
  }

  MPI_File_close(&fh);
  return ok;
}


//--------------------------------------------------
// explicit template class instantiation
template class mpiio::FieldsWriter<3>;
