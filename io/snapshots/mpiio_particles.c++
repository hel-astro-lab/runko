#include <mpi.h>

#include "io/snapshots/mpiio_particles.h"
#include "io/snapshots/mpiio_header.h"
#include "core/pic/tile.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>


// Particle field names and count are defined in mpiio_header.h:
//   mpiio::prtcl_field_names[12], mpiio::num_prtcl_fields


/// Write a 512-byte particle snapshot header.
///
/// Reuses the same binary layout as write_header but with particle-specific
/// metadata and field names.
static int write_prtcl_header(
  MPI_File fh,
  int64_t n_prtcls_global,
  int32_t species,
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
  auto puti64 = [&](int offset, int64_t val) {
    std::memcpy(buf + offset, &val, 8);
  };

  put32(0, mpiio::magic);
  put32(4, mpiio::version);
  put32(8, static_cast<uint32_t>(mpiio::header_size));
  put32(12, static_cast<uint32_t>(num_fields));

  // n_prtcls_global stored as int64 at offset 16 (occupies nx+ny slots)
  puti64(16, n_prtcls_global);
  puti32(24, 0);           // nz (unused)
  puti32(28, species);     // stride slot repurposed as species
  puti32(32, 0);           // Nx (unused)
  puti32(36, 0);           // Ny (unused)
  puti32(40, 0);           // Nz (unused)
  puti32(44, 0);           // NxMesh (unused)
  puti32(48, 0);           // NyMesh (unused)
  puti32(52, 0);           // NzMesh (unused)
  puti32(56, lap);
  put32(60, 4);            // dtype_size = sizeof(float)

  // Field names: 12 entries of 16 chars each, starting at offset 64
  for (int f = 0; f < num_fields; f++) {
    std::strncpy(buf + 64 + f * 16, mpiio::prtcl_field_names[f], 15);
  }

  MPI_Status status;
  return MPI_File_write_at(fh, 0, buf, 512, MPI_BYTE, &status);
}


//--------------------------------------------------

template<>
mpiio::ParticlesWriter<3>::ParticlesWriter(
  const std::string& prefix,
  int64_t n_prtcls,
  int species)
  : prefix_(prefix),
    n_prtcls_(n_prtcls),
    species_(species),
    prtcl_buf_(0)
{}


//--------------------------------------------------

/// Per-tile sampling info computed during write().
struct TileInfo {
  uint64_t cid;
  int64_t count;
  int64_t stride;
  int64_t sampled;
};


template<>
bool mpiio::ParticlesWriter<3>::write(corgi::Grid<3>& grid, int lap)
{
  const auto sp = static_cast<runko::index_t>(species_);

  // 1. Count particles per local tile
  std::vector<TileInfo> tiles;
  int64_t local_total = 0;

  for (auto cid : grid.get_local_tiles()) {
    auto* pic_tile = dynamic_cast<pic::Tile<3>*>(&grid.get_tile(cid));
    if (!pic_tile) continue;
    if (sp >= pic_tile->number_of_species()) continue;

    const auto count = static_cast<int64_t>(pic_tile->number_of_particles(sp));
    tiles.push_back({ cid, count, 1, 0 });
    local_total += count;
  }

  // 2. Global total particle count
  int64_t global_total = 0;
  MPI_Allreduce(
    &local_total, &global_total, 1, MPI_INT64_T, MPI_SUM,
    MPI_Comm(grid.comm));

  // 3. Per-tile quotas and strides
  const double frac = (global_total > 0)
    ? std::min(1.0, static_cast<double>(n_prtcls_) / static_cast<double>(global_total))
    : 0.0;

  for (auto& t : tiles) {
    if (t.count == 0) continue;
    auto quota = static_cast<int64_t>(std::round(static_cast<double>(t.count) * frac));
    if (quota <= 0)    { t.sampled = 0; continue; }
    if (quota > t.count) quota = t.count;
    t.stride  = t.count / quota;
    t.sampled = t.count / t.stride;
  }

  // 4-6. Compute local/global sampled counts and rank offset
  int64_t local_sampled = 0;
  for (const auto& t : tiles) local_sampled += t.sampled;

  int64_t rank_offset = 0;
  MPI_Exscan(
    &local_sampled, &rank_offset, 1, MPI_INT64_T, MPI_SUM,
    MPI_Comm(grid.comm));
  if (grid.comm.rank() == 0) rank_offset = 0;

  int64_t global_sampled = 0;
  MPI_Allreduce(
    &local_sampled, &global_sampled, 1, MPI_INT64_T, MPI_SUM,
    MPI_Comm(grid.comm));

  // 7. Open file and write header
  const std::string filename =
    prefix_ + "/prtcls_" + std::to_string(species_)
    + "_" + std::to_string(lap) + ".bin";

  MPI_File fh;
  int rc = MPI_File_open(
    MPI_Comm(grid.comm), filename.c_str(),
    MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &fh);
  if (rc != MPI_SUCCESS) return false;

  // Pre-allocate the file to its full size so that independent writes
  // don't leave a sparse file (some MPI implementations handle this poorly).
  const MPI_Offset total_size = mpiio::header_size
    + static_cast<MPI_Offset>(num_prtcl_fields_) * global_sampled
      * static_cast<MPI_Offset>(sizeof(float));
  MPI_File_set_size(fh, total_size);

  if (grid.comm.rank() == 0) {
    rc = write_prtcl_header(
      fh, global_sampled, static_cast<int32_t>(species_), lap, num_prtcl_fields_);
    if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
  }

  // 8. Pack all local sampled particles into prtcl_buf_
  const auto n_local = static_cast<runko::index_t>(local_sampled);
  prtcl_buf_.invalidating_resize(n_local);

  if (local_sampled > 0) {
    const auto buf_mds = prtcl_buf_.mds();
    runko::index_t tile_buf_offset = 0;

    for (const auto& t : tiles) {
      if (t.sampled == 0) continue;

      auto& tile = dynamic_cast<pic::Tile<3>&>(grid.get_tile(t.cid));
      const auto& container = tile.particles(sp);

      const auto n_s      = static_cast<runko::index_t>(t.sampled);
      const auto stride_t = static_cast<runko::index_t>(t.stride);

      // Gather sampled local positions for field interpolation
      auto sampled_pos = runko::VecList<float>(n_s);
      {
        const auto sp_mds  = sampled_pos.mds();
        const auto src_pos = container.pos_mds();

        tyvi::mdgrid_work {}
          .for_each_index(sampled_pos, [=](const auto idx) {
            const auto si = static_cast<runko::index_t>(idx[0] * stride_t);
            sp_mds[idx][0] = src_pos[si][0];
            sp_mds[idx][1] = src_pos[si][1];
            sp_mds[idx][2] = src_pos[si][2];
          })
          .wait();
      }

      // Interpolate E, B at sampled positions
      auto [E_interp, B_interp] = tile.interpolate_fields_at(sampled_pos);

      // Pack all 12 fields into prtcl_buf_ 
      using vt = pic::ParticleContainer::value_type;
      const auto mx = static_cast<vt>(tile.mins[0]);
      const auto my = static_cast<vt>(tile.mins[1]);
      const auto mz = static_cast<vt>(tile.mins[2]);

      const auto src_pos = container.pos_mds();
      const auto src_vel = container.vel_mds();
      const auto E_mds   = E_interp.mds();
      const auto B_mds   = B_interp.mds();

      const auto dst_range = std::array<runko::index_t, 2>{
        tile_buf_offset, tile_buf_offset + n_s };
      const auto dst_sub = std::submdspan(buf_mds, dst_range);

      tyvi::mdgrid_work {}
        .for_each_index(dst_sub, [=](const auto idx) {
          const auto si = static_cast<runko::index_t>(idx[0] * stride_t);

          dst_sub[idx][0]  = src_pos[si][0] + mx;   // global x
          dst_sub[idx][1]  = src_pos[si][1] + my;   // global y
          dst_sub[idx][2]  = src_pos[si][2] + mz;   // global z
          dst_sub[idx][3]  = src_vel[si][0];        // ux
          dst_sub[idx][4]  = src_vel[si][1];        // uy
          dst_sub[idx][5]  = src_vel[si][2];        // uz
          dst_sub[idx][6]  = E_mds[idx][0];         // ex
          dst_sub[idx][7]  = E_mds[idx][1];         // ey
          dst_sub[idx][8]  = E_mds[idx][2];         // ez
          dst_sub[idx][9]  = B_mds[idx][0];         // bx
          dst_sub[idx][10] = B_mds[idx][1];         // by
          dst_sub[idx][11] = B_mds[idx][2];         // bz
        })
        .wait();

      tile_buf_offset += n_s;
    }
  }

  // 9. Sync to staging if GPU
#if defined(TYVI_BACKEND_CPU)
  const float* write_ptr = (local_sampled > 0) ? prtcl_buf_.span().data() : nullptr;
#elif defined(TYVI_BACKEND_HIP)
  if (local_sampled > 0) {
    tyvi::mdgrid_work {}.sync_to_staging(prtcl_buf_).wait();
  }
  const float* write_ptr = (local_sampled > 0) ? prtcl_buf_.staging_span().data() : nullptr;
#endif

  // 10. Write 12 fields via independent MPI_File_write_at
  const MPI_Offset hdr_size    = mpiio::header_size;
  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(global_sampled) * static_cast<MPI_Offset>(sizeof(float));

  float dummy = 0.0f;
  for (int f = 0; f < num_prtcl_fields_; f++) {
    const MPI_Offset file_offset =
      hdr_size + static_cast<MPI_Offset>(f) * field_bytes
      + static_cast<MPI_Offset>(rank_offset) * static_cast<MPI_Offset>(sizeof(float));

    const auto buf_offset = static_cast<int64_t>(f) * local_sampled;

    MPI_Status status;
    rc = MPI_File_write_at(
      fh, file_offset,
      (local_sampled > 0) ? write_ptr + buf_offset : &dummy,
      static_cast<int>(local_sampled), MPI_FLOAT, &status);
    if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
  }

  // 11. Close
  rc = MPI_File_close(&fh);
  return rc == MPI_SUCCESS;
}
