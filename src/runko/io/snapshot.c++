// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/io/snapshot.h"

#include "mpi.h"
#include "pika/mpi.hpp"
#include "runko/actions/pic.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/coords.h"
#include "runko/emf/yee_lattice.h"
#include "runko/mdgrid_common.h"
#include "runko/particles_common.h"
#include "runko/tools/config_parser.h"
#include "runko/tools/math.h"
#include "thrust/copy.h"
#include "thrust/count.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/memory.h"
#include "tyvi/execution.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace runko {

namespace ta   = tyvi::actions;
namespace te   = tyvi::exec;
namespace pmpi = pika::mpi::experimental;

namespace {

/// Write the 512-byte binary header to an MPI file handle.
/// Returns MPI error code (MPI_SUCCESS on success).
///
/// Header layout (all little-endian):
///   [0:4]     uint32  magic
///   [4:8]     uint32  version
///   [8:12]    uint32  header_size
///   [12:16]   uint32  num_fields
///   [16:20]   int32   nx  (global output dims after stride)
///   [20:24]   int32   ny
///   [24:28]   int32   nz
///   [28:32]   int32   stride
///   [32:36]   int32   Nx  (tile grid dims)
///   [36:40]   int32   Ny
///   [40:44]   int32   Nz
///   [44:48]   int32   NxMesh
///   [48:52]   int32   NyMesh
///   [52:56]   int32   NzMesh
///   [56:60]   int32   lap
///   [60:64]   uint32  dtype_size = 4
///   [64:64+num_fields*16]  char[16]*num_fields  field names (null-padded)
///   [remainder:512]  reserved zeros
auto
  write_header(
    MPI_File fh,
    int32_t nx,
    int32_t ny,
    int32_t nz,
    int32_t stride,
    int32_t Nx,
    int32_t Ny,
    int32_t Nz,
    int32_t NxMesh,
    int32_t NyMesh,
    int32_t NzMesh,
    int32_t lap,
    int32_t num_fields)
{
  return te::just() | te::then([=] {
           auto buf        = std::array<char, header_size> {};
           auto buf_ptr_at = [&](const auto n) {
             return std::ranges::next(buf.data(), n);
           };

           auto put32 = [&](int offset, std::uint32_t val) {
             std::memcpy(buf_ptr_at(offset), &val, 4);
           };

           auto puti32 = [&](int offset, std::int32_t val) {
             std::memcpy(buf_ptr_at(offset), &val, 4);
           };

           put32(0, magic);
           put32(4, version);
           put32(8, static_cast<std::uint32_t>(header_size));
           put32(12, static_cast<std::uint32_t>(num_fields));

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
           for(int f = 0; f < num_emf_fields; f++) {
             const auto name = std::string_view { emf_field_names[f] };
             const auto n    = std::ranges::min(name.size(), 15uz);
             std::memcpy(buf_ptr_at(64 + f * 16), name.data(), n);
           }
           for(int s = 0; s < num_fields - num_emf_fields; s++) {
             const auto name = std::format("n{}", s);
             const auto n    = std::ranges::min(name.size(), 15uz);
             std::memcpy(buf_ptr_at(64 + (num_emf_fields + s) * 16), name.data(), n);
           }

           MPI_Status status;
           MPI_File_write_at(fh, 0, buf.data(), 512, MPI_BYTE, &status);
           //
           // return te::just(fh, 0, buf.data(), header_size, MPI_BYTE) |
           // pmpi::transform_mpi(MPI_File_iwrite_at);
         });
}

/// Number of particle species in the config (q0/m0, q1/m1, ...).
int
  count_species(const toolbox::ConfigParser& config)
{
  auto nspecies = 0;
  while(config.get<double>(std::format("q{}", nspecies)) and
        config.get<double>(std::format("m{}", nspecies))) {
    ++nspecies;
  }
  return nspecies;
}

/// Snapshot parameters parsed from the config.
struct emf_snapshot_params {
  int Nx, Ny, Nz;
  int NxMesh, NyMesh, NzMesh;
  int stride;
  int nspecies;       // number of particle species (capped at max_species)
  int nxt, nyt, nzt;  // per-tile output size after stride
  int nx, ny, nz;     // global output size

  int num_fields() const { return num_emf_fields + nspecies; }

  static emf_snapshot_params from_config(const toolbox::ConfigParser& config)
  {
    const auto tiles = toolbox::get_extent_list(config, "n_tiles", 3);
    const auto cells = toolbox::get_extent_list(config, "n_cells_per_tile", 3);

    const auto stride = config.get<std::ptrdiff_t>("io_grid_stride").value_or(1);
    if(stride <= 0) {
      throw std::runtime_error {
        std::format("emf_snapshot: io_grid_stride has to be positive (got {})", stride)
      };
    }

    const auto nspecies = std::min(count_species(config), max_species);

    auto p = emf_snapshot_params {};

    p.Nx       = checked_cast<int>(tiles[0]);
    p.Ny       = checked_cast<int>(tiles[1]);
    p.Nz       = checked_cast<int>(tiles[2]);
    p.NxMesh   = checked_cast<int>(cells[0]);
    p.NyMesh   = checked_cast<int>(cells[1]);
    p.NzMesh   = checked_cast<int>(cells[2]);
    p.stride   = checked_cast<int>(stride);
    p.nspecies = nspecies;

    p.nxt = std::max(1, p.NxMesh / p.stride);
    p.nyt = std::max(1, p.NyMesh / p.stride);
    p.nzt = std::max(1, p.NzMesh / p.stride);

    p.nx = p.Nx * p.nxt;
    p.ny = p.Ny * p.nyt;
    p.nz = p.Nz * p.nzt;

    return p;
  }
};

/// Output directory from argument or from config (io_outdir).
std::string
  resolve_outdir(
    const toolbox::ConfigParser& config,
    const std::optional<std::string>& outdir)
{
  if(outdir) { return outdir.value(); }

  const auto x = config.get<std::string>("io_outdir");
  if(not x or x.value().empty()) { return "runko_output"; }
  if(x.value() == "auto") {
    throw std::runtime_error {
      "snapshot: io_outdir = \"auto\" can not be resolved, pass the output "
      "directory explicitly."
    };
  }
  return x.value();
}

/// Pack tile field data into tile_buf using for_each_index kernels.
///
/// Stride subsampling for E and B, volume-sum for J
/// and number density deposited directly at coarse resolution.
void
  pack_tile(
    const emf::YeeLattice& yee,
    const pic::particle_containers* particles,
    const std::array<double, 3> tile_mins,
    const emf_snapshot_params& p,
    runko::IOFieldGrid<float>& tile_buf)
{
  const auto Emds = yee.nonhalo_submds(yee.mds_E());
  const auto Bmds = yee.nonhalo_submds(yee.mds_B());
  const auto Jmds = yee.nonhalo_submds(yee.mds_J());

  const int stride   = p.stride;
  const int nf       = p.num_fields();
  const auto buf_mds = tile_buf.mds();

  tyvi::mdgrid_work w {};

  // E and B fields: subsample by stride-hopping; zero density slots
  w.for_each_index(buf_mds, [=](const auto idx) {
    const auto iz = idx[0], iy = idx[1], ix = idx[2];

    const auto si = ix * static_cast<std::size_t>(stride);
    const auto sj = iy * static_cast<std::size_t>(stride);
    const auto sk = iz * static_cast<std::size_t>(stride);

    // ex, ey, ez
    buf_mds[iz, iy, ix][0] = Emds[si, sj, sk][0];
    buf_mds[iz, iy, ix][1] = Emds[si, sj, sk][1];
    buf_mds[iz, iy, ix][2] = Emds[si, sj, sk][2];

    // bx, by, bz
    buf_mds[iz, iy, ix][3] = Bmds[si, sj, sk][0];
    buf_mds[iz, iy, ix][4] = Bmds[si, sj, sk][1];
    buf_mds[iz, iy, ix][5] = Bmds[si, sj, sk][2];

    // density slots: zero n0..n{nspecies-1} (overwritten below if there are particles)
    for(int s = 0; s < nf - num_emf_fields; s++)
      buf_mds[iz, iy, ix][num_emf_fields + s] = 0.0f;
  });

  // J fields: volume-sum over stride^3 cells
  w.for_each_index(buf_mds, [=](const auto idx) {
    const auto iz = idx[0], iy = idx[1], ix = idx[2];

    const auto ustride = static_cast<std::size_t>(stride);

    float sjx = 0.0f, sjy = 0.0f, sjz = 0.0f;
    for(int kk = 0; kk < stride; kk++)
      for(int jj = 0; jj < stride; jj++)
        for(int ii = 0; ii < stride; ii++) {
          const auto si = ix * ustride + static_cast<std::size_t>(ii);
          const auto sj = iy * ustride + static_cast<std::size_t>(jj);
          const auto sk = iz * ustride + static_cast<std::size_t>(kk);
          sjx += Jmds[si, sj, sk][0];
          sjy += Jmds[si, sj, sk][1];
          sjz += Jmds[si, sj, sk][2];
        }

    buf_mds[iz, iy, ix][6] = sjx;
    buf_mds[iz, iy, ix][7] = sjy;
    buf_mds[iz, iy, ix][8] = sjz;
  });

  if(not particles) {
    w.wait();
    return;
  }

  // Number density: deposit particles directly at coarse resolution.
  // Density slots are already zeroed in the EB kernel above.
  using vt              = pic::ParticleContainer::value_type;
  const auto mx         = static_cast<vt>(tile_mins[0]);
  const auto my         = static_cast<vt>(tile_mins[1]);
  const auto mz         = static_cast<vt>(tile_mins[2]);
  const auto inv_stride = vt { 1 } / static_cast<vt>(stride);
  // particles appended after pack_outgoing may sit exactly on the upper tile face
  const auto last = std::array { static_cast<vt>(p.nxt - 1),
                                 static_cast<vt>(p.nyt - 1),
                                 static_cast<vt>(p.nzt - 1) };

  for(int s = 0; s < p.nspecies; s++) {
    const auto species = static_cast<std::size_t>(s);
    if(not particles->contains(species)) { continue; }

    const auto& container = particles->at(species);
    const auto pos_mds    = container.pos_mds();
    const auto ids_mds    = container.ids_mds();
    const auto fi         = static_cast<std::size_t>(num_emf_fields + s);

    w.for_each_index(pos_mds, [=](const auto idx) {
      if(ids_mds[idx][] == runko::dead_prtc_id) { return; }
      const auto px = pos_mds[idx][0] - mx;
      const auto py = pos_mds[idx][1] - my;
      const auto pz = pos_mds[idx][2] - mz;

      const auto ci =
        static_cast<std::size_t>(sstd::min(sstd::floor(px * inv_stride), last[0]));
      const auto cj =
        static_cast<std::size_t>(sstd::min(sstd::floor(py * inv_stride), last[1]));
      const auto ck =
        static_cast<std::size_t>(sstd::min(sstd::floor(pz * inv_stride), last[2]));

      auto* const n = &thrust::raw_reference_cast(buf_mds[ck, cj, ci][fi]);
      sstd::atomic_add(n, vt { 1 });
    });
  }
  w.wait();
}

/// Write the 512-byte particle snapshot header to an MPI file handle.
///
/// Header layout is the same as for fields snapshot, except:
///   [16:24]   int64   n_prtcls (global number of written particles)
///   [24:28]   int32   unused (zero)
///   [28:32]   int32   species
///   [32:56]   int32   unused (zeros)
///   [64:64+num_prtcl_fields*16]  char[16]*num_prtcl_fields  field names
void
  write_prtcl_header(
    MPI_File fh,
    const std::int64_t n_prtcls,
    const std::int32_t species,
    const std::int32_t lap)
{
  auto buf        = std::array<char, header_size> {};
  auto buf_ptr_at = [&](const auto n) { return std::ranges::next(buf.data(), n); };

  auto put = [&](int offset, const auto val) {
    std::memcpy(buf_ptr_at(offset), &val, sizeof(val));
  };

  put(0, magic);
  put(4, version);
  put(8, static_cast<std::uint32_t>(header_size));
  put(12, static_cast<std::uint32_t>(num_prtcl_fields));
  put(16, n_prtcls);
  put(28, species);
  put(56, lap);
  put(60, std::uint32_t { 4 });  // dtype_size = sizeof(float)

  for(int f = 0; f < num_prtcl_fields; f++) {
    const auto name = std::string_view { prtcl_field_names[f] };
    const auto n    = std::ranges::min(name.size(), 15uz);
    std::memcpy(buf_ptr_at(64 + f * 16), name.data(), n);
  }

  MPI_Status status;
  if(
    MPI_SUCCESS !=
    MPI_File_write_at(fh, 0, buf.data(), header_size, MPI_BYTE, &status)) {
    throw std::runtime_error { "prtcl_snapshot: writing header failed!" };
  }
}

/// Every stride:th particle slot holding a live particle is sampled.
auto
  is_sampled(const pic::ParticleContainer& container, const runko::index_t stride)
{
  return [=, ids_mds = container.ids_mds()](const runko::index_t n) {
    return n % stride == 0u and ids_mds[n][] != runko::dead_prtc_id;
  };
}

std::int64_t
  count_sampled(const pic::ParticleContainer& container, const runko::index_t stride)
{
  const auto slots_begin = thrust::counting_iterator<runko::index_t>(0uz);

  return static_cast<std::int64_t>(thrust::count_if(
    tyvi::mdgrid_work {}.on_this(),
    slots_begin,
    slots_begin + container.ssize(),
    is_sampled(container, stride)));
}

/// Sampling of particles of one species in one local tile.
struct prtcl_tile_sample {
  const pic::ParticleContainer* container;
  const emf::YeeLattice* yee;
  std::array<emf::YeeLattice::value_type, 3> lattice_origo;
  std::int64_t count;         // live particles
  runko::index_t stride = 1;  // every stride:th slot is sampled
  std::int64_t sampled  = 0;  // live particles in sampled slots
};

/// Write sampled particles of given species to {outdir}/prtcls_{species}_{lap}.bin.
///
/// Writes a sub-sample of approximately n_prtcls particles in total,
/// with per-tile quotas proportional to the tile's particle count.
void
  write_prtcl_species(
    runko::simulation_context& sim,
    const long lap,
    const std::size_t species,
    const std::int64_t n_prtcls,
    const std::string& outdir)
{
  using yee_value_type = emf::YeeLattice::value_type;

  // 1. Count live particles in local tiles.
  auto samples             = std::vector<prtcl_tile_sample> {};
  std::int64_t local_total = 0;

  for(auto&& [id, yee, idx]: sim.view_tiles<
                             emf::YeeLattice,
                             runko::cartesian_index<3>,
                             runko::local_tile_tag>()) {
    const auto* particles = sim.tiles.try_get<pic::particle_containers>(id);
    if(not particles or not particles->contains(species)) { continue; }

    const auto& container = particles->at(species);
    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);
    const auto lattice_origo =
      std::array { static_cast<yee_value_type>(gc.mins()[0]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[1]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[2]) - emf::halo_size };

    const auto count = count_sampled(container, 1);
    samples.push_back({ &container, &yee, lattice_origo, count });
    local_total += count;
  }

  // 2. Per-tile quotas and strides based on the global particle count.
  std::int64_t global_total = 0;
  if(
    MPI_SUCCESS != MPI_Allreduce(
                     &local_total,
                     &global_total,
                     1,
                     MPI_INT64_T,
                     MPI_SUM,
                     MPI_COMM_WORLD)) {
    throw std::runtime_error { std::format(
      "prtcl_snapshot: total particle count MPI_Allreduce failed") };
  }

  const auto frac =
    global_total > 0
      ? std::min(1.0, static_cast<double>(n_prtcls) / static_cast<double>(global_total))
      : 0.0;

  std::int64_t local_sampled = 0;
  for(auto& t: samples) {
    const auto quota = std::min(
      t.count,
      static_cast<std::int64_t>(std::round(static_cast<double>(t.count) * frac)));
    if(quota <= 0) { continue; }

    t.stride  = static_cast<runko::index_t>(t.count / quota);
    t.sampled = count_sampled(*t.container, t.stride);
    local_sampled += t.sampled;
  }

  // 3. Location of this rank's particles in the file.
  std::int64_t rank_offset = 0;
  if(
    MPI_SUCCESS !=
    MPI_Exscan(&local_sampled, &rank_offset, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD)) {
    throw std::runtime_error { std::format(
      "prtcl_snapshot: rank offset MPI_Exscan failed") };
  }

  int comm_rank = -1;
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank)) {
    throw std::runtime_error { std::format("prtcl_snapshot: MPI_Comm_rank failed") };
  }

  // MPI_Exscan leaves receive buffer undefined on rank 0.
  if(comm_rank == 0) { rank_offset = 0; }

  std::int64_t global_sampled = 0;

  if(
    MPI_SUCCESS != MPI_Allreduce(
                     &local_sampled,
                     &global_sampled,
                     1,
                     MPI_INT64_T,
                     MPI_SUM,
                     MPI_COMM_WORLD)) {
    throw std::runtime_error { std::format(
      "prtcl_snapshot: gloabl_sampled MPI_Allreduce failed") };
  }

  // 4. Pack sampled particles and their interpolated fields to SoA buffer:
  //    field f is contiguous at offset f * local_sampled.
  auto prtcl_buf =
    runko::PrtclFieldList<float>(static_cast<std::size_t>(local_sampled));
  const auto buf_mds = prtcl_buf.mds();

  auto sampled_slots   = runko::device_scalar_segments<runko::index_t>();
  auto tile_buf_offset = 0uz;

  const auto w = tyvi::mdgrid_work {};
  for(const auto& t: samples) {
    if(t.sampled == 0) { continue; }
    const auto n_s = static_cast<std::size_t>(t.sampled);

    sampled_slots.resize(n_s);
    {
      const auto slots_begin = thrust::counting_iterator<runko::index_t>(0uz);
      const auto out         = sampled_slots.component_view<>();

      const auto sampled_end = thrust::copy_if(
        w.on_this(),
        slots_begin,
        slots_begin + t.container->ssize(),
        out.begin(),
        is_sampled(*t.container, t.stride));

      if(sampled_end != out.end()) {
        throw std::logic_error {
          "prtcl_snapshot: number of sampled slots changed between counting and "
          "copying."
        };
      }
    }

    const auto slots_mds    = sampled_slots.mds();
    const auto pos_mds      = t.container->pos_mds();
    const auto vel_mds      = t.container->vel_mds();
    const auto interpolator = t.yee->interpolate_EB_linear_1st(t.lattice_origo);

    const auto dst_range = std::array { tile_buf_offset, tile_buf_offset + n_s };
    const auto dst_mds   = std::submdspan(buf_mds, dst_range);

    w.for_each_index(dst_mds, [=](const auto idx) {
      using vt     = pic::ParticleContainer::value_type;
      const auto n = static_cast<std::size_t>(slots_mds[idx[0]][]);

      // particle positions are already global
      dst_mds[idx][0] = pos_mds[n][0];
      dst_mds[idx][1] = pos_mds[n][1];
      dst_mds[idx][2] = pos_mds[n][2];
      dst_mds[idx][3] = vel_mds[n][0];
      dst_mds[idx][4] = vel_mds[n][1];
      dst_mds[idx][5] = vel_mds[n][2];

      const auto EB = interpolator(toolbox::Vec3<vt>(pos_mds[n]));

      dst_mds[idx][6]  = EB.E[0];
      dst_mds[idx][7]  = EB.E[1];
      dst_mds[idx][8]  = EB.E[2];
      dst_mds[idx][9]  = EB.B[0];
      dst_mds[idx][10] = EB.B[1];
      dst_mds[idx][11] = EB.B[2];
    });

    tile_buf_offset += n_s;
  }
  w.wait();

  // On CPU the device buffer is host-accessible; on GPU use staging.
#if defined(TYVI_BACKEND_CPU)
  const float* write_ptr = prtcl_buf.span().data();
#elif defined(TYVI_BACKEND_HIP)
  tyvi::mdgrid_work {}.sync_to_staging(prtcl_buf).wait();
  const float* write_ptr = prtcl_buf.staging_span().data();
#endif

  // 5. Write to file.
  const auto filename = std::format("{}/prtcls_{}_{}.bin", outdir, species, lap);

  MPI_File fh;
  if(
    MPI_SUCCESS != MPI_File_open(
                     MPI_COMM_WORLD,
                     filename.c_str(),
                     MPI_MODE_CREATE | MPI_MODE_WRONLY,
                     MPI_INFO_NULL,
                     &fh)) {
    throw std::runtime_error {
      std::format("prtcl_snapshot: could not open {}", filename)
    };
  }

  auto throw_n_close = [&](const std::string_view what) {
    MPI_File_close(&fh);
    throw std::runtime_error {
      std::format("prtcl_snapshot: {} failed for {}", what, filename)
    };
  };

  const MPI_Offset field_bytes =
    static_cast<MPI_Offset>(global_sampled) * static_cast<MPI_Offset>(sizeof(float));

  // Pre-allocate the file to its full size
  const MPI_Offset total_size = header_size + num_prtcl_fields * field_bytes;
  if(MPI_SUCCESS != MPI_File_set_size(fh, total_size)) {
    throw_n_close("MPI_File_set_size");
  }

  if(comm_rank == 0) {
    write_prtcl_header(
      fh,
      global_sampled,
      checked_cast<std::int32_t>(species),
      checked_cast<std::int32_t>(lap));
  }

  if(local_sampled > 0) {
    for(int f = 0; f < num_prtcl_fields; f++) {
      const MPI_Offset file_offset =
        header_size + f * field_bytes +
        static_cast<MPI_Offset>(rank_offset) * static_cast<MPI_Offset>(sizeof(float));

      const auto buf_offset = f * local_sampled;

      MPI_Status status;
      if(
        MPI_SUCCESS != MPI_File_write_at(
                         fh,
                         file_offset,
                         std::ranges::next(write_ptr, buf_offset),
                         checked_cast<int>(local_sampled),
                         MPI_FLOAT,
                         &status)) {
        throw_n_close("MPI_File_write_at");
      }
    }
  }

  if(MPI_SUCCESS != MPI_File_close(&fh)) {
    throw std::runtime_error {
      std::format("prtcl_snapshot: could not close {}", filename)
    };
  }
}


/// Spectra snapshot parameters parsed from the config.
struct spectra_snapshot_params {
  int Nx, Ny, Nz;
  int NxMesh, NyMesh, NzMesh;
  int stride;
  int nbins;
  float umin, umax;
  int nspecies;  // number of particle species (capped at max_spectra_species)
  int nxt;       // per-tile x output size after stride
  int nx;        // global x output size

  int num_fields() const { return nspecies * num_spectra_per_species; }

  static spectra_snapshot_params from_config(const toolbox::ConfigParser& config)
  {
    const auto tiles = toolbox::get_extent_list(config, "n_tiles", 3);
    const auto cells = toolbox::get_extent_list(config, "n_cells_per_tile", 3);

    const auto stride =
      config.get<std::ptrdiff_t>("io_spectra_stride")
        .or_else([&] { return config.get<std::ptrdiff_t>("io_grid_stride"); })
        .value_or(1);
    if(stride <= 0) {
      throw std::runtime_error { std::format(
        "spectra_snapshot: io_spectra_stride has to be positive (got {})",
        stride) };
    }

    const auto nbins = config.get<std::ptrdiff_t>("io_n_spectra_bins").value_or(200);
    if(nbins <= 0) {
      throw std::runtime_error { std::format(
        "spectra_snapshot: io_n_spectra_bins has to be positive (got {})",
        nbins) };
    }

    const auto umin = config.get<double>("io_spectra_umin").value_or(1e-4);
    const auto umax = config.get<double>("io_spectra_umax").value_or(1e3);
    if(not(0.0 < umin and umin < umax)) {
      throw std::runtime_error { std::format(
        "spectra_snapshot: 0 < io_spectra_umin < io_spectra_umax does not hold "
        "(got {} and {})",
        umin,
        umax) };
    }

    auto p = spectra_snapshot_params {};

    p.Nx       = checked_cast<int>(tiles[0]);
    p.Ny       = checked_cast<int>(tiles[1]);
    p.Nz       = checked_cast<int>(tiles[2]);
    p.NxMesh   = checked_cast<int>(cells[0]);
    p.NyMesh   = checked_cast<int>(cells[1]);
    p.NzMesh   = checked_cast<int>(cells[2]);
    p.stride   = checked_cast<int>(stride);
    p.nbins    = checked_cast<int>(nbins);
    p.umin     = static_cast<float>(umin);
    p.umax     = static_cast<float>(umax);
    p.nspecies = std::min(count_species(config), max_spectra_species);

    p.nxt = std::max(1, p.NxMesh / p.stride);
    p.nx  = p.Nx * p.nxt;

    return p;
  }
};

/// Write the 512-byte spectra snapshot header to an MPI file handle.
///
/// Header layout is the same as for fields snapshot (with ny = Ny, nz = Nz), and:
///   [64:64+num_fields*16]  char[16]*num_fields  field names (s0_u, s0_bx, ...)
///   [256:260] int32   nbins
///   [260:264] float   umin
///   [264:268] float   umax
void
  write_spectra_header(
    MPI_File fh,
    const spectra_snapshot_params& p,
    const std::int32_t lap)
{
  auto buf        = std::array<char, header_size> {};
  auto buf_ptr_at = [&](const auto n) { return std::ranges::next(buf.data(), n); };

  auto put = [&](int offset, const auto val) {
    std::memcpy(buf_ptr_at(offset), &val, sizeof(val));
  };

  put(0, magic);
  put(4, version);
  put(8, static_cast<std::uint32_t>(header_size));
  put(12, static_cast<std::uint32_t>(p.num_fields()));
  put(16, std::int32_t { p.nx });
  put(20, std::int32_t { p.Ny });
  put(24, std::int32_t { p.Nz });
  put(28, std::int32_t { p.stride });
  put(32, std::int32_t { p.Nx });
  put(36, std::int32_t { p.Ny });
  put(40, std::int32_t { p.Nz });
  put(44, std::int32_t { p.NxMesh });
  put(48, std::int32_t { p.NyMesh });
  put(52, std::int32_t { p.NzMesh });
  put(56, lap);
  put(60, std::uint32_t { 4 });  // dtype_size = sizeof(float)

  for(int s = 0; s < p.nspecies; s++) {
    for(int q = 0; q < num_spectra_per_species; q++) {
      const auto name = std::format("s{}_{}", s, spectra_suffixes[q]);
      const auto f    = s * num_spectra_per_species + q;
      const auto n    = std::ranges::min(name.size(), 15uz);
      std::memcpy(buf_ptr_at(64 + f * 16), name.data(), n);
    }
  }

  put(256, std::int32_t { p.nbins });
  put(260, p.umin);
  put(264, p.umax);

  MPI_Status status;
  if(
    MPI_SUCCESS !=
    MPI_File_write_at(fh, 0, buf.data(), header_size, MPI_BYTE, &status)) {
    throw std::runtime_error { "spectra_snapshot: writing header failed!" };
  }
}

/// Histogram particles of one tile into tile_buf ([nxt][nbins] with field
/// species * 4 + {0=u, 1=bx, 2=by, 3=bz}).
///
/// The 4 spectra per species are:
///   u      = sqrt(ux^2 + uy^2 + uz^2)  - log10 bins [umin, umax]
///   beta_x = ux / gamma                - linear bins [-1, +1]
///   beta_y = uy / gamma                - linear bins [-1, +1]
///   beta_z = uz / gamma                - linear bins [-1, +1]
/// where gamma = sqrt(1 + u^2). Out-of-range values go to the boundary bins.
void
  histogram_tile(
    const pic::particle_containers* particles,
    const double tile_xmin,
    const spectra_snapshot_params& p,
    runko::SpectraGrid<float>& tile_buf)
{
  const auto buf_mds = tile_buf.mds();

  tyvi::mdgrid_work w {};

  w.for_each_index(buf_mds, [=](const auto idx) {
    for(int f = 0; f < num_spectra_per_species * max_spectra_species; f++) {
      buf_mds[idx][f] = 0.0f;
    }
  });

  if(not particles) {
    w.wait();
    return;
  }

  using vt              = pic::ParticleContainer::value_type;
  const auto mx         = static_cast<vt>(tile_xmin);
  const auto inv_stride = vt { 1 } / static_cast<vt>(p.stride);
  const auto log_umin   = static_cast<vt>(std::log10(p.umin));
  const auto inv_dlog =
    static_cast<vt>(p.nbins) / static_cast<vt>(std::log10(p.umax) - std::log10(p.umin));
  const auto inv_dbeta  = static_cast<vt>(p.nbins) / vt { 2 };
  const auto last_bin_f = static_cast<vt>(p.nbins - 1);
  // particles appended after pack_outgoing may sit exactly on the upper tile face
  const auto last_x_f = static_cast<vt>(p.nxt - 1);

  for(int s = 0; s < p.nspecies; s++) {
    const auto species = static_cast<std::size_t>(s);
    if(not particles->contains(species)) { continue; }

    const auto& container = particles->at(species);
    const auto pos_mds    = container.pos_mds();
    const auto vel_mds    = container.vel_mds();
    const auto ids_mds    = container.ids_mds();
    const auto base       = static_cast<std::size_t>(s * num_spectra_per_species);

    w.for_each_index(pos_mds, [=](const auto idx) {
      if(ids_mds[idx][] == runko::dead_prtc_id) { return; }

      const auto px = pos_mds[idx][0] - mx;
      const auto ix =
        static_cast<std::size_t>(sstd::min(sstd::floor(px * inv_stride), last_x_f));

      const auto ux = vel_mds[idx][0];
      const auto uy = vel_mds[idx][1];
      const auto uz = vel_mds[idx][2];

      const auto u2        = ux * ux + uy * uy + uz * uz;
      const auto u_mag     = sstd::sqrt(u2);
      const auto inv_gamma = vt { 1 } / sstd::sqrt(vt { 1 } + u2);

      const auto bin = [=](const vt raw) {
        return static_cast<std::size_t>(sstd::clamp(raw, vt { 0 }, last_bin_f));
      };

      const auto ib_u  = bin(sstd::floor((sstd::log10(u_mag) - log_umin) * inv_dlog));
      const auto ib_bx = bin(sstd::floor((ux * inv_gamma + vt { 1 }) * inv_dbeta));
      const auto ib_by = bin(sstd::floor((uy * inv_gamma + vt { 1 }) * inv_dbeta));
      const auto ib_bz = bin(sstd::floor((uz * inv_gamma + vt { 1 }) * inv_dbeta));

      const auto deposit = [=](const std::size_t ib, const std::size_t q) {
        auto* const n = &thrust::raw_reference_cast(buf_mds[ix, ib][base + q]);
        sstd::atomic_add(n, 1.0f);
      };

      deposit(ib_u, 0);
      deposit(ib_bx, 1);
      deposit(ib_by, 2);
      deposit(ib_bz, 3);
    });
  }
  w.wait();
}

}  // namespace

tyvi::actions::sexpr_sender
  emf_snapshot(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir_arg)
{
  return te::just() | te::then([&sim, lap, outdir_arg = std::move(outdir_arg)]() {
           const auto p  = emf_snapshot_params::from_config(sim.config);
           auto tile_buf = runko::IOFieldGrid<float>(
             static_cast<std::size_t>(p.nzt),
             static_cast<std::size_t>(p.nyt),
             static_cast<std::size_t>(p.nxt));
           const auto outdir = resolve_outdir(sim.config, outdir_arg);

           const auto nf = p.num_fields();

           std::filesystem::create_directories(outdir);
           const auto filename = std::format("{}/flds_{}.bin", outdir, lap);


           MPI_File fh;
           if(
             MPI_SUCCESS != MPI_File_open(
                              MPI_COMM_WORLD,
                              filename.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                              MPI_INFO_NULL,
                              &fh)) {
             throw std::runtime_error {
               std::format("emf_snapshot: could not open {}", filename)
             };
           }

           auto throw_n_close = [&](const std::string_view what) {
             MPI_File_close(&fh);
             throw std::runtime_error {
               std::format("emf_snapshot: {} failed for {}", what, filename)
             };
           };

           int comm_rank = -1;
           if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank)) {
             throw_n_close("MPI_Comm_rank");
           }


           const MPI_Offset hdr_size    = header_size;
           const MPI_Offset field_bytes = static_cast<MPI_Offset>(p.nx) * p.ny * p.nz *
                                          static_cast<MPI_Offset>(sizeof(float));
           const int tile_elems         = p.nxt * p.nyt * p.nzt;

           // Pre-allocate the file to its full size
           const MPI_Offset total_size = hdr_size + nf * field_bytes;
           if(MPI_SUCCESS != MPI_File_set_size(fh, total_size)) {
             throw_n_close("MPI_File_set_size");
           }

           // Rank 0 writes the 512-byte header.
           if(comm_rank == 0) {
             tyvi::this_thread::sync_wait(write_header(
               fh,
               p.nx,
               p.ny,
               p.nz,
               p.stride,
               p.Nx,
               p.Ny,
               p.Nz,
               p.NxMesh,
               p.NyMesh,
               p.NzMesh,
               checked_cast<std::int32_t>(lap),
               nf));
           }


           for(auto&& [id, yee, idx]: sim.view_tiles<
                                      emf::YeeLattice,
                                      runko::cartesian_index<3>,
                                      runko::local_tile_tag>()) {
             const auto* particles = sim.tiles.try_get<pic::particle_containers>(id);
             const auto tile_mins =
               runko::global_coordinates(sim, idx.template as<double>().data).mins();

             pack_tile(yee, particles, tile_mins, p, tile_buf);

      // On CPU the device buffer is host-accessible; on GPU use staging.
#if defined(TYVI_BACKEND_CPU)
             const float* write_ptr = tile_buf.span().data();
#elif defined(TYVI_BACKEND_HIP)
             // GPU: copy packed data from device to host staging buffer
             tyvi::mdgrid_work {}.sync_to_staging(tile_buf).wait();
             const float* write_ptr = tile_buf.staging_span().data();
#endif

             const auto ti = static_cast<MPI_Offset>(idx[0]);
             const auto tj = static_cast<MPI_Offset>(idx[1]);
             const auto tk = static_cast<MPI_Offset>(idx[2]);

             // Write each field for this tile
             for(int f = 0; f < nf; f++) {
               const MPI_Offset field_base = hdr_size + f * field_bytes;
               // Write row by row (each row of nxt floats is contiguous in file)
               for(int ks = 0; ks < p.nzt; ks++) {
                 for(int js = 0; js < p.nyt; js++) {
                   const MPI_Offset file_offset =
                     field_base + ((tk * p.nzt + ks) * p.ny * p.nx +
                                   (tj * p.nyt + js) * p.nx + ti * p.nxt) *
                                    static_cast<MPI_Offset>(sizeof(float));

                   const auto buf_offset =
                     f * tile_elems + ks * p.nyt * p.nxt + js * p.nxt;


                   MPI_Status status;
                   const auto rc = MPI_File_write_at(
                     fh,
                     file_offset,
                     write_ptr + buf_offset,
                     p.nxt,
                     MPI_FLOAT,
                     &status);
                   if(rc != MPI_SUCCESS) {
                     throw std::runtime_error {
                       "emf_snapshot: MPI_File_write_at failed!"
                     };
                   }

                   //
                   // senders.push_back(
                   // te::just(
                   // fh,
                   // file_offset,
                   // std::ranges::next(write_ptr, buf_offset),
                   // 1,
                   // MPI_FLOAT) |
                   // pmpi::transform_mpi(MPI_File_iwrite_at));
                 }
               }
             }
           }
           if(MPI_SUCCESS != MPI_File_close(&fh)) {
             throw std::runtime_error {
               std::format("emf_snapshot: could not close {}", filename)
             };
           }
           return ta::null;


           //            return te::when_all_vector(std::move(senders)) |
           //                   te::then([fh, filename = std::move(filename)]() mutable
           //                   {
           //                     if(MPI_SUCCESS != MPI_File_close(&fh)) {
           //                       throw std::runtime_error {
           //                         std::format("emf_snapshot: could not close {}",
           //                         filename)
           //                       };
           //                     }
           //
           //                     return ta::null;
           //                   });
         });
}

tyvi::actions::sexpr_sender
  prtcl_snapshot(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir_arg,
    std::optional<long> n_prtcls_arg)
{
  return te::just() |
         te::then([&sim, lap, outdir_arg = std::move(outdir_arg), n_prtcls_arg]() {
           const auto n_prtcls = n_prtcls_arg.or_else(
             [&] { return sim.config.get<long>("io_n_sampled_prtcls"); });

           if(not n_prtcls or n_prtcls.value() < 0) {
             throw std::runtime_error {
               "prtcl_snapshot: number of sampled particles has to be given "
               "(io_n_sampled_prtcls) and non-negative."
             };
           }

           const auto outdir = resolve_outdir(sim.config, outdir_arg);
           std::filesystem::create_directories(outdir);

           const auto nspecies = static_cast<std::size_t>(count_species(sim.config));
           for(auto s = 0uz; s < nspecies; ++s) {
             write_prtcl_species(sim, lap, s, n_prtcls.value(), outdir);
           }

           return ta::null;
         });
}

tyvi::actions::sexpr_sender
  spectra_snapshot(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir_arg)
{
  return te::just() | te::then([&sim, lap, outdir_arg = std::move(outdir_arg)]() {
           const auto p      = spectra_snapshot_params::from_config(sim.config);
           const auto outdir = resolve_outdir(sim.config, outdir_arg);
           const auto nf     = p.num_fields();

           std::filesystem::create_directories(outdir);
           const auto filename = std::format("{}/pspectra_{}.bin", outdir, lap);

           auto tile_buf = runko::SpectraGrid<float>(
             static_cast<std::size_t>(p.nxt),
             static_cast<std::size_t>(p.nbins));

           MPI_File fh;
           if(
             MPI_SUCCESS != MPI_File_open(
                              MPI_COMM_WORLD,
                              filename.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                              MPI_INFO_NULL,
                              &fh)) {
             throw std::runtime_error {
               std::format("spectra_snapshot: could not open {}", filename)
             };
           }

           auto throw_n_close = [&](const std::string_view what) {
             MPI_File_close(&fh);
             throw std::runtime_error {
               std::format("spectra_snapshot: {} failed for {}", what, filename)
             };
           };

           int comm_rank = -1;
           if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank)) {
             throw_n_close("MPI_Comm_rank");
           }

           // Each field is a (Nz, Ny, nx, nbins) array.
           const MPI_Offset field_bytes = static_cast<MPI_Offset>(p.Nz) * p.Ny * p.nx *
                                          p.nbins *
                                          static_cast<MPI_Offset>(sizeof(float));
           const int tile_row_elems     = p.nxt * p.nbins;

           // Pre-allocate the file to its full size
           const MPI_Offset total_size = header_size + nf * field_bytes;
           if(MPI_SUCCESS != MPI_File_set_size(fh, total_size)) {
             throw_n_close("MPI_File_set_size");
           }

           if(comm_rank == 0) {
             write_spectra_header(fh, p, checked_cast<std::int32_t>(lap));
           }

           for(auto&& [id, idx]:
               sim.view_tiles<runko::cartesian_index<3>, runko::local_tile_tag>()) {
             const auto* particles = sim.tiles.try_get<pic::particle_containers>(id);
             const auto tile_xmin =
               runko::global_coordinates(sim, idx.template as<double>().data).mins()[0];

             histogram_tile(particles, tile_xmin, p, tile_buf);

      // On CPU the device buffer is host-accessible; on GPU use staging.
#if defined(TYVI_BACKEND_CPU)
             const float* write_ptr = tile_buf.span().data();
#elif defined(TYVI_BACKEND_HIP)
             tyvi::mdgrid_work {}.sync_to_staging(tile_buf).wait();
             const float* write_ptr = tile_buf.staging_span().data();
#endif

             const auto ti = static_cast<MPI_Offset>(idx[0]);
             const auto tj = static_cast<MPI_Offset>(idx[1]);
             const auto tk = static_cast<MPI_Offset>(idx[2]);

             // One contiguous block of nxt * nbins floats per tile per field.
             for(int f = 0; f < nf; f++) {
               const MPI_Offset file_offset = header_size + f * field_bytes +
                                              ((tk * p.Ny + tj) * p.nx + ti * p.nxt) *
                                                p.nbins *
                                                static_cast<MPI_Offset>(sizeof(float));

               const auto buf_offset = f * tile_row_elems;

               MPI_Status status;
               if(
                 MPI_SUCCESS != MPI_File_write_at(
                                  fh,
                                  file_offset,
                                  std::ranges::next(write_ptr, buf_offset),
                                  tile_row_elems,
                                  MPI_FLOAT,
                                  &status)) {
                 throw_n_close("MPI_File_write_at");
               }
             }
           }

           if(MPI_SUCCESS != MPI_File_close(&fh)) {
             throw std::runtime_error {
               std::format("spectra_snapshot: could not close {}", filename)
             };
           }

           return ta::null;
         });
}

}  // namespace runko
