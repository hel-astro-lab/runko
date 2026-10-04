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
#include "thrust/memory.h"
#include "tyvi/execution.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <array>
#include <bit>
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

    // Count particle species from config (q0/m0, q1/m1, ...)
    auto nspecies = 0;
    while(nspecies < max_species and
          config.get<double>(std::format("q{}", nspecies)) and
          config.get<double>(std::format("m{}", nspecies))) {
      ++nspecies;
    }

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
      "emf_snapshot: io_outdir = \"auto\" can not be resolved, pass the output "
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
  prtcl_snapshot(runko::simulation_context&, long)
{ return te::just(ta::null); }
tyvi::actions::sexpr_sender
  spectra_snapshot(runko::simulation_context&, long)
{ return te::just(ta::null); }

}  // namespace runko
