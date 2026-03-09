#include <mpi.h>

#include "io/snapshots/mpiio_fields.h"
#include "io/snapshots/mpiio_header.h"
#include "core/emf/tile.h"

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

  // Field names: 10 entries of 16 chars each, starting at offset 64
  for (int f = 0; f < num_fields; f++) {
    std::strncpy(buf + 64 + f * 16, field_names[f], 15);
  }

  // [224:256] already zeroed

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
    stride_(stride)
{
  nxt_ = std::max(1, NxMesh_ / stride_);
  nyt_ = std::max(1, NyMesh_ / stride_);
  nzt_ = std::max(1, NzMesh_ / stride_);

  nx_ = Nx_ * nxt_;
  ny_ = Ny_ * nyt_;
  nz_ = Nz_ * nzt_;

  tile_buf_.resize(num_fields * nxt_ * nyt_ * nzt_);
}


template<>
void mpiio::FieldsWriter<3>::read_tile(emf::Tile<3>& tile)
{
  const auto [Emds, Bmds, Jmds] = tile.view_EBJ_on_host();

  const int tile_elems = nxt_ * nyt_ * nzt_;

  // E and B fields: subsample by stride-hopping
  for (int ks = 0; ks < nzt_; ks++)
    for (int js = 0; js < nyt_; js++)
      for (int is = 0; is < nxt_; is++) {
        const int idx = ks*nyt_*nxt_ + js*nxt_ + is;

        // ex, ey, ez
        tile_buf_[0*tile_elems + idx] =
          Emds[is*stride_, js*stride_, ks*stride_][0];
        tile_buf_[1*tile_elems + idx] =
          Emds[is*stride_, js*stride_, ks*stride_][1];
        tile_buf_[2*tile_elems + idx] =
          Emds[is*stride_, js*stride_, ks*stride_][2];

        // bx, by, bz
        tile_buf_[3*tile_elems + idx] =
          Bmds[is*stride_, js*stride_, ks*stride_][0];
        tile_buf_[4*tile_elems + idx] =
          Bmds[is*stride_, js*stride_, ks*stride_][1];
        tile_buf_[5*tile_elems + idx] =
          Bmds[is*stride_, js*stride_, ks*stride_][2];
      }

  // J fields: volume-sum over stride^3 cells
  const int jx_off = 6 * tile_elems;
  const int jy_off = 7 * tile_elems;
  const int jz_off = 8 * tile_elems;

  for (int ks = 0; ks < nzt_; ks++)
    for (int js = 0; js < nyt_; js++)
      for (int is = 0; is < nxt_; is++) {
        float sjx = 0.0f, sjy = 0.0f, sjz = 0.0f;
        for (int kstride = 0; kstride < stride_; kstride++)
          for (int jstride = 0; jstride < stride_; jstride++)
            for (int istride = 0; istride < stride_; istride++) {
              sjx += Jmds[is*stride_ + istride,
                          js*stride_ + jstride,
                          ks*stride_ + kstride][0];
              sjy += Jmds[is*stride_ + istride,
                          js*stride_ + jstride,
                          ks*stride_ + kstride][1];
              sjz += Jmds[is*stride_ + istride,
                          js*stride_ + jstride,
                          ks*stride_ + kstride][2];
            }
        const int idx = ks*nyt_*nxt_ + js*nxt_ + is;
        tile_buf_[jx_off + idx] = sjx;
        tile_buf_[jy_off + idx] = sjy;
        tile_buf_[jz_off + idx] = sjz;
      }

  // rho: write zeros (rho not in EBJ, will be added later)
  const int rho_off = 9 * tile_elems;
  std::memset(tile_buf_.data() + rho_off, 0, tile_elems * sizeof(float));
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
    read_tile(tile);

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
            tile_buf_.data() + buf_offset,
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
