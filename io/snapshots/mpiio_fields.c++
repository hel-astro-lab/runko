#include <mpi.h>

#include "io/snapshots/mpiio_fields.h"
#include "io/snapshots/mpiio_header.h"
#include "core/emf/tile.h"

#include <algorithm>
#include <cstring>
#include <string>


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

  // J fields: volume-average over stride^3 cells
  // Clear jx, jy, jz portions first
  const int jx_off = 6 * tile_elems;
  const int jy_off = 7 * tile_elems;
  const int jz_off = 8 * tile_elems;
  std::memset(tile_buf_.data() + jx_off, 0, 3 * tile_elems * sizeof(float));

  for (int ks = 0; ks < nzt_; ks++)
    for (int kstride = 0; kstride < stride_; kstride++)
      for (int js = 0; js < nyt_; js++)
        for (int jstride = 0; jstride < stride_; jstride++)
          for (int is = 0; is < nxt_; is++)
            for (int istride = 0; istride < stride_; istride++) {
              const int idx = ks*nyt_*nxt_ + js*nxt_ + is;
              tile_buf_[jx_off + idx] +=
                Jmds[is*stride_ + istride,
                     js*stride_ + jstride,
                     ks*stride_ + kstride][0];
              tile_buf_[jy_off + idx] +=
                Jmds[is*stride_ + istride,
                     js*stride_ + jstride,
                     ks*stride_ + kstride][1];
              tile_buf_[jz_off + idx] +=
                Jmds[is*stride_ + istride,
                     js*stride_ + jstride,
                     ks*stride_ + kstride][2];
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
  MPI_File_open(
    MPI_Comm(grid.comm), filename.c_str(),
    MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &fh);

  // Rank 0 writes the 256-byte header
  if (grid.comm.rank() == 0) {
    mpiio::write_header(
      fh,
      nx_, ny_, nz_,
      stride_,
      Nx_, Ny_, Nz_,
      NxMesh_, NyMesh_, NzMesh_,
      lap);
  }

  // Ensure header is written before any rank writes field data
  MPI_Barrier(MPI_Comm(grid.comm));

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
          MPI_File_write_at(
            fh, file_offset,
            tile_buf_.data() + buf_offset,
            nxt_, MPI_FLOAT, &status);
        }
      }
    }
  }

  MPI_File_close(&fh);
  return true;
}


//--------------------------------------------------
// explicit template class instantiation
template class mpiio::FieldsWriter<3>;
