#include "test_prtcl_writer.h"

#include <mpi4cpp/mpi.h>
#include <mpi.h>
#include "../tools/mesh.h"
#include "../tools/ezh5/src/ezh5.hpp"
#include "namer.h"
#include "../tools/fastlog.h"


using namespace mpi4cpp;


using ezh5::File;


template<size_t D>
h5io::TestPrtclWriter<D>::TestPrtclWriter(
        const std::string& prefix, 
        int Nx, int NxMesh,
        int Ny, int NyMesh,
        int Nz, int NzMesh,
        int ppc, int n_local_tiles,
        int n_test_particles_approx) : fname(prefix)
{
  //fname = prefix + "-" + to_string(lap) + extension;

  // total number of particles
  int n_particles = Nx*NxMesh*Ny*NyMesh*Nz*NzMesh*ppc;

  // corresponding thinning factor is then (approximately)
  int stride = n_particles/n_test_particles_approx;

  // how many particles would this particular rank have
  int my_n_prtcls = (n_local_tiles*NxMesh*NyMesh*NzMesh*ppc)/stride;

  // global minimum number of particles per rank
  int min_n_prtcls;

  int n_comm_size, my_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &n_comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Allreduce(&my_n_prtcls, &min_n_prtcls, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  // cutoff id is then (to satisfy min_n_prtcls in our search)
  int cutoff_id = stride*min_n_prtcls;

  // calculate true number of test particles
  int n_test_particles = min_n_prtcls*n_comm_size;


  // get maximum global prtcl number per tile
  //int n_prtcls_per_rank = n_tiles * NxMesh*NyMesh*NzMesh*ppc;

  //n_comm_size
  //n_test_particles = n_comm_size*(n_prtcls_per_rank/stride)
  //stride = n_comm_size*n_prtcls_per_rank/n_test_particles_approx;
  //n_prtcls_per_rank;

  if(my_rank == 0) {
    std::cout << 
         " n_particles: " << n_particles
      << " stride: " << stride
      << " my_n_prtcls: " << my_n_prtcls
      << " min_n_prtcls: " << min_n_prtcls
      << " cutoff_id: " << cutoff_id
      << " n_test_particles: " << n_test_particles
      << "\n";
  }

  /// how many time steps to save
  int nt=1;

  /// total number of (test) particles
  int np = min_n_prtcls;

  /// how many ranks/categories of particles (second id parameter)
  int nr = n_comm_size;

  // 8 variables per particle to save
  for(size_t i=0; i<8; i++) arrs.emplace_back(nt, np, nr);
  for(size_t i=0; i<8; i++) rbuf.emplace_back(nt, np, nr);

  // +2 int arrays
  for(size_t i=0; i<2; i++) arrs2.emplace_back(nt, np, nr);
  for(size_t i=0; i<2; i++) rbuf2.emplace_back(nt, np, nr);
}




template<>
inline void h5io::TestPrtclWriter<2>::read_tiles(
    corgi::Node<2>& grid)
{
  // clear target arrays
  for(auto& arr : arrs) arr.clear();

  // target arrays
  auto& xloc = arrs[0];
  auto& yloc = arrs[1];
  auto& zloc = arrs[2];
  auto& ux   = arrs[3];
  auto& uy   = arrs[4];
  auto& uz   = arrs[5];
  auto& wgt  = arrs[6];

  auto& ids   = arrs2[0];
  auto& procs = arrs2[1];

  int ip;
  int tstep = 1;

  // read my local tiles
  for(auto cid : grid.get_local_tiles() ){
    auto& tile = dynamic_cast<pic::Tile<2>&>(grid.get_tile( cid ));
    auto& container = tile.get_container( ispc );

    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );

    double* ch;
    ch = &( container.wgt(0) );

    // reference the ids 
    int* idn[2];
    for(int i=0; i<2; i++) idn[i] = &( container.id(i,0) );


    // loop and search over all particles
    int ir;
    for(size_t n=0; n<container.size(); n++) {
      if(idn[0][n] % stride == 0) {
        ip = idn[0][n]; // get id
        ir = idn[1][n]; // get proc

        // save particle to correct position
        xloc( tstep, ip, ir) = loc[0][n];
        yloc( tstep, ip, ir) = loc[1][n];
        zloc( tstep, ip, ir) = loc[2][n];
        ux(   tstep, ip, ir) = vel[0][n];
        uy(   tstep, ip, ir) = vel[1][n];
        uz(   tstep, ip, ir) = vel[2][n];
        wgt(  tstep, ip, ir) = ch[n];
        ids(  tstep, ip, ir) = idn[0][n];
        procs(tstep, ip, ir) = idn[1][n];
      }
    }


  }
}



template<size_t D>
inline void h5io::TestPrtclWriter<D>::mpi_reduce_snapshots(
    corgi::Node<D>& grid)
{
  /* based on https://gist.github.com/rmcgibbo/7178576
   */

  for(auto& arr : rbuf) arr.clear();

  int tag = 0;
  const int size = grid.comm.size(); //MPI::COMM_WORLD.Get_size();
  const int rank = grid.comm.rank(); //MPI::COMM_WORLD.Get_rank();
  const int lastpower = 1 << fastlog2(size);

  // each of the ranks greater than the last power of 2 less than size
  // need to downshift their data, since the binary tree reduction below
  // only works when N is a power of two.
  
  std::vector<mpi::request> reqs;
  for (int i = lastpower; i < size; i++) {
    if (rank == i) {
      for(size_t els=0; els<arrs.size(); els++) {
        reqs.push_back( grid.comm.isend(i-lastpower, tag+els, arrs[els].data(), arrs[els].size()) );
      }
    }
  }

  for (int i = 0; i < size-lastpower; i++) {
    if (rank == i) {
      for(size_t els=0; els<arrs.size(); els++) {
        grid.comm.recv(i+lastpower, tag+els, rbuf[els].data(), rbuf[els].size());
        arrs[els] += rbuf[els];
      }
    }
  }
  mpi::wait_all(reqs.begin(), reqs.end());


  for (int d = 0; d < fastlog2(lastpower); d++) {
    for (int k = 0; k < lastpower; k += 1 << (d + 1)) {
      const int receiver = k;
      const int sender = k + (1 << d);
      if (rank == receiver) {

        for(size_t els=0; els<arrs.size(); els++) {
          grid.comm.recv(sender, tag+els, rbuf[els].data(), rbuf[els].size());
          arrs[els] += rbuf[els];
        }
      }
      else if (rank == sender) {
        for(size_t els=0; els<arrs.size(); els++) {
          grid.comm.send(receiver, tag+els, arrs[els].data(), arrs[els].size());
        }
      }
    }
  }


}



//template<size_t D>
//inline bool h5io::TestPrtclWriter<D>::write(
//    corgi::Node<D>& grid, int lap)
//{
//  read_tiles(grid);
//  mpi_reduce_snapshots(grid);
//
//  if( grid.comm.rank() == 0 ) {
//    // build filename
//    std::string full_filename = 
//      fname + 
//      "/flds" +
//      /*std::to_string(grid.comm.rank()) + */
//      "_" +
//      std::to_string(lap) +
//      extension;
//    //std::cout << "QW: " << full_filename << std::endl;
//
//    // open file and write
//    File file(full_filename, H5F_ACC_TRUNC);
//    file["Nx"] = arrs[0].Nx;
//    file["Ny"] = arrs[0].Ny;
//    file["Nz"] = arrs[0].Nz;
//
//    file["ex"] = arrs[0].serialize();
//    file["ey"] = arrs[1].serialize();
//    file["ez"] = arrs[2].serialize();
//
//    file["bx"] = arrs[3].serialize();
//    file["by"] = arrs[4].serialize();
//    file["bz"] = arrs[5].serialize();
//
//    file["jx"] = arrs[6].serialize();
//    file["jy"] = arrs[7].serialize();
//    file["jz"] = arrs[8].serialize();
//
//    file["rho"]= arrs[9].serialize();
//  }
//
//  return true;
//}






//--------------------------------------------------
// explicit template class instantiations
template class h5io::TestPrtclWriter<2>;
