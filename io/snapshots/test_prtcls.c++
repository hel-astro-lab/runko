#include <mpi4cpp/mpi.h>

#include "io/snapshots/test_prtcls.h"
#include "tools/ezh5/src/ezh5.hpp"
#include "pic/particle.h"
#include "pic/tile.h"

using ezh5::File;

template<size_t D>
h5io::TestPrtclWriter<D>::TestPrtclWriter(
        const std::string& prefix, 
        int Nx_in, int NxMesh_in,
        int Ny_in, int NyMesh_in,
        int Nz_in, int NzMesh_in,
        int ppc_in, int n_local_tiles_in,
        int n_test_particles_approx_in) :
  SnapshotWriter<D>{prefix}
{

  // enforce long accuracy due to persistent overflows
  auto Nx = static_cast<long>(Nx_in);
  auto Ny = static_cast<long>(Ny_in);
  auto Nz = static_cast<long>(Nz_in);

  auto NxMesh = static_cast<long>(NxMesh_in);
  auto NyMesh = static_cast<long>(NyMesh_in);
  auto NzMesh = static_cast<long>(NzMesh_in);

  auto ppc = static_cast<long>(ppc_in);
  auto n_local_tiles = static_cast<long>(n_local_tiles_in);
  auto n_test_particles_approx = static_cast<long>(n_test_particles_approx_in);


  //fname = prefix + "-" + to_string(lap) + extension;
    
  // internal mpi calls
  int n_comm_size_in, my_rank_in;
  MPI_Comm_size(MPI_COMM_WORLD, &n_comm_size_in);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_in);

  auto n_comm_size = static_cast<long>(n_comm_size_in);
  auto my_rank     = static_cast<long>(my_rank_in);

  // total number of particles
  long n_particles = Nx*NxMesh*Ny*NyMesh*Nz*NzMesh*ppc;

  // corresponding thinning factor is then (approximately)
  stride = n_particles/n_test_particles_approx;

  // avoid under/overflows
  if(stride > NxMesh*NyMesh*NzMesh*ppc) {
    if(my_rank == 0) {
      std::cout << "WARNING TestPrtclWriter underflow; fallbacking to minimum number of particles\n";
    }
    stride = NxMesh*NyMesh*NzMesh*ppc;
  }
  //stride = stride > NxMesh*NyMesh*NzMesh*ppc ? NxMesh*NyMesh*NzMesh*ppc : stride;

  if( stride <= 0 ) {
    if(my_rank == 0) {
      std::cout << "WARNING TestPrtclWriter overflow; fallbacking to minimum number of particles\n";
    }
    stride = NxMesh*NyMesh*NzMesh*ppc;
  }

  // how many particles would this particular rank have
  long my_n_prtcls = (n_local_tiles*NxMesh*NyMesh*NzMesh*ppc)/stride;

  // global minimum number of particles per rank
  long min_n_prtcls;

  MPI_Allreduce(&my_n_prtcls, &min_n_prtcls, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);

  // cutoff id is then (to satisfy min_n_prtcls in our search)
  cutoff_id = stride*min_n_prtcls;

  // calculate true number of test particles
  long n_test_particles = min_n_prtcls*n_comm_size;


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


  //nt=1;
  np = min_n_prtcls;
  nr = n_comm_size;

  // 7+6 (loc, vel, wgt and ex,ey,ez,bz,by,bz), variables per particle to save
  for(size_t i=0; i<=12; i++) arrs.emplace_back(nt, np, nr);
  for(size_t i=0; i<=12; i++) rbuf.emplace_back(nt, np, nr);

  // int arrays
  for(size_t i=0; i<=1; i++) arrs2.emplace_back(nt, np, nr);
  for(size_t i=0; i<=1; i++) rbuf2.emplace_back(nt, np, nr);
}




template<size_t D>
inline void h5io::TestPrtclWriter<D>::read_tiles(
    corgi::Grid<D>& grid)
{
  // clear target arrays
  for(auto& arr : arrs ) arr.clear();
  for(auto& arr : arrs2) arr.clear();

  // target arrays
  auto& xloc = arrs[0];
  auto& yloc = arrs[1];
  auto& zloc = arrs[2];
  auto& ux   = arrs[3];
  auto& uy   = arrs[4];
  auto& uz   = arrs[5];
  auto& wgt  = arrs[6];

  auto& exp  = arrs[7];
  auto& eyp  = arrs[8];
  auto& ezp  = arrs[9];
  auto& bxp  = arrs[10];
  auto& byp  = arrs[11];
  auto& bzp  = arrs[12];

  auto& ids   = arrs2[0];
  auto& procs = arrs2[1];

  int ip, ir;
  int tstep = 0;

  // read my local tiles
  for(auto cid : grid.get_local_tiles() ){
    auto& tile = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    auto& container = tile.get_container( ispc );
    int nparts = container.size();

    float_p* loc[3];
    for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

    float_p* vel[3];
    for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

    float_p* ch;
    ch = &( container.wgt(0) );

    float_p *ex, *ey, *ez, *bx, *by, *bz;
    ex = &( container.Epart[0*nparts] );
    ey = &( container.Epart[1*nparts] );
    ez = &( container.Epart[2*nparts] );

    bx = &( container.Bpart[0*nparts] );
    by = &( container.Bpart[1*nparts] );
    bz = &( container.Bpart[2*nparts] );

    // reference the ids 
    int* idn[2];
    for(int i=0; i<2; i++) idn[i] = &( container.id(i,0) );

    // long length id
    long idnl0, idnl1;

    // loop and search over all particles
    for(size_t n=0; n<container.size(); n++) {
      idnl0 = static_cast<long>(idn[0][n]);
      idnl1 = static_cast<long>(idn[1][n]);

      // skip beyond cutoff
      if(idnl0 >= cutoff_id) continue;

      // process only multiples of stride
      if(idnl0 % stride == 0) {
        ip = static_cast<int>(idnl0/stride); // get id
        ir = static_cast<int>(idnl1);        // get proc

        //std::cout << "...prtcl: " <<
        //  ip << " " <<
        //  ir << " " <<
        //  idn[0][n] << " " <<
        //  idn[1][n] << "\n";

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

        exp(  tstep, ip, ir) = ex[n];
        eyp(  tstep, ip, ir) = ey[n];
        ezp(  tstep, ip, ir) = ez[n];

        bxp(  tstep, ip, ir) = bx[n];
        byp(  tstep, ip, ir) = by[n];
        bzp(  tstep, ip, ir) = bz[n];

      }
    }


  }
}



template<size_t D>
inline void h5io::TestPrtclWriter<D>::mpi_reduce_snapshots(
    corgi::Grid<D>& grid)
{
  /* based on https://gist.github.com/rmcgibbo/7178576
   */

  for(auto& arr : rbuf)  arr.clear();
  for(auto& arr : rbuf2) arr.clear();

  int tag = 0;
  const int size = grid.comm.size(); //MPI::COMM_WORLD.Get_size();
  const int rank = grid.comm.rank(); //MPI::COMM_WORLD.Get_rank();
  const int lastpower = 1 << fastlog2(size);

  // each of the ranks greater than the last power of 2 less than size
  // need to downshift their data, since the binary tree reduction below
  // only works when N is a power of two.
  
  std::vector<mpi4cpp::mpi::request> reqs;
  for (int i = lastpower; i < size; i++) {
    if (rank == i) {
      for(size_t els=0; els<arrs.size(); els++) {
        reqs.push_back( grid.comm.isend(i-lastpower, tag+els, arrs[els].data(), arrs[els].size()) );
      }

      for(size_t els=0; els<arrs2.size(); els++) {
        reqs.push_back( grid.comm.isend(i-lastpower, tag+els+20, arrs2[els].data(), arrs2[els].size()) );
      }

    }
  }

  for (int i = 0; i < size-lastpower; i++) {
    if (rank == i) {
      for(size_t els=0; els<arrs.size(); els++) {
        grid.comm.recv(i+lastpower, tag+els, rbuf[els].data(), rbuf[els].size());
        arrs[els] += rbuf[els];
      }

      for(size_t els=0; els<arrs2.size(); els++) {
        grid.comm.recv(i+lastpower, tag+els+20, rbuf2[els].data(), rbuf2[els].size());
        arrs2[els] += rbuf2[els];
      }

    }
  }
  mpi4cpp::mpi::wait_all(reqs.begin(), reqs.end());


  for (int d = 0; d < fastlog2(lastpower); d++) {
    for (int k = 0; k < lastpower; k += 1 << (d + 1)) {
      const int receiver = k;
      const int sender = k + (1 << d);
      if (rank == receiver) {

        for(size_t els=0; els<arrs.size(); els++) {
          grid.comm.recv(sender, tag+els, rbuf[els].data(), rbuf[els].size());
          arrs[els] += rbuf[els];
        }

        for(size_t els=0; els<arrs2.size(); els++) {
          grid.comm.recv(sender, tag+els+20, rbuf2[els].data(), rbuf2[els].size());
          arrs2[els] += rbuf2[els];
        }

      }
      else if (rank == sender) {
        for(size_t els=0; els<arrs.size(); els++) {
          grid.comm.send(receiver, tag+els, arrs[els].data(), arrs[els].size());
        }

        for(size_t els=0; els<arrs2.size(); els++) {
          grid.comm.send(receiver, tag+els+20, arrs2[els].data(), arrs2[els].size());
        }


      }
    }
  }
}



template<size_t D>
inline bool h5io::TestPrtclWriter<D>::write(
    corgi::Grid<D>& grid, int lap)
{
  read_tiles(grid);
  mpi_reduce_snapshots(grid);

  if( grid.comm.rank() == 0 ) {

    // build filename
    std::string full_filename;

    if(ispc == 0) {
        full_filename = 
          fname + "/" +
          file_name + 
          "_" +
          std::to_string(lap) +
          extension;
    } else {
        full_filename = 
          fname + "/" +
          file_name + 
          "-" +
          std::to_string(ispc) +
          "_" +
          std::to_string(lap) +
          extension;
    }

    //std::cout << "QW: " << full_filename << std::endl;

    // open file and write
    File file(full_filename, H5F_ACC_TRUNC);
    file["Nx"] = arrs[0].Nx;
    file["Ny"] = arrs[0].Ny;
    file["Nz"] = arrs[0].Nz;


    file["x"]   = arrs[0].serialize();
    file["y"]   = arrs[1].serialize();
    file["z"]   = arrs[2].serialize();
    file["vx"]  = arrs[3].serialize();
    file["vy"]  = arrs[4].serialize();
    file["vz"]  = arrs[5].serialize();
    file["wgt"] = arrs[6].serialize();

    file["ex"]  = arrs[7].serialize();
    file["ey"]  = arrs[8].serialize();
    file["ez"]  = arrs[9].serialize();

    file["bx"]  = arrs[10].serialize();
    file["by"]  = arrs[11].serialize();
    file["bz"]  = arrs[12].serialize();

    file["id"]   = arrs2[0].serialize();
    file["proc"] = arrs2[1].serialize();

  }

  return true;
}




//--------------------------------------------------
// explicit template class instantiations
template class h5io::TestPrtclWriter<1>;
template class h5io::TestPrtclWriter<2>;
template class h5io::TestPrtclWriter<3>;
