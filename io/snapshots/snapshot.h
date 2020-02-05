#pragma once

#include <vector>
#include <array>
#include <string>

#include <mpi4cpp/mpi.h>
#include "../../definitions.h"
#include "../../corgi/corgi.h"
#include "../../tools/mesh.h"
#include "../namer.h"
#include "../../tools/fastlog.h"


namespace h5io { 


template<size_t D>
class SnapshotWriter {

  public:

    /// general file extension to be appended to file names
    const string extension = ".h5";
    
    /// Object to handle file names and extensions
    std::string fname;

    /// meshes
    std::vector< toolbox::Mesh<real_short> > arrs;

    /// mpi receive buffer
    std::vector< toolbox::Mesh<real_short> > rbuf;


    /// data stride length
    int stride = 1;

    /// constructor that creates a name and opens the file handle
    SnapshotWriter( const std::string& prefix ) : fname{prefix} { }

    // NOTE: modify these 2 functions to make your own snapshot io

    /// read tile meshes into memory
    virtual void read_tiles(corgi::Grid<D>& grid) = 0;

    /// write hdf5 file
    virtual bool write(corgi::Grid<D>& grid, int lap) = 0;

    /// communicate snapshots with a B-tree cascade to rank 0
    virtual void mpi_reduce_snapshots(corgi::Grid<D>& grid)
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

      std::vector<mpi4cpp::mpi::request> reqs;
      for (int i = lastpower; i < size; i++) {
        if (rank == i) {
          for(size_t els=0; els<arrs.size(); els++) {
            reqs.push_back( 
                grid.comm.isend(
                  i-lastpower, 
                  tag+els, 
                  arrs[els].data(), 
                  arrs[els].size()) 
                );
          }
        }
      }

      for (int i = 0; i < size-lastpower; i++) {
        if (rank == i) {
          for(size_t els=0; els<arrs.size(); els++) {
            grid.comm.recv(
                i+lastpower, 
                tag+els, 
                rbuf[els].data(), 
                rbuf[els].size()
                );
            arrs[els] += rbuf[els];
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
              grid.comm.recv(
                  sender, 
                  tag+els, 
                  rbuf[els].data(), 
                  rbuf[els].size()
                  );
              arrs[els] += rbuf[els];
            }
          }
          else if (rank == sender) {
            for(size_t els=0; els<arrs.size(); els++) {
              grid.comm.send(
                  receiver, 
                  tag+els, 
                  arrs[els].data(), 
                  arrs[els].size()
                  );
            }
          }
        }
      }
  }
};



}
