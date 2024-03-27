#pragma once

#include <vector>
#include <string>

#include "io/snapshots/snapshot.h"
#include "corgi/corgi.h"
#include "tools/mesh.h"


namespace h5io { 


/// IO object for handling and storing sub-selection of particles 
//
// Particles are selected evenly from the initial condition. This is done
// by collecting their ids and then tagging those that are multiples of
// some stride value. Full history of these pre-identified particles is
// then saved to disk.
template<size_t D>
class TestPrtclWriter :
  public SnapshotWriter<D>
{

  public:

    using SnapshotWriter<D>::fname;
    using SnapshotWriter<D>::extension;
    using SnapshotWriter<D>::arrs;
    using SnapshotWriter<D>::rbuf;

  public:

    std::vector< toolbox::Mesh<int,0> > arrs2;
    std::vector< toolbox::Mesh<int,0> > rbuf2;


    /// general file name used for outputs
    const string file_name = "test-prtcls";

    /// do not consider particles beyond this id
    long cutoff_id;

    /// how many time steps to save
    int nt=1;

    /// total number of (test) particles
    long np;

    /// how many ranks/categories of particles (second id parameter)
    long nr;

    /// particle species/container to read/write
    int ispc = 0;

    /// data stride length
    long stride = 1;

    /// constructor that creates a name and creates uniform sampling of particles
    TestPrtclWriter(
        const std::string& prefix, 
        int Nx, int NxMesh,
        int Ny, int NyMesh,
        int Nz, int NzMesh,
        int ppc, int n_local_tiles,
        int n_test_particles_approx);


    /// read tile meshes into memory
    void read_tiles(corgi::Grid<D>& grid) override;

    /// write hdf5 file
    bool write(corgi::Grid<D>& grid, int lap) override;

    /// communicate snapshots with a B-tree cascade to rank 0
    // NOTE: this is modified to send 2 sets of arrays because we need
    // float + int separately
    void mpi_reduce_snapshots(corgi::Grid<D>& grid) override;

};


} // end of namespace h5io

