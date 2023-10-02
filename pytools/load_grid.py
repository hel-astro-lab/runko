# -*- coding: utf-8 -*- 

import numpy as np
import h5py as h5

import pycorgi
import pyrunko



def balance_mpi(n, conf, comm_size=None, mpi_master_mode=False):
    if conf.oneD:
        return balance_mpi_1D(n, comm_size=comm_size)
    if conf.twoD:
        return balance_mpi_2D(n, comm_size=comm_size)
    elif conf.threeD:
        return balance_mpi_3D(n, comm_size=comm_size, mpi_master_mode=mpi_master_mode)



#load nodes to be in stripe formation (splitted in X=horizontal direction)
def balance_mpi_1D(n, comm_size):
    if n.rank() == 0: #only master initializes; then sends
        if comm_size == None:
            comm_size = n.size()
        for i in range(n.get_Nx()):
            val = int(np.floor( i / n.get_Nx() * comm_size ))
            n.set_mpi_grid(i, val)
    n.bcast_mpi_grid()
    return

# load nodes using 2D Hilbert curve
def balance_mpi_2D(n, comm_size=None):

    if n.rank() == 0:  # only master initializes; then sends
        if comm_size == None:
            comm_size = n.size()

        nx = n.get_Nx()
        ny = n.get_Ny()

        m0 = np.log2(nx)
        m1 = np.log2(ny)

        if not (m0.is_integer()):
            raise ValueError("Nx is not power of 2 (i.e. 2^m)")

        if not (m1.is_integer()):
            raise ValueError("Ny is not power of 2 (i.e. 2^m)")

        # print('Generating hilbert with 2^{} {}'.format(m0,m1))
        hgen = pyrunko.tools.twoD.HilbertGen(int(m0), int(m1))

        igrid = np.zeros((nx, ny), int)
        grid = np.zeros((nx, ny))  # , int)

        for i in range(nx):
            for j in range(ny):
                grid[i, j] = hgen.hindex(i, j)

        # print(grid)
        hmin, hmax = np.min(grid), np.max(grid)
        for i in range(nx):
            for j in range(ny):
                igrid[i, j] = np.floor(comm_size * grid[i, j] / (hmax + 1))

        # check that nodes get about same work load
        # y = np.bincount(igrid.flatten())
        # ii = np.nonzero(y)[0]
        # print(list(zip(ii,y[ii])))

        # print("grid:")
        for i in range(nx):
            for j in range(ny):
                val = igrid[i, j]
                n.set_mpi_grid(i, j, val)
                    # print("({},{}) = {}".format(i,j,val))

    # broadcast calculation to everybody
    n.bcast_mpi_grid()

    return


# load nodes using 3D Hilbert curve
def balance_mpi_3D(n, comm_size=None, mpi_master_mode=False):

    if n.rank() == 0:  # only master initializes; then sends

        if comm_size == None:
            comm_size = n.size()

        # master mode does not allocate any tiles to rank =0
        if mpi_master_mode:
            comm_size -= 1

        nx = n.get_Nx()
        ny = n.get_Ny()
        nz = n.get_Nz()

        m0 = np.log2(nx)
        m1 = np.log2(ny)
        m2 = np.log2(nz)

        if not (m0.is_integer()):
            raise ValueError("Nx is not power of 2 (i.e. 2^m)")

        if not (m1.is_integer()):
            raise ValueError("Ny is not power of 2 (i.e. 2^m)")

        if not (m2.is_integer()):
            raise ValueError("Nz is not power of 2 (i.e. 2^m)")

        # print('Generating hilbert with 2^{} {}'.format(m0,m1))
        hgen = pyrunko.tools.threeD.HilbertGen(int(m0), int(m1), int(m2))

        igrid = np.zeros((nx, ny, nz), int)
        grid = np.zeros((nx, ny, nz))  # , int)

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    grid[i, j, k] = hgen.hindex(i, j, k)

        # print(grid)
        hmin, hmax = np.min(grid), np.max(grid)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    igrid[i, j, k] = np.floor(comm_size * grid[i, j, k] / (hmax + 1))

        if mpi_master_mode:
            igrid[:,:,:] += 1 # offset by one so that rank=0 is skipped

        # check that nodes get about same work load
        y = np.bincount(igrid.flatten())
        ii = np.nonzero(y)[0]
        #print("work load:")
        #print(list(zip(ii,y[ii])))

        # print("grid:")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    val = igrid[i, j, k]
                    n.set_mpi_grid(i, j, k, val)
                    # print("({},{}) = {}".format(i,j,val))

    # broadcast calculation to everybody
    n.bcast_mpi_grid()

    return


# load nodes in a recursive "catepillar track" type grid.
def load_catepillar_track_mpi(
        n, 
        nx_track_len,
        conf,
        comm_size=None):


    #if True:  # only master initializes; then sends
    if n.rank() == 0:  # only master initializes; then sends
        if comm_size == None:
            comm_size = n.size()

        #nx = 2**3
        nx = nx_track_len
        ny = n.get_Ny()
        nz = n.get_Nz()

        m0 = np.log2(nx)
        m1 = np.log2(ny)
        m2 = np.log2(nz)

        if not (m0.is_integer()):
            raise ValueError("Nx is not power of 2 (i.e. 2^m)")

        if not (m1.is_integer()):
            raise ValueError("Ny is not power of 2 (i.e. 2^m)")
        
        if conf.threeD and not( m2.is_integer() ):
            raise ValueError("Nz is not power of 2 (i.e. 2^m)")

        # print('Generating hilbert with 2^{} {}'.format(m0,m1))
        if conf.twoD:
            hgen = pyrunko.tools.twoD.HilbertGen(int(m0), int(m1))
        elif conf.threeD:
            hgen = pyrunko.tools.threeD.HilbertGen(int(m0), int(m1), int(m2) )
        grid = np.zeros((nx, ny, nz))  # , int)

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if conf.twoD:
                        grid[i, j, k] = hgen.hindex(i, j)
                    elif conf.threeD:
                        grid[i, j, k] = hgen.hindex(i, j, k)

        # print(grid)
        hmin, hmax = np.min(grid), np.max(grid)

        # create alternating true grid
        nxt = n.get_Nx()
        igrid = np.zeros((nxt, ny, nz), int)

        for i in range(nxt):
            for j in range(ny):
                for k in range(nz):
                    ic = i % nx
                    igrid[i, j, k] = np.floor(comm_size * grid[ic, j, k] / (hmax + 1))

        # check that nodes get about same work load
        y = np.bincount(igrid.flatten())
        ii = np.nonzero(y)[0]

        #print('rank load statistics')
        #print(list(zip(ii,y[ii])))

        for i in range(nxt):
            for j in range(ny):
                for k in range(nz):
                    val = igrid[i, j, k]

                    if conf.twoD:
                        n.set_mpi_grid(i, j, val)
                    elif conf.threeD:
                        n.set_mpi_grid(i, j, k, val)

    # broadcast calculation to everybody
    n.bcast_mpi_grid()

    return



#make random starting order
def load_mpi_randomly(n):
    np.random.seed(0)
    if n.rank() == 0:
        for i in range(n.get_Nx()):
            for j in range(n.get_Ny()):
                val = np.random.randint( n.size() )
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()
    return


#load nodes to be in stripe formation (splitted in X=horizontal direction)
def load_mpi_x_strides(n, conf):
    if n.rank() == 0: #only master initializes; then sends
        stride = np.zeros( (n.get_Nx()), int)
        dx = np.float64(n.get_Nx()) / np.float64(n.size() ) 
        for i in range(n.get_Nx()):
            val = int( i/dx )
            stride[i] = val

        for j in range(n.get_Ny()):
            for i in range(n.get_Nx()):
                val = stride[i]
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()
    return


#load nodes to be in stripe formation (splitted in Y=vertical direction)
def load_mpi_y_strides(n, conf):
    if n.rank() == 0: #only master initializes; then sends
        stride = np.zeros( (n.get_Ny()), int)
        dy = np.float64(n.get_Ny()) / np.float64(n.size() ) 
        for j in range(n.get_Ny()):
            val = int( j/dy )
            stride[j] = val

        for i in range(n.get_Nx()):
            for j in range(n.get_Ny()):
                val = stride[j]
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()
    return


# get mpi grid from the grid object
# NOTE: every rank returns their version. In general this should be the
# same for all because grids are kept up-to-date with mpi broadcasting from
# the master rank.
def get_mpi_grid(grid, conf):
    nx = grid.get_Nx()
    ny = grid.get_Ny()
    nz = grid.get_Nz()

    #mpi_grid = np.zeros((nx, ny, nz), int64)
    #mpi_grid = np.zeros((ny, nx, nz), int64)
    mpi_grid = np.zeros((nx, ny, nz), int)

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if conf.twoD:
                    val = grid.get_mpi_grid(i, j)
                elif conf.threeD:
                    val = grid.get_mpi_grid(i, j, k)

                mpi_grid[i,j,k] = val 

    #print('slice in get mpi:', mpi_grid[:,0,0])
    return mpi_grid


# multi-rank call to get tile ids from each rank; 
# NOTE: when ran with multiple ranks this only returns the local
#       list of tile ids.
def get_tile_id_grid(grid, conf):

    # get tile ids from everybody
    nx = grid.get_Nx()
    ny = grid.get_Ny()
    nz = grid.get_Nz()
    tile_grid = np.zeros( (nx,ny,nz), int )

    for cid in grid.get_tile_ids():
        tile = grid.get_tile( cid )

        if conf.twoD:
            (i, j) = tile.index
            k = 0
        elif conf.threeD:
            (i, j, k) = tile.index

        # test that there are no error tiles
        assert cid+1 > 0

        # finally, set id to correct location
        tile_grid[i,j,k] = cid + 1

    #TODO gather back to master info from all ranks

    return tile_grid



# save a snapshot of the mpi tile grid configuration to disk as h5 file
def save_mpi_grid_to_disk(outdir, lap, grid, conf):

    # get tile id grid; need to call with every rank
    #tid_grid = get_tile_id_grid(grid, conf)

    #only master saves rank performs this
    if grid.rank() == 0: 
        mpi_grid = get_mpi_grid(grid, conf)
        fname = outdir + '/mpi_{}.h5'.format(str(lap)) #.rjust(4,'0'))
        f5 = h5.File(fname,'w')

        dset1 = f5.create_dataset('mpi_grid',  data=mpi_grid)
        #dset2 = f5.create_dataset('tile_grid', data=tid_grid)

        f5.close()

