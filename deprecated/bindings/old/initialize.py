from initialize_pic import spatialLoc #XXX make this general/not pic-specific

import numpy as np

import pycorgi
import pyrunko 


#make random starting order
def loadMpiRandomly(n):
    np.random.seed(0)
    if n.master:
        for i in range(n.get_Nx()):
            for j in range(n.get_Ny()):
                val = np.random.randint( n.size() )
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()


#load nodes to be in stripe formation (splitted in X=horizontal direction)
def loadMpiXStrides(n):
    if n.master: #only master initializes; then sends
        stride = np.zeros( (n.get_Nx()), np.int64)
        dx = np.float(n.get_Nx()) / np.float(n.size() ) 
        for i in range(n.get_Nx()):
            val = np.int( i/dx )
            stride[i] = val

        for j in range(n.get_Ny()):
            for i in range(n.get_Nx()):
                val = stride[i]
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()


#load nodes to be in stripe formation (splitted in Y=vertical direction)
def loadMpiYStrides(n):
    if n.master: #only master initializes; then sends
        stride = np.zeros( (n.get_Ny()), np.int64)
        dy = np.float(n.get_Ny()) / np.float(n.Nrank) 
        for j in range(n.get_Ny()):
            val = np.int( j/dy )
            stride[j] = val

        for i in range(n.get_Nx()):
            for j in range(n.get_Ny()):
                val = stride[j]
                n.set_mpi_grid(i, j, val)
    n.bcast_mpi_grid()


# load nodes using 2D Hilbert curve
def loadMpi2D(n, comm_size=None):
    if n.master: #only master initializes; then sends

        if comm_size == None:
            comm_size = n.size()

        nx = n.get_Nx()
        ny = n.get_Ny()

        m0 = np.log2(nx)
        m1 = np.log2(ny)

        if not(m0.is_integer()):
            raise ValueError('Nx is not power of 2 (i.e. 2^m)')

        if not(m1.is_integer()):
            raise ValueError('Ny is not power of 2 (i.e. 2^m)')

        #print('Generating hilbert with 2^{} {}'.format(m0,m1))
        hgen = pyrunko.tools.twoD.HilbertGen(np.int(m0), np.int(m1))

        igrid = np.zeros( (nx, ny), np.int64)
        grid  = np.zeros( (nx, ny) ) #, np.int64)

        for i in range(nx):
            for j in range(ny):
                grid[i,j] = hgen.hindex(i,j)
        #print(grid)
        hmin, hmax = np.min(grid), np.max(grid)
        for i in range(nx):
            for j in range(ny):
                igrid[i,j] = np.floor( comm_size*grid[i,j]/(hmax+1) ) 

        #check that nodes get about same work load
        #y = np.bincount(igrid.flatten())
        #ii = np.nonzero(y)[0]
        #print(list(zip(ii,y[ii])))

        #print("grid:")
        for i in range(nx):
            for j in range(ny):
                val = igrid[i,j] 
                n.set_mpi_grid(i, j, val)
                #print("({},{}) = {}".format(i,j,val))


    n.bcast_mpi_grid()


# load nodes using 3D Hilbert curve
def loadMpi3D(n, comm_size=None):
    if n.master: #only master initializes; then sends
        if comm_size == None:
            comm_size = n.size()

        nx = n.get_Nx()
        ny = n.get_Ny()
        nz = n.get_Nz()

        m0 = np.log2(nx)
        m1 = np.log2(ny)
        m2 = np.log2(nz)

        if not(m0.is_integer()):
            raise ValueError('Nx is not power of 2 (i.e. 2^m)')

        if not(m1.is_integer()):
            raise ValueError('Ny is not power of 2 (i.e. 2^m)')

        if not(m2.is_integer()):
            raise ValueError('Nz is not power of 2 (i.e. 2^m)')


        #print('Generating hilbert with 2^{} {}'.format(m0,m1))
        hgen = pyrunko.tools.threeD.HilbertGen( np.int(m0), np.int(m1), np.int(m2) )

        igrid = np.zeros( (nx, ny, nz), np.int64)
        grid  = np.zeros( (nx, ny, nz) ) #, np.int64)

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    grid[i,j,k] = hgen.hindex(i,j,)

        #print(grid)
        hmin, hmax = np.min(grid), np.max(grid)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    igrid[i,j,k] = np.floor( comm_size*grid[i,j,k]/(hmax+1) ) 

        #check that nodes get about same work load
        #y = np.bincount(igrid.flatten())
        #ii = np.nonzero(y)[0]
        #print(list(zip(ii,y[ii])))

        #print("grid:")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    val = igrid[i,j,k] 
                    n.set_mpi_grid(i, j, k, val)
                    #print("({},{}) = {}".format(i,j,val))

    #broadcast calculation to everybody
    n.bcast_mpi_grid()



def initialize_tile(c, i, j, n, conf):

    #initialize tile dimensions 
    c.cfl = conf.cfl
    c.dx = conf.dx

    # initialize analysis tiles ready for incoming simulation data
    for ip in range(conf.Nspecies):
        c.add_analysis_species()

    #set bounding box of the tile 
    mins = spatialLoc(n, [i,j], [0,0,0], conf)
    maxs = spatialLoc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
    c.set_tile_mins(mins[0:1])
    c.set_tile_maxs(maxs[0:1])
    c.threshold = conf.clipThreshold


#load tiles into each grid
def loadTiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

            if n.get_mpi_grid(i,j) == n.rank():
                c = pyrunko.vlv.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                initialize_tile(c, i, j, n, conf)

                #add it to the grid
                n.add_tile(c, (i,)) 




#def initializeEmptyVelocityMesh(mesh, conf):
#    mesh.Nblocks = [conf.Nvx, conf.Nvy, conf.Nvz]
#
#    pmins = [conf.vxmin, conf.vymin, conf.vzmin]
#    pmaxs = [conf.vxmax, conf.vymax, conf.vzmax]
#    mesh.z_fill( pmins, pmaxs )
#
#
#
## create empty vmesh object according to conf specifications
#def createEmptyVelocityMesh(conf):
#    mesh = ptools.VeloMesh()
#    initializeEmptyVelocityMesh(mesh, conf)
#
#    return mesh
#









