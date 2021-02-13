from mpi4py import MPI

import unittest

import sys
import numpy as np

import pycorgi
import pyrunko

import gc
#gc.set_debug(gc.DEBUG_LEAK|gc.DEBUG_STATS)

class Conf:

    Nx = 3
    Ny = 3
    Nz = 1

    one = False
    twoD = False
    threeD = False

    NxMesh = 5
    NyMesh = 5
    NzMesh = 5

    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

    zmin = 0.0
    zmax = 1.0

    qe = 1.0

    #def __init__(self):
    #    print("initialized...")


#load tiles into each grid
def loadTiles1D(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #if n.get_mpi_grid(i) == n.rank:
            c = pyrunko.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.add_tile(c, (i,) ) 


#load tiles into each grid
def loadTiles2D(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #if n.get_mpi_grid(i,j) == n.rank:
            c = pyrunko.fields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.add_tile(c, (i,j) ) 


#load tiles into each grid
def loadTiles3D(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            for k in range(n.get_Nz()):
                #print("putting", i,j,k)
                #if n.get_mpi_grid(i,j) == n.rank:
                c = pyrunko.fields.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                c.lengths = [conf.Nx, conf.Ny, conf.Nz]

                n.add_tile(c, (i,j,k) ) 
        
                #c = pyrunko.fields.threeD.make_and_add_tile( n, 
                #        conf.NxMesh, conf.NyMesh, conf.NzMesh, (i,j,k))



def wrap(ii, N):
    return (N + ii % N) % N
    #if ii < 0:
    #    return N-1
    #if ii > N:
    #    return 0
    #return ii


class FLD_inits(unittest.TestCase):

    def test_propagators_1d(self):
        conf = Conf()
        conf.oneD = True
        conf.NyMesh = 1 #force 1D
        conf.NzMesh = 1 #

        fdtd2 = pyrunko.fields.oneD.FDTD2()
        tile = pyrunko.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        fdtd2.push_e(tile)
        fdtd2.push_half_b(tile)
        
    def test_propagators_2d(self):
        conf = Conf()
        conf.twoD = True
        conf.NzMesh = 1 #force 2D

        fdtd2 = pyrunko.fields.twoD.FDTD2()
        tile = pyrunko.fields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        fdtd2.push_e(tile)
        fdtd2.push_half_b(tile)





class Communications(unittest.TestCase):

    """ Load tiles, put running values in to them,
        update boundaries, check that every tile
        has correct boundaries.
    """
    def test_updateBoundaries(self):

        conf = Conf()
        conf.oneD = True
        conf.NyMesh = 1 #force 1D
        conf.NzMesh = 1 #

        grid = pycorgi.oneD.Grid(conf.Nx, 1, 1)
        grid.set_grid_lims(conf.xmin, conf.xmax)

        loadTiles1D(grid, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(grid.get_Nx()):
            #if n.get_mpi_grid(i,j) == n.rank:
            if True:
                c = grid.get_tile(i)
                yee = c.get_yee(0)

                for q in range(conf.NxMesh):
                    for r in range(conf.NyMesh):
                        for s in range(conf.NzMesh):
                            yee.ex[q,r,s] = val
                            yee.ey[q,r,s] = val
                            yee.ez[q,r,s] = val
                                        
                            yee.bx[q,r,s] = val
                            yee.by[q,r,s] = val
                            yee.bz[q,r,s] = val
                                        
                            yee.jx[q,r,s] = val
                            yee.jy[q,r,s] = val
                            yee.jz[q,r,s] = val
                            val += 1

        data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))
        

        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i,) = c.index

            yee = c.get_yee(0)

            for q in range(conf.NxMesh):
                for r in range(conf.NyMesh):
                    for s in range(conf.NzMesh):
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 0] = yee.ex[q,r,s]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 1] = yee.ey[q,r,s]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 2] = yee.ez[q,r,s]
                                                                                     
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 3] = yee.bx[q,r,s]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 4] = yee.by[q,r,s]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 5] = yee.bz[q,r,s]
                                                                                     
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 6] = yee.jx[q,r,s]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 7] = yee.jy[q,r,s]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + r, 0*conf.NzMesh + s, 8] = yee.jz[q,r,s]

        #print("r=0-------")
        #print(data[:,:,0,0])
        #print("r=1-------")
        #print(data[:,:,1,0])
        #print("r=2-------")
        #print(data[:,:,2,0])

        #update boundaries
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            c.update_boundaries(grid)

        ref = np.zeros(( conf.Nx*conf.Ny*conf.Nz, 
            conf.Nx*conf.NxMesh, 
            conf.Ny*conf.NyMesh, 
            conf.Nz*conf.NzMesh))

        m = 0
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i,) = c.index
            yee = c.get_yee(0)

            for s in range(-1, conf.NzMesh+1, 1):
                for r in range(-1, conf.NyMesh+1, 1):
                    for q in range(-1, conf.NxMesh+1, 1):
                        #print("q k r ({},{},{})".format(q,k,r))
                        qq = wrap( i*conf.NxMesh + q, conf.Nx*conf.NxMesh )
                        rr = wrap( 0*conf.NyMesh + r, conf.Ny*conf.NyMesh )
                        ss = wrap( 0*conf.NzMesh + s, conf.Nz*conf.NzMesh )
                        ref[m, qq, rr, ss] = yee.ex[q,r,s]
            m += 1

        #print("cid = 0")
        #print(ref[0,:,:,0])
        #print(ref[0,:,:,1])
        #print(ref[0,:,:,2])

        # check every tile separately
        for m in range(conf.Nx*conf.Ny):
            for i in range(conf.Nx*conf.NxMesh):
                for j in range(conf.Ny*conf.NyMesh):
                    for k in range(conf.Nz*conf.NzMesh):
                        if ref[m,i,j,k] == 0:
                            continue

                        #loop over ex, ey, ez,...
                        for tp in range(9):
                            self.assertEqual( ref[m,i,j,k], data[i,j,k,tp] )
                


    def test_updateBoundaries2D(self):

        conf = Conf()
        conf.twoD = True

        conf.NxMesh = 4
        conf.NyMesh = 4
        conf.NzMesh = 1 #force 2D

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(grid, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                #if n.get_mpi_grid(i,j) == n.rank:
                if True:
                    c = grid.get_tile(i,j)
                    yee = c.get_yee(0)

                    for s in range(conf.NzMesh):
                        for r in range(conf.NyMesh):
                            for q in range(conf.NxMesh):
                                yee.ex[q,r,s] = val
                                yee.ey[q,r,s] = val
                                yee.ez[q,r,s] = val

                                yee.bx[q,r,s] = val
                                yee.by[q,r,s] = val
                                yee.bz[q,r,s] = val

                                yee.jx[q,r,s] = val
                                yee.jy[q,r,s] = val
                                yee.jz[q,r,s] = val
                                val += 1

        data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))

        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j) = c.index

            yee = c.get_yee(0)
            for q in range(conf.NxMesh):
                for r in range(conf.NyMesh):
                    for s in range(conf.NzMesh):
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 0] = yee.ex[q,r,s]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 1] = yee.ey[q,r,s]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 2] = yee.ez[q,r,s]
                                                                                                        
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 3] = yee.bx[q,r,s]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 4] = yee.by[q,r,s]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 5] = yee.bz[q,r,s]
                                                                                                        
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 6] = yee.jx[q,r,s]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 7] = yee.jy[q,r,s]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + r, 0*conf.NzMesh + s, 8] = yee.jz[q,r,s]

        #print("r=0-------")
        #print(data[:,:,0,0])
        #print("r=1-------")
        #print(data[:,:,1,0])
        #print("r=2-------")
        #print(data[:,:,2,0])

        #update boundaries
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            c.update_boundaries(grid)
        #for i in [1]:
        #    for j in [1]:
        #        c = grid.get_tile(i,j)
        #        c.update_boundaries_2d(grid)

        ref = np.zeros((conf.Nx*conf.Ny*conf.Nz, conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh))


        m = 0
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j) = c.index
            yee = c.get_yee(0)

            for s in range(-3, conf.NzMesh+3, 1):
                for r in range(-3, conf.NyMesh+3, 1):
                    for q in range(-3, conf.NxMesh+3, 1):
                        #print("q k r ({},{},{})".format(q,k,r))

                        qq = wrap( i*conf.NxMesh + q, conf.Nx*conf.NxMesh )
                        rr = wrap( j*conf.NyMesh + r, conf.Ny*conf.NyMesh )
                        ss = wrap( 0*conf.NzMesh + s, conf.Nz*conf.NzMesh )
                        #print( ref[m, qq, kk, rr]  )
                        #print(  yee.ex[q,k,r] )
                        ref[m, qq, rr, ss] = yee.ex[q,r,s]
            m += 1

    
        #print("--------------------------------------------------")
        #print("cid = 0")
        #print(ref[4,:,:,0])
        #print(ref[4,:,:,1])
        #print(ref[4,:,:,2])

        # check every tile separately
        for m in range(conf.Nx*conf.Ny):
            for i in range(conf.Nx*conf.NxMesh):
                for j in range(conf.Ny*conf.NyMesh):
                    for k in range(conf.Nz*conf.NzMesh):
                        if ref[m,i,j,k] == 0:
                            continue

                        #loop over ex, ey, ez,...
                        for tp in range(9):
                            self.assertEqual( ref[m,i,j,k], data[i,j,k,tp] )
                

        #print("--------------------------------------------------")
        #check halo regions of the middle tile

        c = grid.get_tile(1,1)
        yee = c.get_yee(0)
        arr = np.zeros((conf.NxMesh+6, conf.NyMesh+6))
        for r in range(-3, conf.NyMesh+3, 1):
            for q in range(-3, conf.NxMesh+3, 1):
                arr[q+3, r+3] = yee.ex[q,r,0]

        ref2 = data[1:-1, 1:-1, 0,0] #strip outer boundaries away
        nx,ny=np.shape(ref2)
        for i in range(nx):
            for j in range(ny):
                self.assertEqual(ref2[i,j], arr[i,j])

    # testing a spesific seg fault with loading of yee lattices. This same test fails with 3D
    def test_2D_tile_memory_bug(self):
        conf = Conf()
        conf.twoD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.NxMesh = 40
        conf.NyMesh = 40

        #print("mem bug --------")
        #print("grid")
        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, 1)
        #print("lims")
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        #print("load")
        tile = pyrunko.fields.twoD.Tile(conf.NxMesh, conf.NyMesh, 1)
        #print("ref count tile 2d:", sys.getrefcount(tile))
        grid.add_tile(tile, (0,0) ) 
        #print("ref count tile 2d:", sys.getrefcount(tile))

        #print("getting 000")
        c = grid.get_tile(0,0)
        yee = c.get_yee()

        #print("mem bug +++++++")

    def test_3D_mesh_memory_bug(self):
        conf = Conf()
        conf.threeD = True

        conf.NxMesh = 10
        conf.NyMesh = 10
        conf.NzMesh = 10

        #print("\n mesh mem bug --------")
        #print("separate pure python mesh")
        abc = pyrunko.tools.Mesh_H3(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        #print("ref count abc:", sys.getrefcount(abc))
        nx = abc.size()
        #print("size:", nx)

        H = 3
        self.assertEqual(nx, (conf.NxMesh+H*2)*(conf.NyMesh+H*2)*(conf.NzMesh+H*2))


    # Captures memory bug with mesh initialization; when internal meshes in YeeLattice
    # are too big, compiler tires to over-optimize. Then some stuff never gets allocated.
    # This is fixed now by a copy-and-swap algorithm in toolbox::Mesh.
    def test_3D_tile_memory_bug(self):

        conf = Conf()
        conf.threeD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3
        conf.NxMesh = 10
        conf.NyMesh = 10
        conf.NzMesh = 10

        #print("\n mem bug --------")
        #print("grid")
        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        #print("lims")
       #grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

        #Second create mechanism with automatic tying of pointer to grid lifetime
        if False:
            #print("create2")
            tile = pyrunko.fields.threeD.make_and_add_tile(
                    grid, 
                    conf.NxMesh, conf.NyMesh, conf.NzMesh,
                    (1,0,0)
                    )
        #First create mechanism with python-organized object lifetime
        else:
            #First default create mechanism
            # NOTE: leads to segfaulting with 3D tiles. Most likely py GC cleans tiles
            # aggressively. FIXED.
            #print("create")
            tile = pyrunko.fields.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            #print("ref count tile:", sys.getrefcount(tile))
            grid.add_tile(tile, (1,0,0) ) 
            #print("ref count tile:", sys.getrefcount(tile))
            #print("end of create")

        #tile = grid.get_tile(0,0,1)
        tile.set_tile_mins([1.0,1.0,1.0])
        tile.set_tile_maxs([2.0,2.0,2.0])
        #tile.add_analysis_species()

        #print("creating yee2")
        #yee2 = pyrunko.fields.YeeLattice(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        #print(sys.getrefcount(yee2))
        #print("ending yee2")


        if True:
            #print("getting yee p0")
            yee0 = tile.get_yee()
            #yee0 = tile.yee
            #yee0 = tile.get_yeeptr()
            #print("ref count yee0:", sys.getrefcount(yee0))
            #print("getting ex")
            ex0 = yee0.ex
            #ex0 = yee0.get_ex()
            #ex0 = yee0.ex2
            #print("getting size")
            nt = ex0.size()
            #print("size was", nt)

            nx = ex0.Nx
            ny = ex0.Ny
            nz = ex0.Nz
            #print("nx ny nz ", nx,ny,nz)

            self.assertEqual(nx, conf.NxMesh)
            self.assertEqual(ny, conf.NyMesh)
            self.assertEqual(nz, conf.NzMesh)
            #print("getting corner")
            val = ex0[1,1,1]

            #print("setting same corner")
            ex0[1,1,1] = 2.0
            #print("end of getting yee p0")

        # now re-ask for the tile to test handling of multiple copies of same object
        if True:
            #print("---phase 2---")
            #print("fetching another copy of the same tile")
            #print("getting 000")
            #c = grid.get_tile(0,0,1)
            c = grid.get_tile(1,0,0)
            #print(c.cid)
            #print(c.mins) 
            #print(c.maxs) 
            #print(c.index) 

            #print("ref count tile:", sys.getrefcount(c))

            #print("getting yee p1")
            yee = c.get_yee()
            #print("ref count yee:", sys.getrefcount(yee))
            #yee = c.yee
            #yee = c.get_yeeptr()
            #print("getting ex")
            ex = yee.ex
            #print("ref count ex:", sys.getrefcount(ex))
            #ex = yee.get_ex()

            nt = ex.size()
            #print("size was", nt)
            nx = ex.Nx
            ny = ex.Ny
            nz = ex.Nz
            #print("nx ny nz ", nx,ny,nz)


            #print("getting corner")
            val = ex[1,1,1]
            #print("val = ", val)
            #print("getting size")
            nt = ex.size()
            #print("size from yee", nt)

            #print("size from tile directly", c.yee.ex.nx)
            #self.assertEqual(nx, conf.NxMesh)

        #analysis = c.get_analysis()

        #now ask even more references
        if False:
            c1= grid.get_tile(0,0,1)
            yee1 = c1.get_yee()
            ex1 = c1.ex

            c2= grid.get_tile(0,0,1)
            yee2 = c2.get_yee()
            ex2 = c2.ex

            c3= grid.get_tile(0,0,1)
            yee3 = c3.get_yee()
            ex3 = c3.ex

        #print("mem bug +++++++")

    def test_tile_indices2D(self):

        conf = Conf()
        conf.twoD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 1 

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(grid, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                #if n.get_mpi_grid(i,j) == n.rank:
                if True:
                    c = grid.get_tile(i,j)

                    # Capture 3D indexing bug here where k index
                    # did not match what was expected

                    indx1 = c.index
                    indx2 = c.communication.indices
                    indx3 = c.get_index(grid)

                    #print("what got out 1:",c.cid," <-> ", i,j," vs. indx ", indx1)
                    #print("what got out 2:",c.cid," <-> ", i,j," vs. indx ", indx2)
                    #print("what got out 3:",c.cid," <-> ", i,j," vs. indx ", indx3)

                    self.assertEqual(i, indx1[0])
                    self.assertEqual(j, indx1[1])

                    self.assertEqual(i, indx2[0])
                    self.assertEqual(j, indx2[1])

                    self.assertEqual(i, indx3[0])
                    self.assertEqual(j, indx3[1])


    def test_tile_indices3D(self):

        conf = Conf()
        conf.threeD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 3

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

        loadTiles3D(grid, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        c = grid.get_tile(i,j,k)

                        # Capture 3D indexing bug here where k index
                        # did not match what was expected

                        indx1 = c.index
                        indx2 = c.communication.indices
                        indx3 = c.get_index(grid)

                        #print("what got out 1:",c.cid," <-> ", i,j," vs. indx ", indx1)
                        #print("what got out 2:",c.cid," <-> ", i,j," vs. indx ", indx2)
                        #print("what got out 3:",c.cid," <-> ", i,j," vs. indx ", indx3)

                        self.assertEqual(i, indx1[0])
                        self.assertEqual(j, indx1[1])
                        self.assertEqual(k, indx1[2])

                        self.assertEqual(i, indx2[0])
                        self.assertEqual(j, indx2[1])
                        self.assertEqual(k, indx2[2])

                        self.assertEqual(i, indx3[0])
                        self.assertEqual(j, indx3[1])
                        self.assertEqual(k, indx3[2])


    def test_updateBoundaries3D(self):

        conf = Conf()
        conf.threeD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3
        conf.NxMesh = 10
        conf.NyMesh = 10
        conf.NzMesh = 10

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
        loadTiles3D(grid, conf)

        c = grid.get_tile(1,1,1)
        yee = c.get_yee() #FIXME: yee does not exist here

        orig = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh))

        #insert running values into Yee
        val = 1.0
        for k in range(grid.get_Nz()):
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        c = grid.get_tile(i,j,k)

                        indx = c.index
                        
                        lens = c.lengths
                        self.assertEqual(lens[0], conf.Nx)
                        self.assertEqual(lens[1], conf.Ny)
                        self.assertEqual(lens[2], conf.Nz)

                        self.assertEqual(i, indx[0])
                        self.assertEqual(j, indx[1])
                        self.assertEqual(k, indx[2])

                        yee = c.get_yee()

                        for s in range(conf.NzMesh):
                            for r in range(conf.NyMesh):
                                for q in range(conf.NxMesh):
                                    ix = i*conf.NxMesh + q 
                                    jy = j*conf.NyMesh + r 
                                    kz = k*conf.NzMesh + s
                                    orig[ix,jy,kz] = val #save original value for comparison

                                    #put same value to yee
                                    yee.ex[q,r,s] = val
                                    yee.ey[q,r,s] = val
                                    yee.ez[q,r,s] = val

                                    yee.bx[q,r,s] = val
                                    yee.by[q,r,s] = val
                                    yee.bz[q,r,s] = val

                                    yee.jx[q,r,s] = val
                                    yee.jy[q,r,s] = val
                                    yee.jz[q,r,s] = val
                                    val += 1


        #update boundaries
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            c.update_boundaries(grid)

        #collect back to this array; default its values to -2 to see errors better
        ref = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))
        ref[:,:,:,:,] = -2.0

        #print("build check boundaries array")
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j, k) = c.get_index(grid)
            yee = c.get_yee()

            for s in range(-3, conf.NzMesh+3, 1):
                if not(s < 0 or s >= conf.NzMesh):
                    continue
                for r in range(-3, conf.NyMesh+3, 1):
                    if not(r < 0 or r >= conf.NyMesh):
                        continue
                    for q in range(-3, conf.NxMesh+3, 1):
                        if not(q < 0 or q >= conf.NxMesh):
                            continue

                        ix = i*conf.NxMesh + q 
                        jy = j*conf.NyMesh + r 
                        kz = k*conf.NzMesh + s

                        qq = wrap( ix, conf.Nx*conf.NxMesh )
                        rr = wrap( jy, conf.Ny*conf.NyMesh )
                        ss = wrap( kz, conf.Nz*conf.NzMesh )

                        ref[qq, rr, ss, 0] = yee.ex[q,r,s]
                        ref[qq, rr, ss, 1] = yee.ey[q,r,s]
                        ref[qq, rr, ss, 2] = yee.ez[q,r,s]

                        ref[qq, rr, ss, 3] = yee.bx[q,r,s]
                        ref[qq, rr, ss, 4] = yee.by[q,r,s]
                        ref[qq, rr, ss, 5] = yee.bz[q,r,s]

                        ref[qq, rr, ss, 6] = yee.jx[q,r,s]
                        ref[qq, rr, ss, 7] = yee.jy[q,r,s]
                        ref[qq, rr, ss, 8] = yee.jz[q,r,s]


        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j, k) = c.index
            yee = c.get_yee(0)
            for s in range(conf.NzMesh):
                for r in range(conf.NyMesh):
                    for q in range(conf.NxMesh):
                        ix = i*conf.NxMesh + q 
                        jy = j*conf.NyMesh + r 
                        kz = k*conf.NzMesh + s

                        val = 0
                        val = orig[ix,jy,kz]

                        #inner regions should be the default -2
                        if (s >= 3 and s < conf.NzMesh-3) or (r >= 3 and r < conf.NyMesh-3) or (q >= 3 and q < conf.NxMesh-3):
                            val = -2.0

                        for iarr in range(9):
                            self.assertEqual( ref[ix,jy,kz,iarr], val )

        #print for debug
        if False:
            large_width = 400
            with np.printoptions(linewidth=large_width):
                print("yee")
                print(ref[:,:,2,0].astype(int))
                print("orig")
                print(orig[:,:,2].astype(int))
        


    def test_neighs3D(self):

        conf = Conf()
        conf.threeD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 3 

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
        loadTiles3D(grid, conf)


        tile = grid.get_tile(1,1,1)
        for ii in [-1, 0, 1]:
            for jj in [-1, 0, 1]:
                for kk in [-1, 0, 1]:
                    if (ii == 0 and jj == 0 and kk == 0): continue

                    yee = tile.get_yee()
                    ind = tile.neighs(ii,jj,kk)


    def test_exchangeCurrents3D(self):

        conf = Conf()
        conf.threeD = True

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3
        conf.NxMesh = 6
        conf.NyMesh = 6
        conf.NzMesh = 6

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
        loadTiles3D(grid, conf)

        orig = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh))

        #insert running values into Yee
        val = 1.0 #base value that is summed 

        for k in range(grid.get_Nz()):
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        c = grid.get_tile(i,j,k)
                        yee = c.get_yee()

                        #loop over only halos
                        for s in range(-3, conf.NzMesh+3, 1):
                            for r in range(-3, conf.NyMesh+3, 1):
                                for q in range(-3, conf.NxMesh+3, 1):

                                    overlapx = 0
                                    overlapy = 0
                                    overlapz = 0

                                    if s < 0 or s > conf.NzMesh-1:
                                        overlapz += 1

                                    if r < 0 or r > conf.NyMesh-1:
                                        overlapy += 1

                                    if q < 0 or q > conf.NxMesh-1:
                                        overlapx += 1

                                    #consider values only outside inner regime 
                                    olap = overlapx + overlapy + overlapz
                                    if olap == 0:
                                        continue

                                    #corresponding true location 
                                    ix = i*conf.NxMesh + q 
                                    jy = j*conf.NyMesh + r 
                                    kz = k*conf.NzMesh + s
                                    qq = wrap( ix, conf.Nx*conf.NxMesh )
                                    rr = wrap( jy, conf.Ny*conf.NyMesh )
                                    ss = wrap( kz, conf.Nz*conf.NzMesh )

                                    #print("ijk {},{},{} | qrs {},{},{} | ii,jj,kk {},{},{} | qq rr ss {},{},{} | multxyz {},{},{}".format(i,j,k, q,r,s, ix,jy,kz, qq,rr,ss, i*conf.NxMesh, j*conf.NyMesh, k*conf.NzMesh))

                                    #add to halos
                                    yee.jx[q,r,s] = 1.0
                                    yee.jy[q,r,s] = 1.0
                                    yee.jz[q,r,s] = 1.0

                                    orig[qq,rr,ss] += 1


        #subtract one because the tile itself never contributes to sums
        orig[:,:,:] -= 1


        #update boundaries
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j, k) = c.get_index(grid)
            c.exchange_currents(grid)

        #collect back to this array; default its values to -2 to see errors better
        ref = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 3))

        #print("build check boundaries array")
        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j, k) = c.get_index(grid)
            yee = c.get_yee()

            #loop over internals where summation is put
            for s in range(conf.NzMesh):
                for r in range(conf.NyMesh):
                    for q in range(conf.NxMesh):
                        ix = i*conf.NxMesh + q 
                        jy = j*conf.NyMesh + r 
                        kz = k*conf.NzMesh + s

                        ref[ix,jy,kz, 0] = yee.jx[q,r,s]
                        ref[ix,jy,kz, 1] = yee.jy[q,r,s]
                        ref[ix,jy,kz, 2] = yee.jz[q,r,s]


        #print for debug
        if False:
            large_width = 400
            with np.printoptions(linewidth=large_width):
                print("yee")
                print(ref[:,:,3,0].astype(int))
                print("orig")
                print(orig[:,:,3].astype(int))

        for cid in grid.get_tile_ids():
            c = grid.get_tile( cid )
            (i, j, k) = c.index
            for s in range(conf.NzMesh):
                for r in range(conf.NyMesh):
                    for q in range(conf.NxMesh):
                        ix = i*conf.NxMesh + q 
                        jy = j*conf.NyMesh + r 
                        kz = k*conf.NzMesh + s

                        val = orig[ix,jy,kz]

                        #print("{},{},{}={}".format(ix,jy,kz,val))

                        for iarr in range(3):
                            self.assertEqual( ref[ix,jy,kz,iarr], val )


        
