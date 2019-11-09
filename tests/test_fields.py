from mpi4py import MPI

import unittest

import sys
import numpy as np

import pycorgi
import pyrunko


class Conf:

    Nx = 3
    Ny = 3
    Nz = 1

    NxMesh = 5
    NyMesh = 5
    NzMesh = 5

    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

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
        conf.NyMesh = 1 #force 1D
        conf.NzMesh = 1 #

        fdtd2 = pyrunko.fields.oneD.FDTD2()
        tile = pyrunko.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        fdtd2.push_e(tile)
        fdtd2.push_half_b(tile)
        
    def test_propagators_2d(self):
        conf = Conf()
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
        conf.NyMesh = 1 #force 1D
        conf.NzMesh = 1 #

        print("creating grid")
        grid = pycorgi.oneD.Grid(conf.Nx, 1, 1)
        print("creating grid")
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

