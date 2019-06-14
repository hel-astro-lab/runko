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

    #def __init__(self):
    #    print("initialized...")


#load tiles into each node
def loadTiles1D(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #if n.get_mpi_grid(i) == n.rank:
            c = pyrunko.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.add_tile(c, (i,) ) 


#load tiles into each node
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

        print("creating node")
        node = pycorgi.oneD.Node(conf.Nx, 1, 1)
        print("creating node")
        node.set_grid_lims(conf.xmin, conf.xmax)

        loadTiles1D(node, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(node.get_Nx()):
            #if n.get_mpi_grid(i,j) == n.rank:
            if True:
                c = node.get_tile(i)
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
        

        for cid in node.get_tile_ids():
            c = node.get_tile( cid )
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
        for cid in node.get_tile_ids():
            c = node.get_tile( cid )
            c.update_boundaries(node)

        ref = np.zeros(( conf.Nx*conf.Ny*conf.Nz, 
            conf.Nx*conf.NxMesh, 
            conf.Ny*conf.NyMesh, 
            conf.Nz*conf.NzMesh))

        m = 0
        for cid in node.get_tile_ids():
            c = node.get_tile( cid )
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

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(node, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(node.get_Nx()):
            for j in range(node.get_Ny()):
                #if n.get_mpi_grid(i,j) == n.rank:
                if True:
                    c = node.get_tile(i,j)
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

        for cid in node.get_tile_ids():
            c = node.get_tile( cid )
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
        for cid in node.get_tile_ids():
            c = node.get_tile( cid )
            c.update_boundaries(node)
        #for i in [1]:
        #    for j in [1]:
        #        c = node.get_tile(i,j)
        #        c.update_boundaries_2d(node)

        ref = np.zeros((conf.Nx*conf.Ny*conf.Nz, conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh))

        m = 0
        for cid in node.get_tile_ids():
            c = node.get_tile( cid )
            (i, j) = c.index
            yee = c.get_yee(0)

            for s in range(-1, conf.NzMesh+1, 1):
                for r in range(-1, conf.NyMesh+1, 1):
                    for q in range(-1, conf.NxMesh+1, 1):
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
                




