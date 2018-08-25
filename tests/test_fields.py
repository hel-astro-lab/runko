import unittest

import sys
import numpy as np

import pycorgi
import pyplasmabox


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
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i) == n.rank:
            c = pyplasmabox.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.addTile(c, (i,) ) 


#load tiles into each node
def loadTiles2D(n, conf):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            c = pyplasmabox.fields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.addTile(c, (i,j) ) 


def wrap(ii, N):
    return (N + ii % N) % N
    #if ii < 0:
    #    return N-1
    #if ii > N:
    #    return 0
    #return ii



class Communications(unittest.TestCase):

    """ Load tiles, put running values in to them,
        update boundaries, check that every tile
        has correct boundaries.
    """
    def test_updateBoundaries(self):

        conf = Conf()
        conf.NyMesh = 1 #force 1D
        conf.NzMesh = 1 #

        node = pycorgi.oneD.Node(conf.Nx, 1, 1)
        node.setGridLims(conf.xmin, conf.xmax)

        loadTiles1D(node, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(node.getNx()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                c = node.getTile(i)
                yee = c.getYee(0)

                for q in range(conf.NxMesh):
                    for k in range(conf.NyMesh):
                        for r in range(conf.NzMesh):
                            yee.ex[q,k,r] = val
                            yee.ey[q,k,r] = val
                            yee.ez[q,k,r] = val

                            yee.bx[q,k,r] = val
                            yee.by[q,k,r] = val
                            yee.bz[q,k,r] = val

                            yee.jx[q,k,r] = val
                            yee.jy[q,k,r] = val
                            yee.jz[q,k,r] = val
                            val += 1

        data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))
        

        for cid in node.getTileIds():
            c = node.getTile( cid )
            (i,) = c.index

            yee = c.getYee(0)

            for q in range(conf.NxMesh):
                for k in range(conf.NyMesh):
                    for r in range(conf.NzMesh):
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 0] = yee.ex[q,k,r]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 1] = yee.ey[q,k,r]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 2] = yee.ez[q,k,r]
                                                  
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 3] = yee.bx[q,k,r]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 4] = yee.by[q,k,r]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 5] = yee.bz[q,k,r]
                                                  
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 6] = yee.jx[q,k,r]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 7] = yee.jy[q,k,r]
                        data[ i*conf.NxMesh + q, 0*conf.NyMesh + k, 0*conf.NzMesh + r, 8] = yee.jz[q,k,r]

        #print("r=0-------")
        #print(data[:,:,0,0])
        #print("r=1-------")
        #print(data[:,:,1,0])
        #print("r=2-------")
        #print(data[:,:,2,0])

        #update boundaries
        for cid in node.getTileIds():
            c = node.getTile( cid )
            c.updateBoundaries(node)

        ref = np.zeros(( conf.Nx*conf.Ny*conf.Nz, 
            conf.Nx*conf.NxMesh, 
            conf.Ny*conf.NyMesh, 
            conf.Nz*conf.NzMesh))

        m = 0
        for cid in node.getTileIds():
            c = node.getTile( cid )
            (i,) = c.index
            yee = c.getYee(0)

            for r in range(-1, conf.NzMesh+1, 1):
                for k in range(-1, conf.NyMesh+1, 1):
                    for q in range(-1, conf.NxMesh+1, 1):
                        #print("q k r ({},{},{})".format(q,k,r))
                        qq = wrap( i*conf.NxMesh + q, conf.Nx*conf.NxMesh )
                        kk = wrap( 0*conf.NyMesh + k, conf.Ny*conf.NyMesh )
                        rr = wrap( 0*conf.NzMesh + r, conf.Nz*conf.NzMesh )
                        ref[m, qq, kk, rr] = yee.ex[q,k,r]
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
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(node, conf)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                #if n.getMpiGrid(i,j) == n.rank:
                if True:
                    c = node.getTile(i,j)
                    yee = c.getYee(0)

                    for q in range(conf.NxMesh):
                        for k in range(conf.NyMesh):
                            for r in range(conf.NzMesh):
                                yee.ex[q,k,r] = val
                                yee.ey[q,k,r] = val
                                yee.ez[q,k,r] = val

                                yee.bx[q,k,r] = val
                                yee.by[q,k,r] = val
                                yee.bz[q,k,r] = val

                                yee.jx[q,k,r] = val
                                yee.jy[q,k,r] = val
                                yee.jz[q,k,r] = val
                                val += 1

        data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))

        for cid in node.getTileIds():
            c = node.getTile( cid )
            (i, j) = c.index

            yee = c.getYee(0)
            for k in range(conf.NyMesh):
                for q in range(conf.NxMesh):
                    for r in range(conf.NzMesh):
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 0] = yee.ex[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 1] = yee.ey[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 2] = yee.ez[q,k,r]

                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 3] = yee.bx[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 4] = yee.by[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 5] = yee.bz[q,k,r]

                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 6] = yee.jx[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 7] = yee.jy[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 8] = yee.jz[q,k,r]

        #print("r=0-------")
        #print(data[:,:,0,0])
        #print("r=1-------")
        #print(data[:,:,1,0])
        #print("r=2-------")
        #print(data[:,:,2,0])

        #update boundaries
        for cid in node.getTileIds():
            c = node.getTile( cid )
            c.updateBoundaries(node)
        #for i in [1]:
        #    for j in [1]:
        #        c = node.getTile(i,j)
        #        c.updateBoundaries2D(node)

        ref = np.zeros((conf.Nx*conf.Ny*conf.Nz, conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh))

        m = 0
        for cid in node.getTileIds():
            c = node.getTile( cid )
            (i, j) = c.index
            yee = c.getYee(0)

            for r in range(-1, conf.NzMesh+1, 1):
                for k in range(-1, conf.NyMesh+1, 1):
                    for q in range(-1, conf.NxMesh+1, 1):
                        #print("q k r ({},{},{})".format(q,k,r))

                        qq = wrap( i*conf.NxMesh + q, conf.Nx*conf.NxMesh )
                        kk = wrap( j*conf.NyMesh + k, conf.Ny*conf.NyMesh )
                        rr = wrap( 0*conf.NzMesh + r, conf.Nz*conf.NzMesh )
                        #print( ref[m, qq, kk, rr]  )
                        #print(  yee.ex[q,k,r] )
                        ref[m, qq, kk, rr] = yee.ex[q,k,r]
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
                




