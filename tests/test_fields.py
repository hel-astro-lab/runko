import unittest

import sys
sys.path.append('python')
import numpy as np

import pyplasmaDev
import pyplasma as plasma


class Conf:

    Nx = 3
    Ny = 3

    NxMesh = 3
    NyMesh = 3

    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

    #def __init__(self):
    #    print("initialized...")


#load cells into each node
def loadCells(n):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            if n.getMpiGrid(i,j) == n.rank:
                c = plasma.PlasmaCell(i, j, n.rank, n.getNx(), n.getNy(), 3, 3)
                #c = plasma.VlasovCell(i, j, n.rank, n.getNx(), n.getNy(), 3, 3)
                n.addCell(c) #TODO load data to cell


def wrap(ii, N):
    if ii < 0:
        ii = N-1
        return N-1
    if ii == N:
        return 0
    return ii


class Communications(unittest.TestCase):

    """ Load cells, put running values in to them,
        update boundaries, check that every cell
        has correct boundaries.
    """
    def test_updateBoundaries(self):

        conf = Conf()

        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadCells(node)

        # lets put values into Yee lattice
        val = 1.0
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                #if n.getMpiGrid(i,j) == n.rank:
                if True:
                    c = node.getCellPtr(i,j)
                    yee = c.getYee(0)

                    for q in range(conf.NxMesh):
                        for k in range(conf.NyMesh):
                            yee.ex[q,k,0] = val
                            val += 1

        data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh))

        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            (i, j) = c.index()

            yee = c.getYee(0)

            for k in range(conf.NyMesh):
                for q in range(conf.NxMesh):
                    data[ i*conf.NxMesh + q, j*conf.NyMesh + k ] = yee.ex[q,k,0]

        #update boundaries
        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            c.updateBoundaries(node)


        ref = np.zeros((conf.Nx*conf.Ny, conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh))

        m = 0
        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            (i, j) = c.index()
            yee = c.getYee(0)

            for k in range(-1, conf.NyMesh+1, 1):
                for q in range(-1, conf.NxMesh+1, 1):
                    qq = wrap( i*conf.NxMesh + q, conf.Nx*conf.NxMesh )
                    kk = wrap( j*conf.NyMesh + k, conf.Ny*conf.NyMesh )
                    ref[m, qq, kk] = yee.ex[q,k,0]
            m += 1

        for m in range(conf.Nx*conf.Ny):
            for i in range(conf.Nx*conf.NxMesh):
                for j in range(conf.Ny*conf.NyMesh):
                    if ref[m,i,j] == 0:
                        continue

                    self.assertEqual( ref[m,i,j], data[i,j] )
                



