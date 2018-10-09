import unittest

import os

import numpy as np
import pycorgi
import pyplasmabox

from visualize import getYee
from visualize import getYee2D
from combine_files import combine_tiles
 
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


# create similar reference array
def fill_ref(node, conf):
    data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))

    # lets put floats into Yee lattice
    val = 1.0
    Nx = node.getNx()
    Ny = node.getNy()
    Nz = node.getNz()

    NxM = conf.NxMesh
    NyM = conf.NyMesh
    NzM = conf.NzMesh

    #print("filling ref")
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            for k in range(node.getNz()):
                #if n.getMpiGrid(i,j) == n.rank:
                if True:
                    for q in range(conf.NxMesh):
                        for r in range(conf.NyMesh):
                            for s in range(conf.NzMesh):
                                #print(i,j,k,q,r,s, "val=",val)
                                # insert 1/val to get more complex floats (instead of ints only)
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 0] = 1.0/val + 1.0 
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 1] = 1.0/val + 2.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 2] = 1.0/val + 3.0

                                data[i*NxM + q, j*NyM + r, k*NzM + s, 3] = 1.0/val + 4.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 4] = 1.0/val + 5.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 5] = 1.0/val + 6.0

                                data[i*NxM + q, j*NyM + r, k*NzM + s, 6] = 1.0/val + 7.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 7] = 1.0/val + 8.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 8] = 1.0/val + 9.0
                                val += 1
    return data


# fill Yee mesh with values
def fill_yee(node, data, conf):

    Nx = node.getNx()
    Ny = node.getNy()
    Nz = node.getNz()

    NxM = conf.NxMesh
    NyM = conf.NyMesh
    NzM = conf.NzMesh

    # lets put ref array into Yee lattice
    #print("filling yee")
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            for k in range(node.getNz()):
                #if n.getMpiGrid(i,j) == n.rank:
                if True:
                    c = node.getTile(i,j,k)
                    yee = c.getYee(0)

                    for q in range(conf.NxMesh):
                        for r in range(conf.NyMesh):
                            for s in range(conf.NzMesh):
                                yee.ex[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 0] 
                                yee.ey[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 1]
                                yee.ez[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 2]

                                yee.bx[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 3]
                                yee.by[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 4]
                                yee.bz[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 5]

                                yee.jx[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 6]
                                yee.jy[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 7]
                                yee.jz[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 8]


class IO(unittest.TestCase):

    def test_write_fields1D(self):

        ##################################################
        # write

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 1
        conf.Nz = 1
        conf.NxMesh = 2
        conf.NyMesh = 1
        conf.NzMesh = 1 
        conf.outdir = "io_test_1D/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        node = pycorgi.oneD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles1D(node, conf)

        ref = fill_ref(node, conf)
        fill_yee(node, ref, conf)

        pyplasmabox.vlv.oneD.writeYee(node, 0, conf.outdir)
        
        ##################################################
        # read using analysis tools

        yee = getYee(node, conf)

        arrs = combine_tiles(conf.outdir+"fields-0_0.h5", "ex", conf)

        Nx = node.getNx()
        Ny = node.getNy()
        Nz = node.getNz()

        NxM = conf.NxMesh
        NyM = conf.NyMesh
        NzM = conf.NzMesh
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    #if n.getMpiGrid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    self.assertAlmostEqual( arrs[i*NxM + q, j*NyM + r, k*NzM + s],  ref[i*NxM + q, j*NyM + r, k*NzM + s, 0], places=6)

        ##################################################

    def test_write_fields2D(self):

        ##################################################
        # write

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 4
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 6
        conf.NzMesh = 1 
        conf.outdir = "io_test_2D/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(node, conf)

        ref = fill_ref(node, conf)
        fill_yee(node, ref, conf)

        pyplasmabox.vlv.twoD.writeYee(node, 0, conf.outdir)
        
        ##################################################
        # read using analysis tools

        yee = getYee2D(node, conf)

        arrs = combine_tiles(conf.outdir+"fields-0_0.h5", "ex", conf)

        Nx = node.getNx()
        Ny = node.getNy()
        Nz = node.getNz()

        NxM = conf.NxMesh
        NyM = conf.NyMesh
        NzM = conf.NzMesh
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    #if n.getMpiGrid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    self.assertAlmostEqual( arrs[i*NxM + q, j*NyM + r, k*NzM + s],  ref[i*NxM + q, j*NyM + r, k*NzM + s, 0], places=6)

        ##################################################









