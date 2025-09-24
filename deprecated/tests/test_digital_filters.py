from mpi4py import MPI

import unittest

import sys
import numpy as np

from math import floor, ceil

# from scipy.signal import convolve2d
# from scipy.signal import convolve

import pycorgi
import pyrunko.pic as pypic
import pyrunko.tools as pytools
import pyrunko.emf as pyfields


def const_field(x, y, z):
    return 1.0, 2.0, 3.0


def linear_ramp(x, y, z):
    return x + y + z, 2 * (x + y + z), 3 * (x + y + z)


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em_tile(tile, conf, ffunc):

    gs = tile.get_grids(0)
    for l in range(conf.NxMesh):
        for m in range(conf.NyMesh):
            for n in range(conf.NzMesh):

                # use indexes directly as coordinates
                valx, valy, valz = ffunc(l, m, n)

                gs.jx[l, m, n] = valx
                gs.jy[l, m, n] = valy
                gs.jz[l, m, n] = valz


# basic Conf file/class for PiC simulation testing
class Conf:

    Nx = 1
    Ny = 1
    Nz = 1

    oneD = False
    twoD = False
    threeD = False

    NxMesh = 5
    NyMesh = 5
    NzMesh = 1

    xmin = 0.0
    xmax = 10.0

    ymin = 0.0
    ymax = 10.0

    cfl = 0.45
    c_omp = 10.0
    ppc = 1

    gamma_e = 0.0
    gamma_i = 0.0

    dx = 1.0
    dy = 1.0
    dz = 1.0

    me = 1
    mi = 1

    Nspecies = 1

    outdir = "out"

    qe = 1.0

    prtcl_types = ['e-', 'e+']

    # def __init__(self):
    #    print("initialized...")

    # update bounding box sizes
    def update_bbox(self):
        self.xmin = 0.0
        self.xmax = self.Nx * self.NxMesh

        self.ymin = 0.0
        self.ymax = self.Ny * self.NyMesh


# three point digital filter
def digi3(i, j):
    digi3a = np.array([[4.0, 2.0], [2.0, 1.0]])
    i2 = abs(i)
    j2 = abs(j)
    if (i2 > 1) or (j2 > 1):
        return 0.0
    else:
        return (0.25 / 3.0) * digi3a[i2, j2]


# five point digital filter
def digi5(i, j):
    digi5a = np.array([[6.0, 4.0, 1], [4.0, 0.0, 0], [1.0, 0.0, 0]])
    i2 = abs(i)
    j2 = abs(j)
    if (i2 > 2) or (j2 > 2):
        return 0.0
    else:
        return (1.0 / 16.0) * digi5a[i2, j2]


def get_js(tile, conf):
    jx = np.zeros((conf.NxMesh, conf.NyMesh, conf.NzMesh))
    jy = np.zeros((conf.NxMesh, conf.NyMesh, conf.NzMesh))
    jz = np.zeros((conf.NxMesh, conf.NyMesh, conf.NzMesh))

    gs = tile.get_grids()

    for q in range(conf.NxMesh):
        for r in range(conf.NyMesh):
            for s in range(conf.NzMesh):
                jx[q, r, s] = gs.jx[q, r, s]
                jy[q, r, s] = gs.jy[q, r, s]
                jz[q, r, s] = gs.jz[q, r, s]

    return jx, jy, jz


class FilterTests(unittest.TestCase):

    # test loading of different filters
    def test_init(self):
        conf = Conf()
        conf.NxMesh = 10
        conf.NyMesh = 15
        conf.NzMesh = 1

        conf.oneD  = False
        conf.twoD  = True
        conf.threeD= False

        flts = []
        flts.append( pyfields.twoD.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        flts.append( pyfields.twoD.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        flts.append( pyfields.twoD.General3pStrided(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        flts.append( pyfields.twoD.Binomial2Strided2(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        flts.append( pyfields.twoD.Compensator2(conf.NxMesh, conf.NyMesh, conf.NzMesh) )


        tile = pyfields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        insert_em_tile(tile, conf, const_field)

        for flt in flts:
            flt.solve(tile)


    def test_normalization(self):

        conf = Conf()
        conf.twoD = True

        tile = pyfields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        insert_em_tile(tile, conf, const_field)

        # get current arrays and sum (only internal values #to avoid garbage on boundaries

        jx, jy, jz = get_js(tile, conf)
        sumx = np.sum(jx[1:-1, 1:-1, 0])
        sumy = np.sum(jy[1:-1, 1:-1, 0])
        sumz = np.sum(jz[1:-1, 1:-1, 0])

        # should be 25,50,75 (now 9,12,27 without boundaries)
        # print("sumx {} sumy {} sumz {}".format(sumx, sumy, sumz))
        # print(jx[:,:,0])
        # print(jy[:,:,0])
        # print(jz[:,:,0])

        flt = pyfields.twoD.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        flt.solve(tile)

        jx1, jy1, jz1 = get_js(tile, conf)
        sumx1 = np.sum(jx1[1:-1, 1:-1, 0])
        sumy1 = np.sum(jy1[1:-1, 1:-1, 0])
        sumz1 = np.sum(jz1[1:-1, 1:-1, 0])

        # print(jx1[:,:,0])
        # print(jy1[:,:,0])
        # print(jz1[:,:,0])
        # print("sumx {} sumy {} sumz {}".format(sumx1, sumy1, sumz1))

        self.assertAlmostEqual(sumx, sumx1, places=5)
        self.assertAlmostEqual(sumy, sumy1, places=5)
        self.assertAlmostEqual(sumz, sumz1, places=5)

    def test_init3D(self):

        conf = Conf()
        conf.NxMesh = 10
        conf.NyMesh = 15
        conf.NzMesh = 20

        conf.oneD = False
        conf.twoD = False
        conf.threeD = True

        flts = []
        flts.append( pyfields.threeD.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        #flts.append( pyfields.threeD.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        #flts.append( pyfields.threeD.General3pStrided(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        #flts.append( pyfields.threeD.Binomial2Strided2(conf.NxMesh, conf.NyMesh, conf.NzMesh) )
        #flts.append( pyfields.threeD.Compensator2(conf.NxMesh, conf.NyMesh, conf.NzMesh) )


        tile = pyfields.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        insert_em_tile(tile, conf, const_field)

        for flt in flts:
            flt.solve(tile)


    def test_normalization3D(self):

        conf = Conf()
        conf.threeD = True
        conf.NxMesh = 6
        conf.NyMesh = 7
        conf.NzMesh = 8

        tile = pyfields.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        insert_em_tile(tile, conf, const_field)

        # get current arrays and sum (only internal values #to avoid garbage on boundaries

        jx, jy, jz = get_js(tile, conf)
        st = 2
        sumx = np.sum(jx[st:-st, st:-st, st:-st])
        sumy = np.sum(jy[st:-st, st:-st, st:-st])
        sumz = np.sum(jz[st:-st, st:-st, st:-st])

        # should be 25,50,75 (now 9,12,27 without boundaries)
        # print("sumx {} sumy {} sumz {}".format(sumx, sumy, sumz))
        # print(jx[:,:,0])
        # print(jy[:,:,0])
        # print(jz[:,:,0])

        flt = pyfields.threeD.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        flt.solve(tile)

        jx1, jy1, jz1 = get_js(tile, conf)
        sumx1 = np.sum(jx1[st:-st, st:-st, st:-st])
        sumy1 = np.sum(jy1[st:-st, st:-st, st:-st])
        sumz1 = np.sum(jz1[st:-st, st:-st, st:-st])

        # print(jx1[:,:,0])
        # print(jy1[:,:,0])
        # print(jz1[:,:,0])
        # print("sumx {} sumy {} sumz {}".format(sumx1, sumy1, sumz1))

        self.assertAlmostEqual(sumx, sumx1, places=5)
        self.assertAlmostEqual(sumy, sumy1, places=5)
        self.assertAlmostEqual(sumz, sumz1, places=5)
