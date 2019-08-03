from mpi4py import MPI

import unittest

import sys
import numpy as np

from math import floor, ceil
#from scipy.signal import convolve2d
#from scipy.signal import convolve

import pycorgi
import pyrunko.pic.twoD as pypic
import pyrunko.tools.twoD as pytools
import pyrunko.fields.twoD as pyfields

from initialize_pic import loadTiles
from initialize_pic import spatialLoc
from injector_pic import inject


def const_field(x, y, z):
    return 1.0

def linear_ramp(x,y,z):
    return x + y + z


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(grid, conf, ffunc):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            c = grid.get_tile(i,j)
            yee = c.get_yee(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        #get x_i+1/2 (Yee lattice so rho_i)
                        xloc0 = spatialLoc(grid, (i,j), (l,  m,n), conf)
                        xloc1 = spatialLoc(grid, (i,j), (l+1,m,n), conf)

                        xmid = 0.5*(xloc0[0] + xloc1[0])
                        ymid = 0.5*(xloc0[1] + xloc1[1])
                        zmid = 0.5*(xloc0[2] + xloc1[2])

                        val = ffunc(xmid, ymid, zmid)

                        yee.jx[l,m,n] = val
                        yee.jy[l,m,n] = val+1.0
                        yee.jz[l,m,n] = val+2.0


# basic Conf file/class for PiC simulation testing
class Conf:

    Nx = 1
    Ny = 1
    Nz = 1

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

    #def __init__(self):
    #    print("initialized...")

    #update bounding box sizes
    def update_bbox(self):
        self.xmin = 0.0
        self.xmax = self.Nx*self.NxMesh

        self.ymin = 0.0
        self.ymax = self.Ny*self.NyMesh


#three point digital filter
def digi3(i,j):
    digi3a = np.array([[ 4., 2.],
                       [ 2., 1.]])
    i2 = abs(i)
    j2 = abs(j)
    if (i2 > 1) or (j2 > 1):
        return 0.0
    else:
        return (0.25/3.0)*digi3a[i2,j2]

#five point digital filter
def digi5(i,j):
    digi5a = np.array([[ 6., 4., 1],
                       [ 4., 0., 0],
                       [ 1., 0., 0]])
    i2 = abs(i)
    j2 = abs(j)
    if (i2 > 2) or (j2 > 2):
        return 0.0
    else:
        return (1.0/16.0)*digi5a[i2,j2]


class FilterTests(unittest.TestCase):

    #test loading of different filters
    def test_init(self):
        conf = Conf()

        flt1 = pyfields.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)




