import unittest

import sys
import numpy as np

import pyplasmabox.rad as pyrad


from initialize_rad import bbodySample
from initialize_rad import rand3Dloc
from initialize_rad import rand3Dvel


try:
    import matplotlib.pyplot as plt
except:
    pass


#make tests deterministic by fixing the RNG seed
np.random.seed(0)



# basic Conf file/class for PiC simulation testing
class Conf:

    Nx = 1
    Ny = 1
    Nz = 1

    NxMesh = 10
    NyMesh = 10
    NzMesh = 1

    xmin = 0.0
    xmax = 10.0

    ymin = 0.0
    ymax = 10.0

    zmin = 0.0
    zmax = 10.0

    cfl = 0.45
    c_omp = 10.0
    ppc = 1

    dx = 1.0
    dy = 1.0
    dz = 1.0

    Nspecies = 1

    outdir = "out"

    #def __init__(self):
    #    print("initialized...")

    #update bounding box sizes
    #
    # NOTE: NxMesh = 5 grid looks like this:
    #
    # xmin      xmax
    #  |         |
    #  v         v
    #  |_|_|_|_|_
    #  0 1 2 3 4 5
    #
    def update_bbox(self):
        self.xmin = 0.0
        self.xmax = self.Nx*self.NxMesh

        self.ymin = 0.0
        self.ymax = self.Ny*self.NyMesh

        self.zmin = 0.0
        self.zmax = self.Nz*self.NzMesh


class RAD(unittest.TestCase):

    def test_communication(self):

        #plt.fig = plt.figure(1, figsize=(3,3))
        #plt.rc('font', family='serif', size=12)
        #plt.rc('xtick')
        #plt.rc('ytick')
        #
        #gs = plt.GridSpec(1, 1)
        #
        #axs = []
        #for ai in range(1):
        #    axs.append( plt.subplot(gs[ai]) )

        conf = Conf()
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.Nx = 3
        conf.Ny = 3
        conf.Ny = 1
        conf.ppc = 1
        conf.update_bbox()

        kT = 1.0 #black body photon field temperature

        Nprtcls = conf.NxMesh * conf.NyMesh * conf.NzMesh * conf.ppc

        container = pyrad.PhotonBlock(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        container.reserve(Nprtcls)

        weight = 1.0

        for ip in range(Nprtcls):
            ene = bbodySample(kT)

            x0 = rand3Dloc(conf)
            u0 = rand3Dvel(1.0)
        
            container.add_particle(x0, u0, weight, ene)
        



