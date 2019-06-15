from mpi4py import MPI
import unittest

import sys
import numpy as np

import pycorgi
import pyrunko.rad as pyrad


import initialize_pic as init_pic
import initialize_rad as init_rad

from initialize_rad import bbodySample
from initialize_rad import rand3Dloc
from initialize_rad import rand3Dvel


do_plots = True

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
    ppc = 1 #particles per cell
    ppt = 1 #photons per tile

    dx = 1.0
    dy = 1.0
    dz = 1.0

    gamma_e = 0.0
    gamma_i = 0.0

    me = 1
    mi = 1

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


class radiation(unittest.TestCase):

    def test_initialization(self):

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

        container = pyrad.PhotonContainer()
        container.reserve(Nprtcls)

        weight = 1.0

        ene_ref = np.zeros(Nprtcls)
        wgt_ref = np.zeros(Nprtcls)
        x0_ref  = np.zeros((Nprtcls,3))
        u0_ref  = np.zeros((Nprtcls,3))

        for ip in range(Nprtcls):
            ene = bbodySample(kT)

            x0 = rand3Dloc(conf)
            u0 = rand3Dvel(1.0)
        
            container.add_particle(x0, u0, weight, ene)
        
            # add also to reference array
            ene_ref[ip]  = ene
            wgt_ref[ip]  = weight
            x0_ref[ip,:] = x0
            u0_ref[ip,:] = u0


        ene  = container.ene()
        wgt  = container.wgt()

        loc0 = container.loc(0)
        loc1 = container.loc(1)
        loc2 = container.loc(2)

        vel0 = container.vel(0)
        vel1 = container.vel(1)
        vel2 = container.vel(2)

        for ip in range(Nprtcls):
            self.assertAlmostEqual( container.ene()[ip],  ene_ref[ip],  places=5)
            self.assertAlmostEqual( container.wgt()[ip],  wgt_ref[ip],  places=5)

            self.assertAlmostEqual( container.loc(0)[ip], x0_ref[ip,0], places=5)
            self.assertAlmostEqual( container.loc(1)[ip], x0_ref[ip,1], places=5)
            self.assertAlmostEqual( container.loc(2)[ip], x0_ref[ip,2], places=5)

            self.assertAlmostEqual( container.vel(0)[ip], u0_ref[ip,0], places=5)
            self.assertAlmostEqual( container.vel(1)[ip], u0_ref[ip,1], places=5)
            self.assertAlmostEqual( container.vel(2)[ip], u0_ref[ip,2], places=5)


    def test_blackbody(self):

        if do_plots:
            try:
                plt.fig = plt.figure(1, figsize=(3,3))
                plt.rc('font', family='serif', size=12)
                plt.rc('xtick')
                plt.rc('ytick')
                
                gs = plt.GridSpec(1, 1)
                 
                axs = []
                for ai in range(1):
                    axs.append( plt.subplot(gs[ai]) )
            except:
                pass

        conf = Conf()
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.ppc = 1
        conf.update_bbox()

        kTbb = 1.5e-2/511.0 # Photon blackbody temperature in units of mc^2; 2e-4 = 0.1 keV

        #Nprtcls = conf.NxMesh * conf.NyMesh * conf.NzMesh * conf.ppc
        Nprtcls = int(1e4)

        container = pyrad.PhotonContainer()
        container.reserve(Nprtcls)

        weight = 1.0
        for ip in range(Nprtcls):
            ene = bbodySample(kTbb)

            x0 = rand3Dloc(conf)
            u0 = rand3Dvel(1.0)
        
            container.add_particle(x0, u0, weight, ene)
        

        nbins = 100
        emin = 1.0e-3
        emax = 2.0e+3
        #Nph = container.size()
        Nph = Nprtcls
        self.assertEqual(Nprtcls, Nph)

        if do_plots:
            try:
                axs[0].set_xlim(emin, emax)
                #axs[0].set_ylim(1e-4, 1.2e-0)
                axs[0].minorticks_on()
                axs[0].set_xscale('log')
                axs[0].set_yscale('log')
            except:
                pass

        ene = np.array(container.ene())
        phhist, edges = np.histogram(ene*511., np.logspace(np.log10(emin), np.log10(emax), nbins))

        if do_plots:
            try:
                axs[0].bar(edges[:-1], phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
            except:
                pass

        prtcl_sum=np.sum(phhist*edges[:-1])/Nph
        print("Energy of photons: {}".format(prtcl_sum))

        bbrad_sum = np.sum(3.0e12*2.*edges**3/(np.exp(edges/kTbb)-1.0))
        print("Blackbody energy: {}".format(bbrad_sum))


        try:
            plt.savefig("blackbody.pdf")
        except:
            pass


    
    # test that tile supports both pic and rad initialization
    def test_initialization(self):

        conf = Conf()
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 1
        conf.Nx = 1
        conf.Ny = 1
        conf.Nz = 1
        conf.ppc = 1
        conf.update_bbox()

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        c = pyrad.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        init_pic.initialize_tile(c, 0, 0, grid, conf)
        init_rad.initialize_tile(c, 0, 0, grid, conf)



