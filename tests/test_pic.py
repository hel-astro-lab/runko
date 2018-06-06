import unittest

import sys
sys.path.append('python')
import numpy as np

import corgi
import pyplasma as plasma
import pypic 

sys.path.append('pic')
from pic import loadCells
from pic import inject
from pic import spatialLoc


from visualize_pic import Particles
from visualize_pic import plot2dParticles
from visualize import saveVisz

try:
    import matplotlib.pyplot as plt
except:
    pass


def filler_no_velocity(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    zz = xloc[2] + np.random.rand(1)

    x0 = [xx, yy, zz]
    u0 = [0.0, 0.0, 0.0]
    return x0, u0


# random number between [xmin, xmax]
def randab(xmin, xmax):
    return xmin + (xmax-xmin)*np.random.rand(1)


def filler(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    zz = xloc[2] + np.random.rand(1)

    ux = randab(-1.0, 1.0)
    uy = randab(-1.0, 1.0)
    uz = randab(-1.0, 1.0)

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0

def const_field(x, y, z):
    return 1.0


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(node, conf, ffunc):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            c = node.getCellPtr(i,j)
            yee = c.getYee(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        #get x_i+1/2 (Yee lattice so rho_i)
                        xloc0 = spatialLoc(node, (i,j), (l,  m,n), conf)
                        xloc1 = spatialLoc(node, (i,j), (l+1,m,n), conf)

                        xmid = 0.5*(xloc0[0] + xloc1[0])
                        ymid = 0.5*(xloc0[1] + xloc1[1])
                        zmid = 0.5*(xloc0[2] + xloc1[2])

                        val = ffunc(xmid, ymid, zmid)

                        yee.ex[l,m,n] = 1.0
                        yee.ey[l,m,n] = 1.0
                        yee.ez[l,m,n] = 1.0

                        yee.bx[l,m,n] = 1.0
                        yee.by[l,m,n] = 1.0
                        yee.bz[l,m,n] = 1.0



# basic Conf file/class for PiC simulation testing
class Conf:

    Nx = 1
    Ny = 1
    Nz = 1

    NxMesh = 4
    NyMesh = 4
    NzMesh = 1

    xmin = 0.0
    xmax = 10.0

    ymin = 0.0
    ymax = 10.0

    cfl = 0.45
    c_omp = 10.0
    ppc = 2

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



class PIC(unittest.TestCase):


    def test_communication(self):

        plt.fig = plt.figure(1, figsize=(3,3))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(1, 1)
        
        axs = []
        for ai in range(1):
            axs.append( plt.subplot(gs[ai]) )



        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3
        conf.update_bbox()

        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, const_field)
        inject(node, filler, conf) #injecting plasma particles


        # push particles couple of times to make them leak into neighboring tiles
        pusher   = pypic.Pusher()
        comm     = pypic.Communicator()


        for lap in range(1):
            plot2dParticles(axs[0], node, conf)
            saveVisz(lap, node, conf)

            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    pusher.solve(cell)

            ##update particle boundaries
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    comm.check_outgoing_particles(cell)




        # count how many particles we now have
        n_particles = 0
        for i in range(conf.Nx):
            for j in range(conf.Ny):
                for k in range(conf.Nz):
                    cid = node.cellId(i,j)
                    c = node.getCellPtr(cid)

                    n_particles += len(c.container.loc(0))

        tot_particles = (conf.Nx*conf.NxMesh *
                        conf.Ny*conf.NyMesh *
                        conf.Nz*conf.NzMesh *
                        conf.ppc)


        # assert that there is equal number of particles as we began with
        self.assertEqual( tot_particles, n_particles )



    def field_interpolation(self):

        conf = Conf()
        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, const_field)
        inject(node, filler_no_velocity, conf) #injecting plasma particles

        #interpolate fields
        fintp = pypic.ParticleFieldInterpolator()
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                fintp.solve(cell)

        #test results
        for i in range(conf.Nx):
            for j in range(conf.Ny):
                cid = node.cellId(i,j)
                c = node.getCellPtr(cid)

                xx = c.container.loc(0)
                yy = c.container.loc(1)
                zz = c.container.loc(2)

                ux = c.container.vel(0)
                uy = c.container.vel(1)
                uz = c.container.vel(2)

                ex = c.container.ex()
                #ey = c.container.ey()
                #ez = c.container.ez()

                #bx = c.container.bx()
                #by = c.container.by()
                #bz = c.container.bz()

                for i, x in enumerate(xx):
                    ex_ref = 1.0
                    self.assertEqual(ex[i], ex_ref)



