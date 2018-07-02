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
from visualize import plot2dYee
from visualize import saveVisz

try:
    import matplotlib.pyplot as plt
except:
    pass


#make tests deterministic by fixing the RNG seed
np.random.seed(0)



def filler_no_velocity(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    #zz = xloc[2] + np.random.rand(1)
    zz = 0.0

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
    #zz = xloc[2] + np.random.rand(1)
    zz = 0.0

    ur = conf.vel
    uc = randab(0.0, 2.0*np.pi) 
    ux = ur*np.sin( uc )
    uy = ur*np.cos( uc ) 
    uz = 0.0

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0

def filler_xvel(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    #zz = xloc[2] + np.random.rand(1)
    zz = 0.0

    ux = randab(-conf.vel, conf.vel)
    uy = 0.0
    uz = 0.0

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0

def zero_field(x,y,z):
    return 0.0

def const_field(x, y, z):
    return 1.0

def linear_field(x, y, z):
    #print("x = {} y = {} z = {}".format(x,y,z))
    #return 1.0*x 
    #return 10.0*y 
    #return 100.0*z 
    #return 1.0*x + 10.0*y 
    return 1.0*x + 10.0*y + 100.0*z


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

                        # get x_i,j,k
                        xloc0 = spatialLoc(node, (i,j), (l,  m,n), conf)

                        #get x_i+1/2, x_j+1/2, x_k+1/2
                        xloc1 = spatialLoc(node, (i,j), (l+1,m,n), conf)
                        yloc1 = spatialLoc(node, (i,j), (l,m+1,n), conf)
                        zloc1 = spatialLoc(node, (i,j), (l,m,n+1), conf)

                        # values in Yee lattice corners
                        xcor = xloc0[0]
                        ycor = xloc0[1]
                        zcor = xloc0[2]

                        # values in Yee lattice mids
                        xmid = 0.5*(xloc0[0] + xloc1[0])
                        ymid = 0.5*(xloc0[1] + yloc1[1])
                        zmid = 0.5*(xloc0[2] + zloc1[2])

                        #val = ffunc(xmid, ymid, zmid)

                        # enforce Yee lattice structure
                        yee.ex[l,m,n] = ffunc(xmid, ycor, zcor)
                        yee.ey[l,m,n] = ffunc(xcor, ymid, zcor)+1.0
                        #yee.ez[l,m,n] = ffunc(xcor, ycor, zmid)+2.0
                        yee.ez[l,m,n] = ffunc(xcor, ycor, zcor)+2.0  #2D hack

                        #yee.bx[l,m,n] = ffunc(xcor, ymid, zmid)+3.0
                        yee.bx[l,m,n] = ffunc(xcor, ymid, zcor)+3.0  #2D hack
                        #yee.by[l,m,n] = ffunc(xmid, ycor, zmid)+4.0 #2D hack
                        yee.by[l,m,n] = ffunc(xmid, ycor, zcor)+4.0
                        yee.bz[l,m,n] = ffunc(xmid, ymid, zcor)+5.0

                        yee.jx[l,m,n] = ffunc(xmid, ymid, zmid)
                        yee.jy[l,m,n] = ffunc(xmid, ymid, zmid)
                        yee.jz[l,m,n] = ffunc(xmid, ymid, zmid)

    


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

    zmin = 0.0
    zmax = 10.0

    cfl = 0.45
    c_omp = 10.0
    ppc = 1


    #dx = 1.0
    #dy = 1.0
    #dz = 1.0

    me = 1
    mi = 1

    Nspecies = 1

    outdir = "out"

    vel = 0.1

    #def __init__(self):
    #    print("initialized...")

    #update bounding box sizes
    #
    # NOTE: NxMesh = 5 grid looks like this:
    #
    # xmin      xmax
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


class PIC(unittest.TestCase):


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
        conf.NxMesh = 2
        conf.NyMesh = 2

        conf.Nx = 3
        conf.Ny = 3
        conf.update_bbox()

        conf.vel = 0.1

        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, const_field)
        inject(node, filler, conf) #injecting plasma particles


        # push particles couple of times to make them leak into neighboring tiles
        pusher   = pypic.Pusher()
        comm     = pypic.Communicator()


        for lap in range(40):
            #plot2dParticles(axs[0], node, conf)
            #saveVisz(lap, node, conf)

            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    pusher.solve(cell)

            #update particle boundaries
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    comm.check_outgoing_particles(cell)

            #copy particles
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    comm.get_incoming_particles(cell, node)

            #delete transferred particles
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    comm.delete_transferred_particles(cell)


        # count how many particles we now have
        n_particles = 0
        for i in range(conf.Nx):
            for j in range(conf.Ny):
                for k in range(conf.Nz):
                    cid = node.cellId(i,j)
                    c = node.getCellPtr(cid)

                    #print("({},{},{}) has {}".format(i,j,k,len(c.container.loc(0))))
                    n_particles += len(c.container.loc(0))

                    #self.assertTrue( 0.0 <= c.container.loc(0) <= conf.xmax )
                    #self.assertTrue( 0.0 <= c.container.loc(1) <= conf.ymax )
                    #self.assertTrue( 0.0 <= c.container.loc(2) <= conf.zmax )

                    for prtcl in range(len(c.container.loc(0))):
                        #print("{} {} {} maxs {} {} {}".format( 
                        #c.container.loc(0)[prtcl], 
                        #c.container.loc(1)[prtcl], 
                        #c.container.loc(2)[prtcl], 
                        #conf.xmax, conf.ymax, conf.zmax))

                        #print("prtcl {} x={} y={} z={} vx={} vy={} vz={}".format(
                        #    prtcl, 
                        #    c.container.loc(0)[prtcl],
                        #    c.container.loc(1)[prtcl],
                        #    c.container.loc(2)[prtcl],
                        #    c.container.vel(0)[prtcl],
                        #    c.container.vel(1)[prtcl],
                        #    c.container.vel(2)[prtcl]))

                        # check location
                        self.assertTrue( 0.0 <= c.container.loc(0)[prtcl] <= conf.xmax )
                        self.assertTrue( 0.0 <= c.container.loc(1)[prtcl] <= conf.ymax )
                        self.assertTrue( 0.0 <= c.container.loc(2)[prtcl] <= conf.zmax )

                        # check velocity 
                        velx = c.container.vel(0)[prtcl]
                        vely = c.container.vel(1)[prtcl]
                        velz = c.container.vel(2)[prtcl]
                        vel = np.sqrt( velx*velx + vely*vely + velz*velz )
                        self.assertAlmostEqual( vel, conf.vel, places=6 )

        tot_particles = (conf.Nx*conf.NxMesh *
                        conf.Ny*conf.NyMesh *
                        conf.Nz*conf.NzMesh *
                        conf.ppc)

        #tot_particles =(conf.NxMesh *
        #                conf.NyMesh *
        #                conf.NzMesh *
        #                conf.ppc)


        # assert that there is equal number of particles as we began with
        self.assertEqual( tot_particles, n_particles )



    def test_const_field_interpolation(self):

        conf = Conf()
        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, const_field)
        inject(node, filler_no_velocity, conf) #injecting plasma particles

        ##update boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.updateBoundaries2D(node)

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
                ey = c.container.ey()
                ez = c.container.ez()

                bx = c.container.bx()
                by = c.container.by()
                bz = c.container.bz()

                for i, x in enumerate(xx):
                    #print(i)
                    ex_ref = 1.0
                    ey_ref = 2.0
                    ez_ref = 3.0
                    self.assertEqual(ex[i], ex_ref)
                    self.assertEqual(ey[i], ey_ref)
                    self.assertEqual(ez[i], ez_ref)

                    bx_ref = 4.0
                    by_ref = 5.0
                    bz_ref = 6.0
                    self.assertEqual(bx[i], bx_ref)
                    self.assertEqual(by[i], by_ref)
                    self.assertEqual(bz[i], bz_ref)


    def test_linear_field_interpolation(self):

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3

        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, linear_field)
        inject(node, filler_no_velocity, conf) #injecting plasma particles

        ##update boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.updateBoundaries2D(node)

        #interpolate fields
        fintp = pypic.ParticleFieldInterpolator()
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                fintp.solve(cell)


        #test results; avoid boundaries because they are cyclic
        for i in range(1,conf.Nx-1):
            for j in range(1,conf.Ny-1):

                cid = node.cellId(i,j)
                c = node.getCellPtr(cid)

                xx = c.container.loc(0)
                yy = c.container.loc(1)
                zz = c.container.loc(2)

                ex = c.container.ex()
                ey = c.container.ey()
                ez = c.container.ez()

                bx = c.container.bx()
                by = c.container.by()
                bz = c.container.bz()

                for i, x in enumerate(xx):
                    #print(i)

                    ref = linear_field(xx[i], yy[i], zz[i]) #exact/true solution
                    ex_ref = ref + 0.0
                    ey_ref = ref + 1.0
                    ez_ref = ref + 2.0

                    #print("asserting {} {} {} vs {}".format(xx[i], yy[i], zz[i], ref))
                    self.assertAlmostEqual(ex[i], ex_ref, places=5)
                    self.assertAlmostEqual(ey[i], ey_ref, places=5)
                    self.assertAlmostEqual(ez[i], ez_ref, places=5)

                    bx_ref = ref + 3.0
                    by_ref = ref + 4.0
                    bz_ref = ref + 5.0
                    self.assertAlmostEqual(bx[i], bx_ref, places=5)
                    self.assertAlmostEqual(by[i], by_ref, places=5)
                    self.assertAlmostEqual(bz[i], bz_ref, places=5)



    def test_filters(self):
        """ filter integration test with rest of the PIC functions"""


        plt.fig = plt.figure(1, figsize=(5,7))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(6, 1)
        
        axs = []
        for ai in range(6):
            axs.append( plt.subplot(gs[ai]) )


        conf = Conf()
        conf.Nx = 10
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 5
        conf.NzMesh = 1
        conf.ppc = 1
        conf.vel = 0.1
        conf.update_bbox()

        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, zero_field)
        inject(node, filler, conf) #injecting plasma particles
        #inject(node, filler_xvel, conf) #injecting plasma particles


        #pusher   = pypic.Pusher()
        #fintp    = pypic.ParticleFieldInterpolator()
        #comm     = pypic.Communicator()
        currint  = pypic.Depositer()
        analyzer = pypic.Analyzator()
        flt     =  pypic.Filter(conf.NxMesh, conf.NyMesh)

        flt.init_gaussian_kernel(2.0, 2.0)

        #for lap in range(0, conf.Nt):
        for lap in [0]:

            #analyze
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    analyzer.analyze(cell)

            ##update boundaries
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    cell.updateBoundaries2D(node)

            #deposit current
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    currint.deposit(cell)

            #exchange currents
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    cell.exchangeCurrents2D(node)

            plot2dParticles(axs[0], node, conf, downsample=0.01)
            plot2dYee(axs[1], node, conf, 'rho')
            plot2dYee(axs[2], node, conf, 'jx')
            plot2dYee(axs[3], node, conf, 'jy')
            plot2dYee(axs[4], node, conf, 'jz')

            #filter
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    flt.get_padded_current(cell, node)
                    flt.fft_image_forward()
                    flt.apply_kernel()
                    flt.fft_image_backward()
                    flt.set_current(cell)

            #cycle new and temporary currents
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    cell.cycleCurrent2D()

            plot2dYee(axs[5], node, conf, 'jx')


            saveVisz(lap, node, conf)

