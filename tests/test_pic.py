from mpi4py import MPI
import unittest

import sys
import numpy as np

import pycorgi
import pyplasmabox.pic.twoD as pypic
import pyplasmabox.tools.twoD as pytools



from initialize_pic import loadTiles
from initialize_pic import spatialLoc
from injector_pic import inject


from visualize_pic import Particles
from visualize_pic import plot2dParticles
from visualize import plotNode
from visualize import plot2dYee
from visualize import getYee2D
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
    zz = 0.5

    ur = conf.vel
    uc = randab(0.0, 2.0*np.pi) 
    us = randab(0.0, 1.0*np.pi) 

    #3D
    #ux = ur*np.sin( uc )*np.sin(us)
    #uy = ur*np.cos( uc )*np.sin(us)
    #uz = ur*np.cos(us)

    #2D
    ux = ur*np.sin( uc )
    uy = ur*np.cos( uc )
    uz = 0.5

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]


    return x0, u0

def filler_xvel(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    #zz = xloc[2] + np.random.rand(1)

    #xx = xloc[0] #+ 0.5
    #yy = xloc[1] #+ 0.5
    zz = 0.5

    x0 = [xx, yy, zz]

    #if (conf.NyMesh < yy < 2.0*conf.NyMesh) and (conf.NxMesh < xx < 2.0*conf.NxMesh):
    #    u0 = [0.0, 1.0, 0.0]
    #else:
    #    u0 = [0.0, 0.0, 0.0]

    u0 = [+1.0, +1.0, +1.0]
    #u0 = [-1.0, -1.0, -1.0]

    return x0, u0

def zero_field(x,y,z):
    return 0.0

def const_field(x, y, z):
    return 1.0

def linear_field(x, y, z):
    #print("x = {} y = {} z = {}".format(x,y,z))
    #return 1.0*x  
    return 1.0*x + 1.0*y
    #return 10.0*y 
    #return 100.0*z 
    #return 1.0*x + 10.0*y 
    return 1.0*x + 10.0*y + 100.0*z


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(node, conf, ffunc):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    for i in range(node.get_Nx()):
        for j in range(node.get_Ny()):
            c = node.get_tile(i,j)
            yee = c.get_yee(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        # get x_i,j,k
                        xloc0 = spatialLoc(node, (i,j), (l,m,n), conf)

                        #get x_i+1/2, x_j+1/2, x_k+1/2
                        xloc1 = spatialLoc(node, (i,j), (l+1,m,  n),   conf)
                        yloc1 = spatialLoc(node, (i,j), (l,  m+1,n),   conf)
                        zloc1 = spatialLoc(node, (i,j), (l,  m,  n+1), conf)

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

    gamma_e = 0.0
    gamma_i = 0.0

    dx = 1.0
    dy = 1.0
    dz = 1.0

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
        conf.NxMesh = 3
        conf.NyMesh = 3

        conf.Nx = 3
        conf.Ny = 3
        conf.Ny = 1
        conf.update_bbox()

        conf.vel = 0.1

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny, conf.Nz)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles(node, conf)
        insert_em(node, conf, const_field)
        inject(node, filler, conf) #injecting plasma particles


        # push particles couple of times to make them leak into neighboring tiles
        pusher   = pypic.BorisPusher()

        for lap in range(40):
            #plot2dParticles(axs[0], node, conf)
            #saveVisz(lap, node, conf)

            for cid in node.get_local_tiles():
                tile = node.get_tile(cid)
                pusher.solve(tile)

            ##################################################
            # communication

            #update particle boundaries
            for cid in node.get_local_tiles():
                tile = node.get_tile(cid)
                tile.check_outgoing_particles()

            # global mpi exchange (independent)
            for cid in node.get_boundary_tiles():
                tile = node.get_tile(cid)
                tile.pack_outgoing_particles()

            # MPI global exchange
            # transfer primary and extra data
            node.send_data(0) #(indepdendent)
            node.send_data(1) #(indepdendent)

            node.recv_data(0) #(indepdendent)
            node.recv_data(1) #(indepdendent)

            node.wait_data(0) #(indepdendent)
            node.wait_data(1) #(indepdendent)

            # global unpacking (independent)
            for cid in node.get_virtual_tiles(): 
                tile = node.get_tile(cid)
                tile.unpack_incoming_particles()
                tile.check_outgoing_particles()

            # transfer local + global
            for cid in node.get_local_tiles():
                tile = node.get_tile(cid)
                tile.get_incoming_particles(node)

            # delete local transferred particles
            for cid in node.get_local_tiles():
                tile = node.get_tile(cid)
                tile.delete_transferred_particles()

            for cid in node.get_virtual_tiles(): 
                tile = node.get_tile(cid)
                tile.delete_all_particles()

        # count how many particles we now have
        n_particles = 0
        for i in range(conf.Nx):
            for j in range(conf.Ny):
                for k in range(conf.Nz):
                    cid = node.id(i,j)
                    c = node.get_tile(cid)

                    container = c.get_container(0)
                    print("({},{},{}) has {}".format(i,j,k,len(container.loc(0))))
                    n_particles += len(container.loc(0))

                    #self.assertTrue( 0.0 <= container.loc(0) <= conf.xmax )
                    #self.assertTrue( 0.0 <= container.loc(1) <= conf.ymax )
                    #self.assertTrue( 0.0 <= container.loc(2) <= conf.zmax )

                    for prtcl in range(len(container.loc(0))):
                        #print("{} {} {} maxs {} {} {}".format( 
                        #container.loc(0)[prtcl], 
                        #container.loc(1)[prtcl], 
                        #container.loc(2)[prtcl], 
                        #conf.xmax, conf.ymax, conf.zmax))

                        #print("prtcl {} x={} y={} z={} vx={} vy={} vz={}".format(
                        #    prtcl, 
                        #    container.loc(0)[prtcl],
                        #    container.loc(1)[prtcl],
                        #    container.loc(2)[prtcl],
                        #    container.vel(0)[prtcl],
                        #    container.vel(1)[prtcl],
                        #    container.vel(2)[prtcl]))

                        # check location
                        self.assertTrue( 0.0 <= container.loc(0)[prtcl] <= conf.xmax )
                        self.assertTrue( 0.0 <= container.loc(1)[prtcl] <= conf.ymax )
                        self.assertTrue( 0.0 <= container.loc(2)[prtcl] <= conf.zmax )

                        # check velocity 
                        velx = container.vel(0)[prtcl]
                        vely = container.vel(1)[prtcl]
                        velz = container.vel(2)[prtcl]
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
        node = pycorgi.twoD.Node(conf.Nx, conf.Ny, conf.Nz)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles(node, conf)
        insert_em(node, conf, const_field)
        inject(node, filler_no_velocity, conf) #injecting plasma particles

        ##update boundaries
        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                tile = node.get_tile(i,j)
                tile.update_boundaries(node)

        #interpolate fields
        fintp = pypic.LinearInterpolator()
        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                tile = node.get_tile(i,j)
                fintp.solve(tile)

        #test results
        for i in range(conf.Nx):
            for j in range(conf.Ny):

                cid = node.id(i,j)
                c = node.get_tile(cid)
                container = c.get_container(0)

                xx = container.loc(0)
                yy = container.loc(1)
                zz = container.loc(2)

                ux = container.vel(0)
                uy = container.vel(1)
                uz = container.vel(2)

                for i, x in enumerate(xx):
                    #print(i)
                    ex_ref = 1.0
                    ey_ref = 2.0
                    ez_ref = 3.0
                    self.assertEqual(container.ex(i), ex_ref)
                    self.assertEqual(container.ey(i), ey_ref)
                    self.assertEqual(container.ez(i), ez_ref)

                    bx_ref = 4.0
                    by_ref = 5.0
                    bz_ref = 6.0
                    self.assertEqual(container.bx(i), bx_ref)
                    self.assertEqual(container.by(i), by_ref)
                    self.assertEqual(container.bz(i), bz_ref)


    def test_linear_field_interpolation(self):

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny, conf.Nz)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles(node, conf)
        insert_em(node, conf, linear_field)
        inject(node, filler_no_velocity, conf) #injecting plasma particles

        ##update boundaries
        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                tile = node.get_tile(i,j)
                tile.update_boundaries(node)

        #interpolate fields
        fintp = pypic.LinearInterpolator()
        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                tile = node.get_tile(i,j)
                fintp.solve(tile)


        #test results; avoid boundaries because they are cyclic
        for i in range(1,conf.Nx-1):
            for j in range(1,conf.Ny-1):

                cid = node.id(i,j)
                c = node.get_tile(cid)
                container = c.get_container(0)

                xx = container.loc(0)
                yy = container.loc(1)
                zz = container.loc(2)

                for i, x in enumerate(xx):
                    #print(i)

                    ref = linear_field(xx[i], yy[i], zz[i]) #exact/true solution
                    ex_ref = ref + 0.0
                    ey_ref = ref + 1.0
                    ez_ref = ref + 2.0

                    #print("asserting {} {} {} vs {}".format(xx[i], yy[i], zz[i], ref))
                    self.assertAlmostEqual(container.ex(i), ex_ref, places=5)
                    self.assertAlmostEqual(container.ey(i), ey_ref, places=5)
                    self.assertAlmostEqual(container.ez(i), ez_ref, places=5)

                    bx_ref = ref + 3.0
                    by_ref = ref + 4.0
                    bz_ref = ref + 5.0
                    self.assertAlmostEqual(container.bx(i), bx_ref, places=5)
                    self.assertAlmostEqual(container.by(i), by_ref, places=5)
                    self.assertAlmostEqual(container.bz(i), bz_ref, places=5)



    def test_filters(self):
        """ filter integration test with rest of the PIC functions"""

        try:
            plt.fig = plt.figure(1, figsize=(5,7))
            plt.rc('font', family='serif', size=12)
            plt.rc('xtick')
            plt.rc('ytick')

            gs = plt.GridSpec(8, 1)

            axs = []
            for ai in range(8):
                axs.append( plt.subplot(gs[ai]) )
        except:
            pass


        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 1
        conf.ppc = 10
        conf.vel = 0.1
        conf.update_bbox()

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny, conf.Nz)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles(node, conf)
        #insert_em(node, conf, linear_field)
        insert_em(node, conf, zero_field)
        inject(node, filler, conf) #injecting plasma particles
        #inject(node, filler_xvel, conf) #injecting plasma particles

        #pusher   = pypic.BorisPusher()
        #fintp    = pypic.LinearInterpolator()
        currint  = pypic.ZigZag()
        analyzer = pypic.Analyzator()

        flt =  pytools.Filter(conf.NxMesh, conf.NyMesh)
        flt.init_gaussian_kernel(1.0, 1.0)


        #for lap in range(0, conf.Nt):
        for lap in range(1):

            #analyze
            for j in range(node.get_Ny()):
                for i in range(node.get_Nx()):
                    tile = node.get_tile(i,j)
                    analyzer.analyze2d(tile)

            #update boundaries
            for j in range(node.get_Ny()):
                for i in range(node.get_Nx()):
                    tile = node.get_tile(i,j)
                    tile.update_boundaries(node)

            #deposit current
            for j in range(node.get_Ny()):
                for i in range(node.get_Nx()):
                    tile = node.get_tile(i,j)
                    currint.solve(tile)

            #exchange currents
            for j in range(node.get_Ny()):
                for i in range(node.get_Nx()):
                    tile = node.get_tile(i,j)
                    tile.exchange_currents(node)

            try:
                plotNode(axs[0], node, conf)
                #plot2dParticles(axs[0], node, conf, downsample=0.1)
                plot2dYee(axs[1], node, conf, 'rho')
                plot2dYee(axs[2], node, conf, 'jx')
                plot2dYee(axs[3], node, conf, 'jy')
                plot2dYee(axs[4], node, conf, 'jz')
            except:
                pass

            yee_ref = getYee2D(node, conf)

            #filter
            for j in range(node.get_Ny()):
                for i in range(node.get_Nx()):
                    print(" i j ({},{})".format(i,j))
                    tile = node.get_tile(i,j)
                    flt.get_padded_current(tile, node)

                    # fourier space filtering
                    flt.fft_image_forward()
                    flt.apply_kernel()
                    flt.fft_image_backward()

                    # direct filtering
                    #for fj in range(1):
                    #    flt.direct_convolve_3point()
                    flt.set_current(tile)

            #cycle new and temporary currents
            for j in range(node.get_Ny()):
                for i in range(node.get_Nx()):
                    tile = node.get_tile(i,j)
                    tile.cycle_current()

            yee = getYee2D(node, conf)

            try:
                plot2dYee(axs[5], node, conf, 'jx')
                plot2dYee(axs[6], node, conf, 'jy')
                plot2dYee(axs[7], node, conf, 'jz')
                #saveVisz(lap, node, conf)
            except:
                pass

            for j in range(conf.Ny*conf.NyMesh):
                for i in range(conf.Nx*conf.NxMesh):
                    #print("({},{})".format(i,j))
                    self.assertAlmostEqual( yee_ref['jx'][i,j], yee['jx'][i,j], places=6 )
                    self.assertAlmostEqual( yee_ref['jy'][i,j], yee['jy'][i,j], places=6 )
                    self.assertAlmostEqual( yee_ref['jz'][i,j], yee['jz'][i,j], places=6 )



    def test_current_exchange(self):
        """test that current exchange is producing hand-build array"""


        #plt.fig = plt.figure(1, figsize=(5,7))
        #plt.rc('font', family='serif', size=12)
        #plt.rc('xtick')
        #plt.rc('ytick')
        #
        #gs = plt.GridSpec(5, 1)
        #
        #axs = []
        #for ai in range(5):
        #    axs.append( plt.subplot(gs[ai]) )


        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 10
        conf.NyMesh = 10
        conf.NzMesh = 1
        conf.ppc = 10
        conf.vel = 0.1
        conf.Nspecies = 2
        conf.me =-1.0        #electron mass-to-charge
        conf.mi = 1.0        #ion mass-to-charge

        conf.update_bbox()

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny, conf.Nz)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles(node, conf)
        #insert_em(node, conf, linear_field)
        insert_em(node, conf, const_field)
        inject(node, filler_xvel, conf) #injecting plasma particles

        #pusher   = pypic.BorisPusher()
        #fintp    = pypic.LinearInterpolator()
        currint  = pypic.ZigZag()
        analyzer = pypic.Analyzator()

        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                c = node.get_tile(i,j)
                yee = c.get_yee(0)
                for l in range(-3, conf.NxMesh+3):
                    for m in range(-3,conf.NyMesh+3):
                        for n in range(-3,conf.NzMesh+3):
                            yee.jx[l,m,n] = 1.0
                            yee.jy[l,m,n] = 1.0
                            yee.jz[l,m,n] = 1.0

        #deposit current
        #for j in range(node.get_Ny()):
        #    for i in range(node.get_Nx()):
        #for j in [1]:
        #    for i in [1]:
        #        tile = node.get_tile(i,j)
        #        currint.solve(tile)

        #exchange currents for the middle one only
        for j in [1]:
            for i in [1]:
                tile = node.get_tile(i,j)
                tile.exchange_currents(node)


        #plotNode(axs[0], node, conf)
        ##plot2dParticles(axs[0], node, conf, downsample=0.1)
        #plot2dYee(axs[1], node, conf, 'rho')
        #plot2dYee(axs[2], node, conf, 'jx')
        #plot2dYee(axs[3], node, conf, 'jy')
        #plot2dYee(axs[4], node, conf, 'jz')
        #saveVisz(-1, node, conf)

        # create reference 3 halo width array by hand
        ref = np.ones((conf.NxMesh, conf.NyMesh))
        ref[0:3,  : ]   = 2
        ref[-3:10,: ]   = 2
        ref[:, 0:3, ]   = 2
        ref[:, -3:10]   = 2

        ref[0:3,  0:3 ] = 4
        ref[7:10, 7:10] = 4
        ref[0:3,  7:10] = 4
        ref[7:10, 0:3 ] = 4
        #print(ref)

        for j in [1]:
            for i in [1]:
                c = node.get_tile(i,j)
                yee = c.get_yee(0)
                for l in range(conf.NxMesh):
                    for m in range(conf.NyMesh):
                        self.assertEqual(ref[l,m], yee.jx[l,m,0] )
                        self.assertEqual(ref[l,m], yee.jy[l,m,0] )
                        self.assertEqual(ref[l,m], yee.jz[l,m,0] )


    def test_current_deposit(self):
        """
        test that current deposit of + and - particles with same x and v
        produce zero current in total.
        """


        try:
            plt.fig = plt.figure(1, figsize=(5,7))
            plt.rc('font', family='serif', size=12)
            plt.rc('xtick')
            plt.rc('ytick')
            gs = plt.GridSpec(5, 1)

            axs = []
            for ai in range(5):
                axs.append( plt.subplot(gs[ai]) )
        except:
            pass


        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 10
        conf.NyMesh = 10
        conf.NzMesh = 1
        conf.ppc = 1
        conf.vel = 0.1
        conf.Nspecies = 2
        conf.me = -1.0        #electron mass-to-charge
        conf.mi =  1.0        #ion mass-to-charge

        conf.update_bbox()

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny, conf.Nz)
        node.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles(node, conf)
        #insert_em(node, conf, const_field)
        inject(node, filler_xvel, conf) #injecting plasma particles

        #pusher   = pypic.BorisPusher()
        #fintp    = pypic.LinearInterpolator()
        currint  = pypic.ZigZag()
        analyzer = pypic.Analyzator()


        #deposit current
        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                tile = node.get_tile(i,j)
                currint.solve(tile)

        #exchange currents for the middle one only
        #for j in [1]:
        #    for i in [1]:
        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                tile = node.get_tile(i,j)
                tile.exchange_currents(node)


        try:
            #plotNode(axs[0], node, conf)
            plot2dParticles(axs[0], node, conf, downsample=0.1)
            plot2dYee(axs[1], node, conf, 'rho')
            plot2dYee(axs[2], node, conf, 'jx')
            plot2dYee(axs[3], node, conf, 'jy')
            plot2dYee(axs[4], node, conf, 'jz')
            #saveVisz(-2, node, conf)
        except:
            pass

        for j in range(node.get_Ny()):
            for i in range(node.get_Nx()):
                c = node.get_tile(i,j)
                yee = c.get_yee(0)
                for l in range(conf.NxMesh):
                    for m in range(conf.NyMesh):
                        self.assertAlmostEqual(yee.jx[l,m,0], 0.0, places=7 )
                        self.assertAlmostEqual(yee.jy[l,m,0], 0.0, places=7 )
                        self.assertAlmostEqual(yee.jz[l,m,0], 0.0, places=7 )
                        #self.assertEqual(yee.jx[l,m,0], 0.0 )
                        #self.assertEqual(yee.jy[l,m,0], 0.0 )
                        #self.assertEqual(yee.jz[l,m,0], 0.0 )


