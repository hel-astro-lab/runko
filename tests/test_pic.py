from mpi4py import MPI
import unittest

import sys
import os
import numpy as np

import pytools  # runko python tools
import pycorgi
import pyrunko


#from visualize_pic import Particles
#from visualize_pic import plot2dParticles
#from visualize import plotNode
#from visualize import plot2dYee
#from visualize import pytools.visualize.get_grids_2D
#from visualize import saveVisz


try:
    import matplotlib.pyplot as plt
except:
    pass


#make tests deterministic by fixing the RNG seed
np.random.seed(0)


def density_profile(xloc, ispcs, conf):
    return conf.ppc


def filler_no_velocity(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    #zz = xloc[2] + np.random.rand()
    zz = 0.0

    x0 = [xx, yy, zz]
    u0 = [0.0, 0.0, 0.0]
    return x0, u0

def filler_no_velocity_3d(xloc, ispcs, conf):
    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    zz = xloc[2] + np.random.rand()

    x0 = [xx, yy, zz]
    u0 = [0.0, 0.0, 0.0]
    return x0, u0


# random number between [xmin, xmax]
def randab(xmin, xmax):
    return xmin + (xmax-xmin)*np.random.rand()


def filler(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    #zz = xloc[2] + np.random.rand()
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
    uz = 0.0

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]

    return x0, u0



def filler3D(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    zz = xloc[2] + np.random.rand()

    ur = conf.vel
    uc = randab(0.0, 2.0*np.pi) 
    us = randab(0.0, 1.0*np.pi) 

    #3D
    ux = ur*np.sin( uc )*np.sin(us)
    uy = ur*np.cos( uc )*np.sin(us)
    uz = ur*np.cos(us)

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]

    return x0, u0

def filler_xvel(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    #zz = xloc[2] + np.random.rand()

    #xx = xloc[0] #+ 0.5
    #yy = xloc[1] #+ 0.5
    zz = 0.0

    x0 = [xx, yy, zz]

    #if (conf.NyMesh < yy < 2.0*conf.NyMesh) and (conf.NxMesh < xx < 2.0*conf.NxMesh):
    #    u0 = [0.0, 1.0, 0.0]
    #else:
    #    u0 = [0.0, 0.0, 0.0]

    u0 = [+1.0, +1.0, +1.0]
    #u0 = [-1.0, -1.0, -1.0]

    return x0, u0


def filler_xvel3D(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    zz = xloc[2] + np.random.rand()
    x0 = [xx, yy, zz]
    u0 = [+1.0, +1.0, +1.0]
    return x0, u0


def zero_field(x,y,z):
    return 0.0

def const_field(x, y, z):
    return 1.0

def linear_field(x, y, z):
    #print("x = {} y = {} z = {}".format(x,y,z))
    return 1.1*x + 1.2*y

def linear_field_3d(x, y, z):
    return 1.0*x + 2.0*y + 3.0*z

def qudratic_field_3d(x, y, z):
    return 1.0*x*x + 2.0*y*y + 3.0*z*z + 4.*x + 5*y + 6*z + 7


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(grid, conf, ffunc, zero_field=False):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    #for k in range(grid.get_Nz()):
    #    for j in range(grid.get_Ny()):
    #        for i in range(grid.get_Nx()):

    for cid in grid.get_local_tiles():
        #c = grid.get_tile(i,j,k)
        c = grid.get_tile(cid)
        gs = c.get_grids()
        i,j,k = pytools.get_index(c, conf)

        for n in range(conf.NzMesh):
            for m in range(conf.NyMesh):
                for l in range(conf.NxMesh):
                    # get x_i,j,k
                    xloc0 = pytools.ind2loc((i,j,k), (l,m,n), conf)

                    #get x_i+1/2, x_j+1/2, x_k+1/2
                    #xloc1 = pytools.ind2loc((i,j,k), (l+1,m,  n),   conf)
                    #yloc1 = pytools.ind2loc((i,j,k), (l,  m+1,n),   conf)
                    #zloc1 = pytools.ind2loc((i,j,k), (l,  m,  n+1), conf)

                    # values in Yee lattice corners
                    #xcor = xloc0[0]
                    #ycor = xloc0[1]
                    #zcor = xloc0[2]

                    # values in Yee lattice mids
                    #xmid = 0.5*(xloc0[0] + xloc1[0])
                    #ymid = 0.5*(xloc0[1] + yloc1[1])
                    #zmid = 0.5*(xloc0[2] + zloc1[2])

                    #val = ffunc(xmid, ymid, zmid)

                    # enforce Yee lattice structure
                    #gs.ex[l,m,n] = ffunc(xmid, ycor, zcor)
                    #gs.ey[l,m,n] = ffunc(xcor, ymid, zcor)
                    #gs.ez[l,m,n] = ffunc(xcor, ycor, zmid)

                    #gs.bx[l,m,n] = ffunc(xcor, ymid, zmid)
                    #gs.by[l,m,n] = ffunc(xmid, ycor, zmid)
                    #gs.bz[l,m,n] = ffunc(xmid, ymid, zcor)

                    #gs.jx[l,m,n] = ffunc(xmid, ymid, zmid)
                    #gs.jy[l,m,n] = ffunc(xmid, ymid, zmid)
                    #gs.jz[l,m,n] = ffunc(xmid, ymid, zmid)

                    #TODO
                    xc = xloc0[0]
                    yc = xloc0[1]
                    zc = xloc0[2]
                    #gs.ex[l,m,n] = ffunc(xc,yc-0.5,zc-0.5)
                    #gs.ey[l,m,n] = ffunc(xc-0.5,yc,zc-0.5)
                    #gs.ez[l,m,n] = ffunc(xc-0.5,yc-0.5,zc)
                    #gs.bx[l,m,n] = ffunc(xc-0.5,yc,zc)
                    #gs.by[l,m,n] = ffunc(xc,yc-0.5,zc)
                    #gs.bz[l,m,n] = ffunc(xc,yc,zc-0.5)

                    gs.ex[l,m,n] = ffunc(xc+0.5,yc,    zc)
                    gs.ey[l,m,n] = ffunc(xc,    yc+0.5,zc)
                    gs.ez[l,m,n] = ffunc(xc,    yc,    zc+0.5)
                    gs.bx[l,m,n] = ffunc(xc,    yc+0.5,zc+0.5)
                    gs.by[l,m,n] = ffunc(xc+0.5,yc,    zc+0.5)
                    gs.bz[l,m,n] = ffunc(xc+0.5,yc+0.5,zc)

                    gs.jx[l,m,n] = ffunc(xc,yc,zc)
                    gs.jy[l,m,n] = ffunc(xc,yc,zc)
                    gs.jz[l,m,n] = ffunc(xc,yc,zc)

                    if not(zero_field):
                        gs.ex[l,m,n] += +0.0
                        gs.ey[l,m,n] += +1.0
                        gs.ez[l,m,n] += +2.0 
                        gs.bx[l,m,n] += +3.0
                        gs.by[l,m,n] += +4.0
                        gs.bz[l,m,n] += +5.0





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

    qe = 1.0
    qi =-1.0

    oneD   = False
    twoD   = False
    threeD = False

    prtcl_types = ['e-', 'e+']

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

        conf.twoD = True
        conf.Nx = 3
        conf.Ny = 3
        conf.Ny = 1
        conf.update_bbox()

        conf.vel = 0.3

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        pytools.pic.load_tiles(grid, conf)
        insert_em(grid, conf, zero_field, zero_field=True)
        pytools.pic.inject(grid, filler, density_profile, conf) #pytools.pic.injecting plasma particles

        # push particles couple of times to make them leak into neighboring tiles
        pusher = pyrunko.pic.twoD.BorisPusher()
        fintp  = pyrunko.pic.twoD.LinearInterpolator()


        for lap in range(20):

            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                fintp.solve(tile)

            #plot2dParticles(axs[0], grid, conf)
            #saveVisz(lap, grid, conf)

            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                pusher.solve(tile)

            ##################################################
            # communication

            #update particle boundaries
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                tile.check_outgoing_particles()

            # global mpi exchange (independent)
            #for cid in grid.get_boundary_tiles():
            #    tile = grid.get_tile(cid)
            #    tile.pack_outgoing_particles()

            ## MPI global exchange
            ## transfer primary and extra data
            #grid.send_data(0) #(indepdendent)
            #grid.send_data(1) #(indepdendent)

            #grid.recv_data(0) #(indepdendent)
            #grid.recv_data(1) #(indepdendent)

            #grid.wait_data(0) #(indepdendent)
            #grid.wait_data(1) #(indepdendent)

            # global unpacking (independent)
            #for cid in grid.get_virtual_tiles(): 
            #    tile = grid.get_tile(cid)
            #    tile.unpack_incoming_particles()
            #    tile.check_outgoing_particles()

            # transfer local + global
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                tile.get_incoming_particles(grid)

            # delete local transferred particles
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                tile.delete_transferred_particles()

            # delete global mpi particles
            #for cid in grid.get_virtual_tiles(): 
            #    tile = grid.get_tile(cid)
            #    tile.delete_all_particles()

        # count how many particles we now have
        n_particles = 0
        #for i in range(conf.Nx):
        #    for j in range(conf.Ny):
        #        for k in range(conf.Nz):
        #    cid = grid.id(i,j)
        for cid in grid.get_local_tiles():
            c = grid.get_tile(cid)
            i,j,k = pytools.get_index(c, conf)

            container = c.get_container(0)
            #print("({},{},{}) has {}".format(i,j,k,len(container.loc(0))))
            n_particles += len(container.loc(0))

            #self.assertTrue( 0.0 <= container.loc(0) <= conf.xmax )
            #self.assertTrue( 0.0 <= container.loc(1) <= conf.ymax )
            #self.assertTrue( 0.0 <= container.loc(2) <= conf.zmax )

            for prtcl in range(len(container.loc(0))):
                #print("{} {} {} maxs {} {} {} id {}/{}".format( 
                #container.loc(0)[prtcl], 
                #container.loc(1)[prtcl], 
                #container.loc(2)[prtcl], 
                #conf.xmax, conf.ymax, conf.zmax, 
                #container.id(0)[prtcl], 
                #container.id(1)[prtcl], 
                #))

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
                self.assertAlmostEqual( vel, conf.vel, places=5 )

        tot_particles = (conf.Nx*conf.NxMesh *
                        conf.Ny*conf.NyMesh *
                        conf.Nz*conf.NzMesh *
                        conf.ppc)

        # assert that there is equal number of particles as we began with
        self.assertEqual( tot_particles, n_particles )


    def test_communication3D(self):

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

        sch  = pytools.Scheduler(print_banner=False) # Set MPI aware task scheduler
        sch.verbose=0

        timer = pytools.Timer()
        timer.verbose = 0  # 0 normal; 1 - debug mode
        sch.timer = timer # remember to update scheduler

        conf = Conf()
        conf.ppc = 2
        conf.threeD = True
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 3
        conf.Nx = 4
        conf.Ny = 4
        conf.Nz = 4
        conf.update_bbox()
        conf.mpi_task_mode = False
        conf.vel = 0.3

        sch.conf = conf # remember to update conf to scheduler


        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
        sch.grid = grid # remember to update scheduler

        pytools.balance_mpi(grid, conf, do_print=False) #balance MPI ranks; Hilbert curve optimization 

        pytools.pic.load_tiles(grid, conf)
        insert_em(grid, conf, zero_field, zero_field=True)
        pytools.pic.inject(grid, filler3D, density_profile, conf) #pytools.pic.injecting plasma particles

        # push particles couple of times to make them leak into neighboring tiles
        sch.pusher  = pyrunko.pic.threeD.BorisPusher()
        sch.fintp   = pyrunko.pic.threeD.LinearInterpolator()
        sch.currint = pyrunko.pic.threeD.ZigZag()

        # start MPI
        # update boundaries
        grid.analyze_boundaries()
        grid.send_tiles()
        grid.recv_tiles()
        MPI.COMM_WORLD.barrier()
        rank = MPI.COMM_WORLD.Get_rank()

        # load virtual mpi halo tiles
        pytools.pic.load_virtual_tiles(grid, conf)

        for lap in range(10):

            # interpolate fields and push particles in x and u
            sch.operate( dict(name='interp_em', solver='fintp',  method='solve', nhood='local', ) )
            sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', ) )

            # --------------------------------------------------
            # particle communication (only local/boundary tiles)

            # local and global particle exchange 
            sch.operate( dict(name='check_outg_prtcls',     solver='tile',  method='check_outgoing_particles',     nhood='local', ) )
            sch.operate( dict(name='pack_outg_prtcls',      solver='tile',  method='pack_outgoing_particles',      nhood='boundary', ) )

            sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p1',                           nhood='all', ) )
            sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p2',                           nhood='all', ) )

            sch.operate( dict(name='unpack_vir_prtcls',     solver='tile',  method='unpack_incoming_particles',    nhood='virtual', ) )
            sch.operate( dict(name='check_outg_vir_prtcls', solver='tile',  method='check_outgoing_particles',     nhood='virtual', ) )
            sch.operate( dict(name='get_inc_prtcls',        solver='tile',  method='get_incoming_particles',       nhood='local', args=[grid,]) )

            sch.operate( dict(name='del_trnsfrd_prtcls',    solver='tile',  method='delete_transferred_particles', nhood='local', ) )
            sch.operate( dict(name='del_vir_prtcls',        solver='tile',  method='delete_all_particles',         nhood='virtual', ) )


            # count how many particles we now have
            n_particles = 0
            #for cid in grid.get_local_tiles():
            #for i in range(conf.Nx):
            #    for j in range(conf.Ny):
            #        for k in range(conf.Nz):
            for cid in grid.get_local_tiles():
                c = grid.get_tile(cid)
                i,j,k = pytools.get_index(c, conf)

                container = c.get_container(0)
                #print("({},{},{}) has {}".format(i,j,k,len(container.loc(0))))
                n_particles += len(container.loc(0))
                #n_particles += len(container.size() )

                #print()
                #print("num prtcls", i,j,k, "len(loc())", len(container.loc(0)), " cont.size:", container.size())
                self.assertTrue( len(container.loc(0)) == container.size() )

                #self.assertTrue( 0.0 <= container.loc(0) <= conf.xmax )
                #self.assertTrue( 0.0 <= container.loc(1) <= conf.ymax )
                #self.assertTrue( 0.0 <= container.loc(2) <= conf.zmax )

                for prtcl in range(len(container.loc(0))):

                    if False: #rank == 0:
                        print("r: {} ----prtcl {} x={} y={} z={} vx={} vy={} vz={} w={} | tile lims {} {} {}".format(
                        rank,
                        prtcl, 
                        container.loc(0)[prtcl],
                        container.loc(1)[prtcl],
                        container.loc(2)[prtcl],
                        container.vel(0)[prtcl],
                        container.vel(1)[prtcl],
                        container.vel(2)[prtcl],
                        container.wgt()[ prtcl],
                        conf.xmax,
                        conf.ymax,
                        conf.zmax),)

                    # check location
                    self.assertTrue( 0.0 <= container.loc(0)[prtcl] <= conf.xmax )
                    self.assertTrue( 0.0 <= container.loc(1)[prtcl] <= conf.ymax )
                    self.assertTrue( 0.0 <= container.loc(2)[prtcl] <= conf.zmax )

                    # check velocity 
                    velx = container.vel(0)[prtcl]
                    vely = container.vel(1)[prtcl]
                    velz = container.vel(2)[prtcl]
                    vel = np.sqrt( velx*velx + vely*vely + velz*velz )
                    self.assertAlmostEqual( vel, conf.vel, places=5 )

            tot_particles = (conf.Nx*conf.NxMesh *
                             conf.Ny*conf.NyMesh *
                             conf.Nz*conf.NzMesh *
                            conf.ppc)

            #print('r: {} np: {}'.format(rank, n_particles))
            n_particles = MPI.COMM_WORLD.reduce( n_particles, op=MPI.SUM, root=0)

            # assert that there is equal number of particles as we began with
            if rank == 0:
                self.assertEqual( tot_particles, n_particles)



    def test_const_field_interpolation_2d(self):

        conf = Conf()
        conf.twoD = True
        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        pytools.pic.load_tiles(grid, conf)
        insert_em(grid, conf, const_field)
        pytools.pic.inject(grid, filler_no_velocity, density_profile, conf) #pytools.pic.injecting plasma particles

        ##update boundaries
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            tile.update_boundaries(grid)

        #interpolate emf
        fintps = []

        fintps.append( pyrunko.pic.twoD.LinearInterpolator() )
        fintps.append( pyrunko.pic.twoD.QuadraticInterpolator() )

        for fintp in fintps:
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                fintp.solve(tile)

            #test results
            for cid in grid.get_local_tiles():
                c = grid.get_tile(cid)
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
                    self.assertAlmostEqual(container.ex(i), ex_ref, places=5)
                    self.assertAlmostEqual(container.ey(i), ey_ref, places=5)
                    self.assertAlmostEqual(container.ez(i), ez_ref, places=5)

                    bx_ref = 4.0
                    by_ref = 5.0
                    bz_ref = 6.0
                    self.assertAlmostEqual(container.bx(i), bx_ref, places=5)
                    self.assertAlmostEqual(container.by(i), by_ref, places=5)
                    self.assertAlmostEqual(container.bz(i), bz_ref, places=5)


    def test_const_field_interpolation_3d(self):
        conf = Conf()

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3

        conf.NxMesh = 6
        conf.NyMesh = 6
        conf.NzMesh = 6

        conf.threeD = True
        conf.update_bbox()

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)


        pytools.pic.load_tiles(grid, conf)
        insert_em(grid, conf, const_field)
        pytools.pic.inject(grid, filler_no_velocity_3d, density_profile, conf) #pytools.pic.injecting plasma particles

        ##update boundaries
        for tile in pytools.tiles_local(grid):
            tile.update_boundaries(grid)

        #interpolate emf
        fintps = []
        fintps.append( pyrunko.pic.threeD.LinearInterpolator() )
        fintps.append( pyrunko.pic.threeD.QuadraticInterpolator() )
        #fintps.append( pyrunko.pic.threeD.CubicInterpolator() )
        fintps.append( pyrunko.pic.threeD.QuarticInterpolator() )

        for fintp in fintps:
            #interpolate
            for tile in pytools.tiles_local(grid):
                fintp.solve(tile)

            #test results
            for tile in pytools.tiles_local(grid):
                container = tile.get_container(0)

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
                    self.assertAlmostEqual(container.ex(i), ex_ref, places=5)
                    self.assertAlmostEqual(container.ey(i), ey_ref, places=5)
                    self.assertAlmostEqual(container.ez(i), ez_ref, places=5)

                    bx_ref = 4.0
                    by_ref = 5.0
                    bz_ref = 6.0
                    self.assertAlmostEqual(container.bx(i), bx_ref, places=5)
                    self.assertAlmostEqual(container.by(i), by_ref, places=5)
                    self.assertAlmostEqual(container.bz(i), bz_ref, places=5)


    def test_linear_field_interpolation_2d(self):

        conf = Conf()
        conf.twoD = True
        conf.Nx = 3
        conf.Ny = 3

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        pytools.pic.load_tiles(grid, conf)
        insert_em(grid, conf, linear_field)
        pytools.pic.inject(grid, filler_no_velocity, density_profile, conf) #pytools.pic.injecting plasma particles

        ##update boundaries
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.update_boundaries(grid)

        #interpolate emf
        fintps = []

        fintps.append( pyrunko.pic.twoD.LinearInterpolator() )
        fintps.append( pyrunko.pic.twoD.QuadraticInterpolator() )

        for fintp in fintps:
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    fintp.solve(tile)


            #test results; avoid boundaries because they are cyclic
            for i in range(1,conf.Nx-1):
                for j in range(1,conf.Ny-1):

                    cid = grid.id(i,j)
                    c = grid.get_tile(cid)
                    container = c.get_container(0)

                    xx = container.loc(0)
                    yy = container.loc(1)
                    zz = container.loc(2)

                    for i, x in enumerate(xx):
                        #print(i)

                        ref = linear_field(xx[i], yy[i], zz[i]) #exact/true solution

                        #print("asserting {} {} {} vs {}".format(xx[i], yy[i], zz[i], ref))
                        self.assertAlmostEqual(container.ex(i), ref+0.0, places=4)
                        self.assertAlmostEqual(container.ey(i), ref+1.0, places=4)
                        self.assertAlmostEqual(container.ez(i), ref+2.0, places=4)
                        self.assertAlmostEqual(container.bx(i), ref+3.0, places=4)
                        self.assertAlmostEqual(container.by(i), ref+4.0, places=4)
                        self.assertAlmostEqual(container.bz(i), ref+5.0, places=4)


    def test_linear_field_interpolation_3d(self):
        conf = Conf()

        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3

        conf.NxMesh = 6
        conf.NyMesh = 6
        conf.NzMesh = 6

        conf.threeD = True
        conf.update_bbox()

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, 
                           conf.ymin, conf.ymax, 
                           conf.zmin, conf.zmax)
        pytools.pic.load_tiles( grid, conf)
        insert_em( grid, conf, linear_field_3d)
        pytools.pic.inject(grid, filler3D, density_profile, conf) #pytools.pic.injecting plasma particles

        ##update boundaries
        for tile in pytools.tiles_local(grid):
            tile.update_boundaries(grid, iarr=[0,1,2])

        #interpolate emf
        fintps = []
        fintps.append( pyrunko.pic.threeD.LinearInterpolator() )
        fintps.append( pyrunko.pic.threeD.QuadraticInterpolator() )
        #fintps.append( pyrunko.pic.threeD.CubicInterpolator() )
        fintps.append( pyrunko.pic.threeD.QuarticInterpolator() )

        for fintp in fintps:
            #interpolate
            for tile in pytools.tiles_local(grid):
                fintp.solve(tile)

            #test results
            #for tile in pytools.tiles_local(grid):
            #avoid boundaries because the test function does not wrap nicely
            for i in range(1,conf.Nx-1):
                for j in range(1,conf.Ny-1):
                    for k in range(1,conf.Nz-1):
                        cid = grid.id(i,j,k)
                        tile = grid.get_tile(cid)

                        container = tile.get_container(0)
                        xx = container.loc(0)
                        yy = container.loc(1)
                        zz = container.loc(2)
                        ux = container.vel(0)
                        uy = container.vel(1)
                        uz = container.vel(2)

                        for i, _ in enumerate(xx):
                            ref = linear_field_3d(xx[i], yy[i], zz[i]) #exact/true solution

                            #print("asserting {} {} {} vs {}".format(xx[i], yy[i], zz[i], ref))
                            self.assertAlmostEqual( container.ex(i), ref+0.0, places=4 )
                            self.assertAlmostEqual( container.ey(i), ref+1.0, places=4 )
                            self.assertAlmostEqual( container.ez(i), ref+2.0, places=4 )
                                                                                      
                            self.assertAlmostEqual( container.bx(i), ref+3.0, places=4 )
                            self.assertAlmostEqual( container.by(i), ref+4.0, places=4 )
                            self.assertAlmostEqual( container.bz(i), ref+5.0, places=4 )


    def skip_test_filters(self):
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
        conf.twoD = True
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 1
        conf.ppc = 10
        conf.vel = 0.1
        conf.update_bbox()

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        pytools.pic.load_tiles(grid, conf)
        #insert_em(grid, conf, linear_field)
        insert_em(grid, conf, zero_field)
        pytools.pic.inject(grid, filler, density_profile, conf) #pytools.pic.injecting plasma particles
        #pytools.pic.inject(grid, filler_xvel, conf) #pytools.pic.injecting plasma particles

        #pusher   = pyrunko.pic.twoD.BorisPusher()
        #fintp    = pyrunko.pic.twoD.LinearInterpolator()
        currint  = pyrunko.pic.twoD.ZigZag()

        flt =  pytools.Filter(conf.NxMesh, conf.NyMesh)
        flt.init_gaussian_kernel(1.0, 1.0)


        #for lap in range(0, conf.Nt):
        for lap in range(1):

            #update boundaries
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    tile.update_boundaries(grid)

            #deposit current
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    currint.solve(tile)

            #exchange currents
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    tile.exchange_currents(grid)

            #try:
            #    #plotNode(axs[0], grid, conf)
            #    #plot2dParticles(axs[0], grid, conf, downsample=0.1)
            #    #plot2dYee(axs[1], grid, conf, 'rho')
            #    #plot2dYee(axs[2], grid, conf, 'jx')
            #    #plot2dYee(axs[3], grid, conf, 'jy')
            #    #plot2dYee(axs[4], grid, conf, 'jz')
            #except:
            #    pass

            grids_refs = pytools.visualize.get_grids_2D(grid, conf)

            #filter
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    #print(" i j ({},{})".format(i,j))
                    tile = grid.get_tile(i,j)
                    flt.get_padded_current(tile, grid)

                    # fourier space filtering
                    flt.fft_image_forward()
                    flt.apply_kernel()
                    flt.fft_image_backward()

                    # direct filtering
                    #for fj in range(1):
                    #    flt.direct_convolve_3point()
                    flt.set_current(tile)

            #cycle new and temporary currents
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    tile.cycle_current()

            gs = pytools.visualize.get_grids_2D(grid, conf)

            #try:
            #    plot2dYee(axs[5], grid, conf, 'jx')
            #    plot2dYee(axs[6], grid, conf, 'jy')
            #    plot2dYee(axs[7], grid, conf, 'jz')
            #    #saveVisz(lap, grid, conf)
            #except:
            #    pass

            for j in range(conf.Ny*conf.NyMesh):
                for i in range(conf.Nx*conf.NxMesh):
                    #print("({},{})".format(i,j))
                    self.assertAlmostEqual( grids_refs['jx'][i,j], gs['jx'][i,j], places=5 )
                    self.assertAlmostEqual( grids_refs['jy'][i,j], gs['jy'][i,j], places=5 )
                    self.assertAlmostEqual( grids_refs['jz'][i,j], gs['jz'][i,j], places=5 )



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
        conf.twoD = True
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

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        pytools.pic.load_tiles(grid, conf)
        #insert_em(grid, conf, linear_field)
        insert_em(grid, conf, const_field)
        pytools.pic.inject(grid, filler_xvel, density_profile, conf) #pytools.pic.injecting plasma particles

        #pusher   = pyrunko.pic.twoD.BorisPusher()
        #fintp    = pyrunko.pic.twoD.LinearInterpolator()
        currint  = pyrunko.pic.twoD.ZigZag()

        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                c = grid.get_tile(i,j)
                gs = c.get_grids(0)
                for l in range(-3, conf.NxMesh+3):
                    for m in range(-3,conf.NyMesh+3):
                        for n in range(-3,conf.NzMesh+3):
                            gs.jx[l,m,n] = 1.0
                            gs.jy[l,m,n] = 1.0
                            gs.jz[l,m,n] = 1.0

        #deposit current
        #for j in range(grid.get_Ny()):
        #    for i in range(grid.get_Nx()):
        #for j in [1]:
        #    for i in [1]:
        #        tile = grid.get_tile(i,j)
        #        currint.solve(tile)

        #exchange currents for the middle one only
        for j in [1]:
            for i in [1]:
                tile = grid.get_tile(i,j)
                tile.exchange_currents(grid)


        #plotNode(axs[0], grid, conf)
        ##plot2dParticles(axs[0], grid, conf, downsample=0.1)
        #plot2dYee(axs[1], grid, conf, 'rho')
        #plot2dYee(axs[2], grid, conf, 'jx')
        #plot2dYee(axs[3], grid, conf, 'jy')
        #plot2dYee(axs[4], grid, conf, 'jz')
        #saveVisz(-1, grid, conf)

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
                c = grid.get_tile(i,j)
                gs = c.get_grids(0)
                for l in range(conf.NxMesh):
                    for m in range(conf.NyMesh):
                        self.assertEqual(ref[l,m], gs.jx[l,m,0] )
                        self.assertEqual(ref[l,m], gs.jy[l,m,0] )
                        self.assertEqual(ref[l,m], gs.jz[l,m,0] )


    def test_current_deposit(self):
        #test that current deposit of + and - particles with same x and v
        #produce zero current in total.


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
        conf.twoD = True
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

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        pytools.pic.load_tiles(grid, conf)
        #insert_em(grid, conf, const_field)
        pytools.pic.inject(grid, filler_xvel, density_profile, conf) #pytools.pic.injecting plasma particles

        #pusher   = pyrunko.pic.twoD.BorisPusher()
        #fintp    = pyrunko.pic.twoD.LinearInterpolator()
        currints = []
        currints.append( pyrunko.pic.twoD.ZigZag() )
        currints.append( pyrunko.pic.twoD.ZigZag_2nd() )
        currints.append( pyrunko.pic.twoD.ZigZag_3rd() )
        currints.append( pyrunko.pic.twoD.ZigZag_4th() )
        
        for currint in currints:

            #deposit current
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    currint.solve(tile)

            #exchange currents for all tiles
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    tile.exchange_currents(grid)

            #try:
            #    #plotNode(axs[0], grid, conf)
            #    #plot2dParticles(axs[0], grid, conf, downsample=0.1)
            #    plot2dYee(axs[1], grid, conf, 'rho')
            #    plot2dYee(axs[2], grid, conf, 'jx')
            #    plot2dYee(axs[3], grid, conf, 'jy')
            #    plot2dYee(axs[4], grid, conf, 'jz')
            #    #saveVisz(-2, grid, conf)
            #except:
            #    pass

            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    c = grid.get_tile(i,j)
                    gs = c.get_grids(0)
                    for l in range(conf.NxMesh):
                        for m in range(conf.NyMesh):
                            self.assertAlmostEqual(gs.jx[l,m,0], 0.0, places=5 )
                            self.assertAlmostEqual(gs.jy[l,m,0], 0.0, places=5 )
                            self.assertAlmostEqual(gs.jz[l,m,0], 0.0, places=5 )
                            #self.assertEqual(gs.jx[l,m,0], 0.0 )
                            #self.assertEqual(gs.jy[l,m,0], 0.0 )
                            #self.assertEqual(gs.jz[l,m,0], 0.0 )

    def test_current_deposit3D(self):
        #test that current deposit of + and - particles with same x and v
        #produce zero current in total.

        conf = Conf()
        conf.threeD = True
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 3
        conf.NxMesh = 6
        conf.NyMesh = 7
        conf.NzMesh = 8

        conf.ppc = 1
        conf.vel = 0.1
        conf.Nspecies = 2
        conf.me = -1.0        #electron mass-to-charge
        conf.mi =  1.0        #ion mass-to-charge

        conf.update_bbox()

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
        pytools.pic.load_tiles(grid, conf)
        pytools.pic.inject(grid, filler_xvel3D, density_profile, conf) #pytools.pic.injecting plasma particles

        currints = []
        currints.append( pyrunko.pic.threeD.ZigZag() )
        currints.append( pyrunko.pic.threeD.ZigZag_2nd() )
        currints.append( pyrunko.pic.threeD.ZigZag_3rd() )
        currints.append( pyrunko.pic.threeD.ZigZag_4th() )

        # loop over different depositers
        for currint in currints:
            #print(currint)

            #deposit current
            #print('deposit')
            for k in range(grid.get_Nz()):
                for j in range(grid.get_Ny()):
                    for i in range(grid.get_Nx()):
                        tile = grid.get_tile(i,j,k)
                        currint.solve(tile)

            #print('exchange')
            #exchange currents for all tiles
            #for k in range(grid.get_Nz()):
            #    for j in range(grid.get_Ny()):
            #        for i in range(grid.get_Nx()):
            #            tile = grid.get_tile(i,j,k)
            #            tile.exchange_currents(grid)

            #print('assert')
            for k in range(grid.get_Nz()):
                for j in range(grid.get_Ny()):
                    for i in range(grid.get_Nx()):
                        c = grid.get_tile(i,j,k)
                        gs = c.get_grids(0)
                        for l in range(conf.NxMesh):
                            for m in range(conf.NyMesh):
                                for n in range(conf.NzMesh):
                                    self.assertAlmostEqual(gs.jx[l,m,n], 0.0, places=5 )
                                    self.assertAlmostEqual(gs.jy[l,m,n], 0.0, places=5 )
                                    self.assertAlmostEqual(gs.jz[l,m,n], 0.0, places=5 )


    def test_test_particle_initialization(self):

        conf = Conf()
        conf.twoD = True

        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.Nx = 3
        conf.Ny = 3
        conf.Ny = 1
        conf.ppc = 1
        conf.update_bbox()

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        # this calls internally test particle addition
        pytools.pic.load_tiles(grid, conf)



    def test_current_array_aligning(self):
        #test current deposition + interpolation to see array aligning in practise
        #this is also closest test we have for a full integration test (non-mpi)

        do_plots = False
        try:
            if do_plots:
                plt.fig = plt.figure(1, figsize=(8,10))
                plt.rc('font', family='serif', size=8)
                plt.rc('xtick')
                plt.rc('ytick')
                
                gs = plt.GridSpec(4, 3)
                gs.update(hspace = 0.0)
                
                axs = []
                for ai in range(12):
                    axs.append( plt.subplot(gs[ai]) )
        except:
            #print()
            pass

        conf = Conf()
        conf.twoD = True

        conf.outdir = "align_test"
        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.ppc = 1
        conf.npasses = 0
        conf.Nt = 2
        conf.update_bbox()

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        pytools.pic.load_tiles(grid, conf)

        #pytools.pic.inject(grid, filler_xvel, conf) #pytools.pic.injecting plasma particles
        #pytools.pic.inject by hand
        #-------------------------------------------------- 
        cid    = grid.id(0,0)
        c      = grid.get_tile(cid) #get cell ptr

        container = c.get_container(0) #ispcs
        container.set_keygen_state(0, 0) #number, rank

        #pytools.pic.inject random tricky points 
        
        #x moving
        x0 = [0.0, 0.0, 0.5]
        u0 = [+0.9, 0.0, 0.0]
        container.add_particle(x0, u0, 1.0)

        x0 = [0.0, 2.99, 0.5]
        u0 = [-0.9, 0.0, 0.0]
        container.add_particle(x0, u0, 1.0)

        #y moving
        x0 = [2.99, 2.99, 0.5]
        u0 = [0.0,-0.9, 0.0]
        container.add_particle(x0, u0, 1.0)

        x0 = [0.0, 0.00, 0.5]
        u0 = [0.0,+0.9, 0.0]
        container.add_particle(x0, u0, 1.0)

        #z moving
        x0 = [2.99, 2.99, 0.5]
        u0 = [0.0, 0.0,-0.9]
        container.add_particle(x0, u0, 1.0)

        x0 = [2.0, 0.0, 0.5]
        u0 = [0.0, 0.0,+0.9]
        container.add_particle(x0, u0, 1.0)

        #x0 = [0.5, 0.5, 0.5]
        #for u0 in [
        #        [ 0.9, 0.0, 0.0],
        #        [-0.9, 0.0, 0.0],
        #        [ 0.0, 0.9, 0.0],
        #        [ 0.0,-0.9, 0.0]]:
        #    container.add_particle(x0, u0, 1.0)

        #-------------------------------------------------- 

        pusher   = pyrunko.pic.twoD.BorisPusher()
        fldprop  = pyrunko.emf.twoD.FDTD2()
        fintp    = pyrunko.pic.twoD.LinearInterpolator()
        currint  = pyrunko.pic.twoD.ZigZag()
        flt      = pyrunko.emf.twoD.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        lap = 0
        for lap in range(lap, conf.Nt):

            #update tile boundaries
            for cid in grid.get_tile_ids():
                tile = grid.get_tile(cid)
                tile.update_boundaries(grid)
    
            #interpolate emf
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                fintp.solve(tile)
    
            #push
            #for cid in grid.get_local_tiles():
            #    tile = grid.get_tile(cid)
            #    pusher.solve(tile)
    
            #push half B
            for cid in grid.get_tile_ids():
                tile = grid.get_tile(cid)
                fldprop.push_half_b(tile)
    
            for cid in grid.get_tile_ids():
                tile = grid.get_tile(cid)
                tile.update_boundaries(grid)
    
            #push full E
            for cid in grid.get_tile_ids():
                tile = grid.get_tile(cid)
                fldprop.push_e(tile)
    
            #current
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                currint.solve(tile)
    
            for cid in grid.get_tile_ids():
                tile = grid.get_tile(cid)
                tile.exchange_currents(grid)
    
            #particle movement
            #for cid in grid.get_local_tiles():
            #    tile = grid.get_tile(cid)
            #    tile.check_outgoing_particles()
    
            #for cid in grid.get_local_tiles():
            #    tile = grid.get_tile(cid)
            #    tile.get_incoming_particles(grid)
    
            #for cid in grid.get_local_tiles():
            #    tile = grid.get_tile(cid)
            #    tile.delete_transferred_particles()
    
            #filtering
            #for fj in range(conf.npasses):
    
            #    #filter each tile
            #    for cid in grid.get_local_tiles():
            #        tile = grid.get_tile(cid)
            #        flt.solve(tile)
    
            #    #get halo boundaries
            #    for cid in grid.get_local_tiles():
            #        tile = grid.get_tile(cid)
            #        tile.update_boundaries(grid)
    
            for cid in grid.get_tile_ids():
                tile = grid.get_tile(cid)
                tile.deposit_current()
    
            #--------------------------------------------------
            # end of cycle // analyzing next
    
            # TODO can this be removed?
            #build analysis statistics
            #for cid in grid.get_local_tiles():
            #    tile = grid.get_tile(cid)


            xx = container.loc(0)
            if do_plots:
                for i in range(len(xx)):
                    print("i = ", i)
                    print("container ex:", container.ex(i))
                    print("container ey:", container.ey(i))
                    print("container ez:", container.ez(i))
                    print("container bx:", container.bx(i))
                    print("container by:", container.by(i))
                    print("container bz:", container.bz(i))
    
            #plot
            gs = pytools.visualize.get_grids_2D(grid, conf)


            if do_plots:
                print("jx:", gs['jx'])
                print("jy:", gs['jy'])
                print("jz:", gs['jz'])
                print("sum of current jx:", np.sum(gs['jx']))
                print("sum of current jy:", np.sum(gs['jy']))
                print("sum of current jz:", np.sum(gs['jz']))
                #plotNode(axs[0], grid, conf)
                #plot2dYee(axs[1],  gs, grid, conf, 'rho', label_title=True)

                #plot2dYee(axs[3],  gs, grid, conf, 'jx' , label_title=True)
                #plot2dYee(axs[4],  gs, grid, conf, 'jy' , label_title=True)
                #plot2dYee(axs[5],  gs, grid, conf, 'jz' , label_title=True)
                #plot2dYee(axs[6],  gs, grid, conf, 'ex' , label_title=True)
                #plot2dYee(axs[7],  gs, grid, conf, 'ey' , label_title=True)
                #plot2dYee(axs[8],  gs, grid, conf, 'ez' , label_title=True)
                #plot2dYee(axs[9],  gs, grid, conf, 'bx' , label_title=True)
                #plot2dYee(axs[10], gs, grid, conf, 'by' , label_title=True)
                #plot2dYee(axs[11], gs, grid, conf, 'bz' , label_title=True)
                #saveVisz(lap, grid, conf)
    
            #assert that arrays are zero (i.e., scheme is charge conserving)
            self.assertTrue( np.abs(np.sum(gs['jx'])) < 1.e7)
            self.assertTrue( np.abs(np.sum(gs['jy'])) < 1.e7)
            self.assertTrue( np.abs(np.sum(gs['jz'])) < 1.e7)







    def test_problematic_prtcls(self):

        # test pic loop behavior when particle is located in wrong container

        conf = Conf()
        conf.twoD = False
        conf.threeD = True
        conf.Nx = 2
        conf.Ny = 2
        conf.Nz = 2

        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 3

        conf.Nt = 3
        conf.update_bbox()

        grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

        pytools.pic.load_tiles(grid, conf)
        insert_em(grid, conf, zero_field)

        #print("============================================================")

        #--------------------------------------------------
        cid    = grid.id(0,0,0)
        c      = grid.get_tile(cid) #get cell ptr

        container = c.get_container(0) #ispcs
        container.set_keygen_state(0, 0) #number, rank

        u0 = [0.0, 0.0, 0.0]

        # problematic prtcl set 1 (at minimum boundary)
        #x0 = [0.0, 1.0, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 0.0, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 1.0, 0.0]
        #container.add_particle(x0, u0, 1.0)


        # problematic prtcl set 2 (at maximum boundary)
        #x0 = [3.0, 1.0, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 3.0, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 1.0, 3.0]
        #container.add_particle(x0, u0, 1.0)


        # problematic prtcl set 3 (completely outside)
        #x0 = [3.1, 1.0, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 3.1, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 1.0, 3.1]
        #container.add_particle(x0, u0, 1.0)


        #x0 = [-0.1, 1.0, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, -0.1, 1.0]
        #container.add_particle(x0, u0, 1.0)

        #x0 = [1.0, 1.0, -0.1]
        #container.add_particle(x0, u0, 1.0)


        # last tile
        cid    = grid.id(1,1,1)
        c      = grid.get_tile(cid) #get cell ptr
        container = c.get_container(0) #ispcs
        container.set_keygen_state(0, 0) #number, rank

        x0 = [6.0, 5.0, 5.0]
        container.add_particle(x0, u0, 1.0)

        x0 = [0.0, 5.0, 5.0]
        container.add_particle(x0, u0, 1.0)


        #--------------------------------------------------
        pusher   = pyrunko.pic.threeD.BorisPusher()
        fldprop  = pyrunko.emf.threeD.FDTD2()
        fintp    = pyrunko.pic.threeD.LinearInterpolator()
        currint  = pyrunko.pic.threeD.ZigZag()
        flt      = pyrunko.emf.threeD.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        lap = 0
        for lap in range(lap, conf.Nt):

            # --------------------------------------------------
            # push B half
            for tile in pytools.tiles_all(grid):
                fldprop.push_half_b(tile)
    
            # --------------------------------------------------
            # comm B
            grid.send_data(2)
            grid.recv_data(2)
            grid.wait_data(2)
    
            # --------------------------------------------------
            # update boundaries
            for tile in pytools.tiles_all(grid):
                tile.update_boundaries(grid)
    
            ##################################################
            # move particles (only locals tiles)
    
            # --------------------------------------------------
            # interpolate emf
            for tile in pytools.tiles_local(grid):
                fintp.solve(tile)
            #print('successful interp')
    
            # --------------------------------------------------
            # push particles in x and u
            #for tile in pytools.tiles_local(grid):
            #    pusher.solve(tile)
    
            ##################################################
            # advance B half
    
            # --------------------------------------------------
            # push B half
            for tile in pytools.tiles_all(grid):
                fldprop.push_half_b(tile)
    
            # --------------------------------------------------
            # comm B
            grid.send_data(1)
            grid.recv_data(1)
            grid.wait_data(1)
    
            # --------------------------------------------------
            # update boundaries
            for tile in pytools.tiles_all(grid):
                tile.update_boundaries(grid)
    
            ##################################################
            # advance E
    
            # --------------------------------------------------
            # push E
            for tile in pytools.tiles_all(grid):
                fldprop.push_e(tile)
    
            # --------------------------------------------------
            # current calculation; charge conserving current deposition
            for tile in pytools.tiles_local(grid):
                currint.solve(tile)
            #print('successful currint')
    
            # --------------------------------------------------
            # clear virtual current arrays for boundary addition after mpi
            for tile in pytools.tiles_virtual(grid):
                tile.clear_current()
    
            # --------------------------------------------------
            # mpi send currents
            grid.send_data(0)
            grid.recv_data(0)
            grid.wait_data(0)
    
            # --------------------------------------------------
            # exchange currents
            for tile in pytools.tiles_all(grid):
                tile.exchange_currents(grid)
    
            ##################################################
            # particle communication (only local/boundary tiles)
    
            # --------------------------------------------------
            # local particle exchange (independent)
            for tile in pytools.tiles_local(grid):
                tile.check_outgoing_particles()
    
            # --------------------------------------------------
            # global mpi exchange (independent)
            for tile in pytools.tiles_boundary(grid):
                tile.pack_outgoing_particles()
    
            # --------------------------------------------------
            # MPI global particle exchange
            # transfer primary and extra data
            grid.send_data(3)
            grid.recv_data(3)
            grid.wait_data(3)
    
            # orig just after send3
            grid.send_data(4)
            grid.recv_data(4)
            grid.wait_data(4)
    
            # --------------------------------------------------
            # global unpacking (independent)
            for tile in pytools.tiles_virtual(grid):
                tile.unpack_incoming_particles()
                tile.check_outgoing_particles()
    
            # --------------------------------------------------
            # transfer local + global
            for tile in pytools.tiles_local(grid):
                tile.get_incoming_particles(grid)
    
            # --------------------------------------------------
            # delete local transferred particles
            for tile in pytools.tiles_local(grid):
                tile.delete_transferred_particles()
    
            # --------------------------------------------------
            # delete all virtual particles (because new prtcls will come)
            for tile in pytools.tiles_virtual(grid):
                tile.delete_all_particles()
    
            # --------------------------------------------------
            # add current to E
            for tile in pytools.tiles_all(grid):
                tile.deposit_current()
    
            # comm E
            grid.send_data(1)
            grid.recv_data(1)
            grid.wait_data(1)
    
            # --------------------------------------------------
            # comm B
            grid.send_data(2)
            grid.recv_data(2)
            grid.wait_data(2)
    
            # --------------------------------------------------
            # update boundaries
            for tile in pytools.tiles_all(grid):
                tile.update_boundaries(grid)
    

        # --------------------------------------------------
        ip = 0
        prtcls = {}
        for tile in pytools.tiles_all(grid):
            container = tile.get_container(0) #ispcs

            xp = container.loc(0)
            yp = container.loc(1)
            zp = container.loc(2)

            for i in range(len(xp)):
                prtcls[ip] = (xp[i],yp[i],zp[i], tile.cid)
                ip += 1

        #for p in prtcls.keys():
        #    print(prtcls[p])




