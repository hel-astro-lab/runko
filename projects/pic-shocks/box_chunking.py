from __future__ import print_function

import numpy as np
from pytools import ind2loc
import pytools


class Chunker:

    def __init__(self, conf):

        #flag to indicate whether we a still moving
        self.moving = True 

        #speed of light
        self.c = conf.cfl

        #end of the box
        self.Lx = conf.Nx*conf.NxMesh 
        
        #location of the right edge
        self.xmin = 15.0 

        #location of the right wall at current time step
        #self.xmax = conf.Nx*conf.NxMesh

        #reflector speed
        self.betabox = conf.betarefl

        # interval we apply the slider
        self.interval = conf.refl_interval

        # lap when we start sliding the box
        self.t0 = conf.refl_lag

        # TODO check box sliding maximum length requirement


    # get updated mpi grid from each node;
    # grid is up-to-date between nodes because it is always broadcasted 
    def get_mpi_grid(self, grid, conf):
        self.mpi_grid = pytools.get_mpi_grid(grid, conf)
        return self.mpi_grid


    #step injector wall forward
    #def step(self):
    #    if self.moving:
    #        dt = self.c*self.interval # cfl x number of time steps
    #        #print('stepping:', self.wloc0, self.wloc1)


    def slide(self, lap, grid, conf):

        # dont slide before t0 lap
        if lap < self.t0:
            return False

        # step box bounds forward
        self.xmin += self.c*self.betabox
        #self.xmax += self.c*self.betabox

        # slide only on multiples of interval
        if not(lap % self.interval == 0):
            return False

        #--------------------------------------------------
        # passed requirements; chunk!

        #print('injecting stripe between', self.wloc0, self.wloc1)
        for tile in pytools.tiles_all(grid):
            do_delete = False

            i,j,k = pytools.get_index(tile, conf)
            mins = tile.mins
            maxs = tile.maxs

            #slide forward
            #mins[0] += self.xmin
            #maxs[0] += self.xmin

            #tile.set_tile_mins(mins)
            #tile.set_tile_maxs(maxs)

            # if right boundary is behind left boundary of slider
            if maxs[0] < self.xmin:
                do_delete = True

            if do_delete:

                #remove prtcls
                tile.delete_all_particles()

                # reset fields
                #g = tile.get_grids(0)

                #btheta = conf.btheta / 180.0 * np.pi
                #bphi = conf.bphi / 180.0 * np.pi
                #beta = conf.beta

                #for n in range(-3, conf.NzMesh + 3):
                #    for m in range(-3, conf.NyMesh + 3):
                #        for l in range(-3, conf.NxMesh + 3):

                #            g.bx[l, m, n] = conf.binit * np.cos(bphi)
                #            g.by[l, m, n] = conf.binit * np.sin(bphi) * np.sin(btheta)
                #            g.bz[l, m, n] = conf.binit * np.sin(bphi) * np.cos(btheta)

                #            g.ex[l, m, n] = 0.0
                #            g.ey[l, m, n] = -beta * g.bz[l, m, n]
                #            g.ez[l, m, n] = beta * g.by[l, m, n]

        # TODO: broadcast updated mpi grid; done via master rank

        #print("Injected total of:", prtcl_tot)
        return True

