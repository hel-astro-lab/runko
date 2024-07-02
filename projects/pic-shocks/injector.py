from __future__ import print_function

from mpi4py import MPI

import numpy as np
from pytools import ind2loc
import pytools


class MovingInjector:

    def __init__(self, conf):

        #flag to indicate whether we a still moving
        self.moving = True 

        #speed of light
        self.c = conf.cfl

        #end of the box
        self.Lx = conf.Nx*conf.NxMesh 
        
        #location of the right wall at previous time step
        self.wloc0 = 15.0 #reflecting wall is always set to 15th cell

        #location of the right wall at current time step
        self.wloc1 = self.wloc0 + conf.inj_startx*conf.c_omp #in skindepth units

        #move injector with almost speed of light
        self.betainj = conf.betainj

        # plasma bulk motion
        self.betaplasma = conf.beta

        # interval we apply the wall
        self.interval = conf.inj_interval

        # width of the EM field damping region; keep it modest
        self.damping_region_width = conf.inj_damping_region

        # finally, correct head position accounting for the plasma motion
        self.wloc1 -= self.betaplasma*self.c*(1*self.interval)

        # get the rank that this injector is being created with
        self.rank = MPI.COMM_WORLD.Get_rank()

        # number of particles injected locally
        self.prtcl_tot = np.zeros(conf.Nspecies, dtype=np.int64)


    #step injector wall forward
    def step(self, lap):
        
        # only step forward on interval periods
        if not(lap % self.interval == 0):
            return

        if self.wloc1 >= self.Lx-10.0:
            #print('Injector stopping at ', self.wloc0, self.wloc1)
            self.moving = False

        if self.moving:
            dt = self.c*self.interval # cfl x number of time steps

            #print('stepping:', self.wloc0, self.wloc1)

            #head of injector
            self.wloc0 = 1.0*self.wloc1 #- self.betaplasma*self.c*(self.interval+1.0)

            # tail of injector
            self.wloc1 += self.betainj*dt

            #correct for plasma bulk motion; 
            # start injecting earlier to close the gap between escaping plasma and injector
            self.wloc0 -= self.betaplasma*dt

            #print('stepping after:', self.wloc0, self.wloc1)



    # TODO: split this into tile-based operations to better reflect standard API
    def inject(self, grid, lap, vel_func, den_func, conf):

        # do not inject if we are no longer moving
        if not(self.moving):
            return 0

        # only inject with modulo interval
        if not(lap % self.interval == 0):
            return 0

        #print('injecting stripe between', self.wloc0, self.wloc1)

        for tile in pytools.tiles_local(grid):
            i,j,k = pytools.get_index(tile, conf)

            # inject particles
            # even species are on their own; odd species are located on
            # top of previous even ones

            tile_xmin = tile.mins[0]
            tile_xmax = tile.maxs[0]

            # load containers
            container1 = tile.get_container(0)
            container2 = tile.get_container(1)

            container1.set_keygen_state(self.prtcl_tot[0], self.rank)
            container2.set_keygen_state(self.prtcl_tot[1], self.rank)

            for l in range(conf.NxMesh):
                xloc0 = ind2loc((i, j, k), (l, 0, 0), conf)

                # injection stripe
                x0 = np.floor(self.wloc0) 
                x1 = np.ceil(self.wloc1)

                if x0 <= xloc0[0] <= x1:
                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            # print(" sub mesh: ({},{},{})".format(l,m,n))
                            xloc = ind2loc((i, j, k), (l, m, n), conf)

                            # calculate how many electron prtcls to inject in this loc
                            ppc = den_func(xloc, 0, conf)

                            for ip in range(ppc):
                                xl0, u0 = vel_func(xloc, 0, conf) #sample electron loc and vel

                                if self.wloc0 <= xl0[0] <= self.wloc1:

                                    #check that we are still inside tile limits
                                    # TODO; is this needed?
                                    if tile_xmin <= xl0[0] <= tile_xmax:

                                        #inject electrons
                                        container1.add_particle(xl0, u0, 1.0)

                                        # inject positrons with new velocity but same position
                                        xl1, u1 = vel_func(xloc, 1, conf)
                                        container2.add_particle(xl0, u1, 1.0)

                                        # increase local counter for book keeping
                                        self.prtcl_tot[0] += 1
                                        self.prtcl_tot[1] += 1


        #print("Injected total of:", prtcl_tot)
        return self.prtcl_tot


    # reset EM fields ahead of the piston
    # TODO: split this into tile-based operations to better reflect standard API
    def damp_em_fields(self, grid, lap, conf):

        #damping region is set to end at the start of the injection location
        # NOTE: this is corrected to account for the plasma motion during
        # the `interval` period we actually apply the injection 

        dt = self.c # cfl x number of time steps
        dt *= lap % self.interval  #time steps taken after injection

        #correct for plasma bulk motion; 
        # start damping earlier to close the gap between escaping plasma and injector
        d0 = int(self.wloc0 - self.betaplasma*dt - self.damping_region_width - 5.0)
        d1 = int(self.wloc0 - self.betaplasma*dt - 5.0 )

        # do not damp the early region when shock is still forming
        if d0 < self.wloc0:
            return

        # loop over local tiles and apply damping if tile crosses damping zone
        for tile in pytools.tiles_local(grid):

            # get index
            ind = tile.index
            if conf.threeD:
                i,j,k = ind
            elif conf.twoD:
                i,j = ind
                k = 0
            elif conf.oneD:
                i = ind[0]
                j, k = 0

            tile_xmin = tile.mins[0]
            tile_xmax = tile.maxs[0]
            
            # gloabl position of the  tile right and tile left boundaries
            xloc = ind2loc((i, j, k), (0, 0, 0), conf)

            # start and end of the active tile
            #xg0 = (i-1)*conf.NxMesh #pad with extra tiles to left of injector
            #xg1 = (i+2)*conf.NxMesh #pad with extra tiles to right of injector

            # check if tile x region is inside damping zone 
            #if not(xg1 < d0 or xg0 > d1):

            # pad damping region with extra tiles around the modified zone
            xg0 = d0 - 2*conf.NxMesh - 5.0
            xg1 = d1 + 3*conf.NxMesh + 5.0

            #   xg0      xg1
            #    |        |  
            #   
            # |     |             case 1
            # xmin  xmax
            #      |     |        case 2
            #      xmin  xmax
            #          |     |    case 3
            #          xmin  xmax
            #if tile_xmin > xg0 and tile_xmax < xg1:
            if tile_xmax > xg0 and tile_xmin < xg1:
                tile.fld1 = d1
                tile.fld2 = d0
                tile.damp_fields()

        return 


