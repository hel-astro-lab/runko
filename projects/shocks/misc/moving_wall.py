from __future__ import print_function

import numpy as np


from initialize_pic import spatialLoc


class MovingWall:

    def __init__(self, conf):
        print("initalizing moving wall")

        #flag to indicate whether we a still moving
        self.moving = True 

        #speed of light
        self.c = conf.cfl

        #end of the box
        self.Lx = conf.Nx*conf.NxMesh 

        #location of the right wall
        self.xloc = 10.0

        #move injector with almost speed of light
        self.betainj = 1.0


    #step injector wall forward
    def step(self):
        if self.moving:
            self.xloc += self.betainj*self.c

        if self.xloc >= self.Lx:
            self.moving = False


    # pre-set tile keygen states
    def set_keygen_state(self, grid):
        rank = grid.rank()
        prtcl_tot = np.zeros(conf.Nspecies,dtype=np.int64)

        #loop over all *local* cells
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                if grid.get_mpi_grid(i,j) == grid.rank():
                    cid    = grid.id(i,j)
                    c      = grid.get_tile(cid) 

                    for ispcs in range(conf.Nspecies):
                        container = c.get_container(ispcs)
                        container.set_keygen_state(prtcl_tot[ispcs], rank)

                        prtcl_tot[ispcs] += conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc


    def inject(self, grid, ffunc, conf):


        rank = grid.rank()
        prtcl_tot = np.zeros(conf.Nspecies,dtype=np.int64)

        #active tile i
        i = self.xloc/conf.NxMesh

        #loop over all *local* cells
        for j in range(grid.get_Ny()):
            if grid.get_mpi_grid(i,j) == grid.rank():
                #print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = grid.id(i,j)
                c      = grid.get_tile(cid) #get cell ptr

                # inject particles
                # even species are on their own; odd species are located on 
                # top of previous even ones
                for ispcs in range(conf.Nspecies):
                    container = c.get_container(ispcs)
                    container.set_keygen_state(prtcl_tot[ispcs], rank)

                    # open and read previously made particle species 
                    #(for location reference)
                    if ispcs % 2 == 1:
                        ref_container = c.get_container(ispcs-1)
                        xxs = ref_container.loc(0)
                        yys = ref_container.loc(1)
                        zzs = ref_container.loc(2)

                        vxs = ref_container.vel(0)
                        vys = ref_container.vel(1)
                        vzs = ref_container.vel(2)

                    ip_mesh = 0
                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))
                                xloc = spatialLoc(grid, (i,j), (l,m,n), conf)

                                for ip in range(conf.ppc):
                                    x0, u0 = ffunc(xloc, ispcs, conf)
                                    if ispcs % 2 == 1:
                                        xx = xxs[ip_mesh]
                                        yy = yys[ip_mesh]
                                        zz = zzs[ip_mesh]
                                        x0 = [xx, yy, zz] #overwrite location

                                    #print("injecting particle sps={} of # {}:th to ({},{},{})".format(
                                    #        ispcs, ip_mesh, x0[0], x0[1], x0[2]))

                                    ip_mesh += 1

                                    container.add_particle(x0, u0, 1.0)
                                    prtcl_tot[ispcs] += 1

        #print("Injected total of:", prtcl_tot)




