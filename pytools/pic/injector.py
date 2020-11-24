# -*- coding: utf-8 -*-

import pyrunko.pic as pypic
import numpy as np
from .tile_initialization import ind2loc


# inject plasma into (individual) cells
def inject(grid, vel_func, den_func, conf):
    rank = grid.rank()

    prtcl_tot = np.zeros(conf.Nspecies, dtype=np.int64)

    # loop over all *local* cells
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            for k in range(grid.get_Nz()):

                # check that we own this tile
                if conf.threeD:
                    pass_flag = grid.get_mpi_grid(i, j, k) == grid.rank()
                elif conf.twoD:
                    pass_flag = grid.get_mpi_grid(i, j) == grid.rank()

                if pass_flag:
                    # print("creating ({},{})".format(i,j))

                    # get cell & its content
                    if conf.threeD: cid  = grid.id(i, j, k)
                    elif conf.twoD: cid  = grid.id(i, j)

                    tile = grid.get_tile(cid)  # get cell ptr

                    # inject particles
                    # even species are on their own; odd species are located on
                    # top of previous even ones
                    for ispcs in range(conf.Nspecies):
                        container = tile.get_container(ispcs)
                        container.set_keygen_state(prtcl_tot[ispcs], rank)

                        # open and read previously made particle species (for location reference)
                        if ispcs > 0:
                            ref_container = tile.get_container(0)
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

                                    # print(" sub mesh: ({},{},{})".format(l,m,n))
                                    xloc = ind2loc((i, j, k), (l, m, n), conf)
                                    
                                    # calculate how many prtcls per species to inject in this loc
                                    ppc = den_func(xloc, ispcs, conf)

                                    #TODO: pair loading factor

                                    for ip in range(ppc):
                                        x0, u0 = vel_func(xloc, ispcs, conf)
                                        if ispcs > 0:
                                            xx = xxs[ip_mesh]
                                            yy = yys[ip_mesh]
                                            zz = zzs[ip_mesh]
                                            x0 = [xx, yy, zz]  # overwrite location

                                        # print("injecting particle sps={} of # {}:th to ({},{},{})".format(
                                        #        ispcs, ip_mesh, x0[0], x0[1], x0[2]))

                                        ip_mesh += 1

                                        container.add_particle(x0, u0, 1.0)
                                        prtcl_tot[ispcs] += 1

    # print("Injected total of:", prtcl_tot)
    return prtcl_tot
