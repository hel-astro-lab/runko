# -*- coding: utf-8 -*-

import pyrunko.pic as pypic
import numpy as np
from .tile_initialization import ind2loc


# dummy function that returns always 1 (for weight initialization)
def unit_w(xloc, ispcs, conf):
    return 1.0


# inject plasma into (individual) cells
def inject(grid, 
           vel_func, 
           den_func, 
           conf, 
           align_species=True, 
           w_func=unit_w,
           ):
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

                    running_seed = 1
                    for ispcs in range(conf.Nspecies):
                        container = tile.get_container(ispcs)
                        container.set_keygen_state(prtcl_tot[ispcs], rank)

                        container.type = conf.prtcl_types[ispcs] # name container

                        # open and read previously made particle species (for location reference)
                        if ispcs == 1:
                            #EPS = 1.0e-5
                            ref_container = tile.get_container(0)
                            xxs = ref_container.loc(0) #+ EPS*np.random.rand(1)
                            yys = ref_container.loc(1) #+ EPS*np.random.rand(1)
                            zzs = ref_container.loc(2) #+ EPS*np.random.rand(1)

                            vxs = ref_container.vel(0)
                            vys = ref_container.vel(1)
                            vzs = ref_container.vel(2)

                        ip_mesh = 0

                        #tot_tiles = conf.Nx*conf.Ny*conf.Nz
                        #np.random.seed(k*conf.Nx*conf.Ny + j*conf.Nx + i + ispcs*tot_tiles) # avoid re-starting the rng cycle

                        xD = 1.0 if conf.NxMesh > 1 else 0.0
                        yD = 1.0 if conf.NyMesh > 1 else 0.0
                        zD = 1.0 if conf.NzMesh > 1 else 0.0

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
                                        w = w_func(xloc, ispcs, conf)

                                        if align_species and ispcs == 1:
                                            xx = xxs[ip_mesh] #+ 1.0e-4 #*np.random.rand(1)
                                            yy = yys[ip_mesh] #+ 1.0e-4 #*np.random.rand(1)
                                            zz = zzs[ip_mesh] #+ 1.0e-4 #*np.random.rand(1)
                                            x0 = [xx, yy, zz]  # overwrite location

                                        #this step done in velocity function
                                        #else:
                                        #    xx = xloc[0] + np.random.rand(1)
                                        #    yy = xloc[1] + np.random.rand(1)
                                        #    zz = xloc[2] + np.random.rand(1)
                                        #    x0 = [xx, yy, zz]  # overwrite location

                                        # print("injecting particle sps={} of # {}:th to ({},{},{})".format(
                                        #        ispcs, ip_mesh, x0[0], x0[1], x0[2]))
                                        ip_mesh += 1
                                        container.add_particle(x0, u0, w)
                                        prtcl_tot[ispcs] += 1

                        # less noisy way of injecting particles
                        #totp_per_tile = conf.NxMesh*conf.NyMesh*conf.NzMesh
                        #xloc0 = ind2loc((i, j, k), (0, 0, 0), conf)
                        #ppc = den_func(xloc0, ispcs, conf)
                        #Nip = totp_per_tile*ppc #total number of prtcls of species in this tile

                        #if align_species and ispcs > 0: #use prev container location
                        #    xxps = np.array(xxs)
                        #    yyps = np.array(yys)
                        #    zzps = np.array(zzs)
                        #else: #create new locations
                        #    xxps = xloc0[0] + conf.NxMesh*np.random.rand(Nip)
                        #    yyps = xloc0[1] + conf.NyMesh*np.random.rand(Nip)
                        #    zzps = xloc0[2] + conf.NzMesh*np.random.rand(Nip)

                        #for ip in range(Nip):
                        #    x0 = [xxps[ip], yyps[ip], zzps[ip] ]
                        #    _, u0 = vel_func(xloc0, ispcs, conf) #NOTE reading only the velocity

                        #    ip_mesh += 1
                        #    container.add_particle(x0, u0, 1.0)
                        #    prtcl_tot[ispcs] += 1



    # print("Injected total of:", prtcl_tot)
    return prtcl_tot
