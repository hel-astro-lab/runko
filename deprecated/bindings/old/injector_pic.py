# -*- coding: utf-8 -*- 
import pycorgi
import pyrunko.pic as pypic

from initialize_pic import spatialLoc

import numpy as np




#inject plasma into (individual) cells
def inject(grid, ffunc, conf, D=2):

    rank = grid.rank()

    prtcl_tot = np.zeros(conf.Nspecies,dtype=np.int64)

    #loop over all *local* cells
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            for k in range(grid.get_Nz()):
                
                if D == 2:
                    pass_flag = grid.get_mpi_grid(i,j) == grid.rank()
                    ind = (i,j)
                elif D == 3:
                    pass_flag = grid.get_mpi_grid(i,j,k) == grid.rank()
                    ind = (i,j,k)

                if pass_flag:
                    #print("creating ({},{})".format(i,j))

                    #get cell & its content
                    if D == 2:
                        cid = grid.id(i,j)
                    if D == 3:
                        cid = grid.id(i,j,k)
                    c = grid.get_tile(cid) #get cell ptr


                    # inject particles
                    # even species are on their own; odd species are located on 
                    # top of previous even ones
                    for ispcs in range(conf.Nspecies):
                        container = c.get_container(ispcs)
                        container.set_keygen_state(prtcl_tot[ispcs], rank)

                        # open and read previously made particle species (for location reference)
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

                                            #vx = vxs[ip_mesh]
                                            #vy = vys[ip_mesh]
                                            #vz = vzs[ip_mesh]
                                            #u0 = [vx, vy, vz] #overwrite velocity

                                        #print("injecting particle sps={} of # {}:th to ({},{},{})".format(
                                        #        ispcs, ip_mesh, x0[0], x0[1], x0[2]))

                                        ip_mesh += 1

                                        container.add_particle(x0, u0, 1.0)
                                        prtcl_tot[ispcs] += 1

    #print("Injected total of:", prtcl_tot)


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(grid, conf, ffunc):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            c = grid.get_tile(i,j)
            gs = c.get_grids(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        # get x_i,j,k
                        xloc0 = spatialLoc(grid, (i,j), (l,m,n), conf)

                        #get x_i+1/2, x_j+1/2, x_k+1/2
                        xloc1 = spatialLoc(grid, (i,j), (l+1,m,  n),   conf)
                        yloc1 = spatialLoc(grid, (i,j), (l,  m+1,n),   conf)
                        zloc1 = spatialLoc(grid, (i,j), (l,  m,  n+1), conf)

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
                        gs.ex[l,m,n] = ffunc(xmid, ycor, zcor)
                        gs.ey[l,m,n] = ffunc(xcor, ymid, zcor)+1.0
                        #gs.ez[l,m,n] = ffunc(xcor, ycor, zmid)+2.0
                        gs.ez[l,m,n] = ffunc(xcor, ycor, zcor)+2.0  #2D hack

                        #gs.bx[l,m,n] = ffunc(xcor, ymid, zmid)+3.0
                        gs.bx[l,m,n] = ffunc(xcor, ymid, zcor)+3.0  #2D hack
                        #gs.by[l,m,n] = ffunc(xmid, ycor, zmid)+4.0 #2D hack
                        gs.by[l,m,n] = ffunc(xmid, ycor, zcor)+4.0
                        gs.bz[l,m,n] = ffunc(xmid, ymid, zcor)+5.0

                        gs.jx[l,m,n] = ffunc(xmid, ymid, zmid)
                        gs.jy[l,m,n] = ffunc(xmid, ymid, zmid)
                        gs.jz[l,m,n] = ffunc(xmid, ymid, zmid)

