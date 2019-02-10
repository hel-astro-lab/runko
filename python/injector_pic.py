# -*- coding: utf-8 -*- 
import pycorgi
import pyplasmabox.pic as pypic

from initialize_pic import spatialLoc

import numpy as np




#inject plasma into (individual) cells
def inject(node, ffunc, conf):

    #loop over all *local* cells
    for i in range(node.get_Nx()):
        for j in range(node.get_Ny()):
            if node.get_mpi_grid(i,j) == node.rank():
                #print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = node.id(i,j)
                c      = node.get_tile(cid) #get cell ptr

                #if not(1 <= i <= 2 and j == 1):
                #    continue

                # inject particles
                # even species are on their own; odd species are located on 
                # top of previous even ones
                for ispcs in range(conf.Nspecies):
                    container = c.get_container(ispcs)

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
                                xloc = spatialLoc(node, (i,j), (l,m,n), conf)

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
                                        #u0 = [vx, vy, vz] #overwrite location

                                    #print("injecting particle sps={} of # {}:th to ({},{},{})".format(
                                    #        ispcs, ip_mesh, x0[0], x0[1], x0[2]))

                                    ip_mesh += 1

                                    container.add_particle(x0, u0, 1.0)


#inject test particles into (individual) cells
def inject_test_particles(node, ffunc, conf):

    #loop over all *local* cells
    for i in range(node.get_Nx()):
        for j in range(node.get_Ny()):
            if node.get_mpi_grid(i,j) == node.rank():
                cid    = node.id(i,j)
                c      = node.get_tile(cid) #get cell ptr

                # inject particles
                for ispcs in range(conf.Nspecies_test):
                    container = c.get_test_container(ispcs)

                    #FIXME to inject more generally

                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                xloc = spatialLoc(node, (i,j), (l,m,n), conf)

                                for ip in range(conf.ppc_test):
                                    x0, u0 = ffunc(xloc, ispcs, conf)
                                    container.add_test_particle(x0, u0)





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
