import numpy as np

import sys, os
sys.path.append('../../python')        #plasma, plasmatools
sys.path.append('../../corgi/pycorgi') #corgi mesh infrastucture

import corgi
import plasmatools as ptools
import pyplasma as plasma

from initialize import createEmptyVelocityMesh

np.random.seed(0)


def fillMesh(mesh, ispcs, x, y, z, fill_function, conf):
    for k in range(mesh.Nblocks[2]):
        for j in range(mesh.Nblocks[1]):
            for i in range(mesh.Nblocks[0]):
                cid = mesh.getBlockID([i,j,k])
                (ux, uy, uz) = mesh.getCenter( cid )

                fval = fill_function(x, y, z, ux, uy, uz, conf, ispcs)

                mesh[i,j,k] = [fval, fval, fval, fval]



def spatialLoc(n, Ncoords, Mcoords, conf):

    #node coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    s, r, q = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #grid spacing
    xmin = n.getXmin()
    ymin = n.getYmin()

    dx = conf.dx
    dy = conf.dy
    dz = conf.dz


    #calculate coordinate extent
    x = xmin + i*NxMesh*dx + s*dx
    y = ymin + j*NyMesh*dy + r*dy
    z = 0.0                + q*dz

    return (x, y, z)



#inject plasma into cells
def inject(n, fill_function, conf, clip=True):

    #loop over all *local* cells
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:

                #get cell & its content
                cid = n.cellId(i,j)
                c = n.getCellPtr(cid) #get cell ptr

                pgrid0 = c.getPlasmaGrid()
                pgrid1 = c.getNewPlasmaGrid()

                pgrid0.qms = [conf.qmE, conf.qmP]
                pgrid1.qms = [conf.qmE, conf.qmP]

                for q in range(conf.NzMesh):
                    for r in range(conf.NyMesh):
                        for s in range(conf.NxMesh):

                            (x, y, z) = spatialLoc(n, (i,j), (s,r,q), conf)


                            #next create mesh for electron population
                            mesh0 = createEmptyVelocityMesh(conf)
                            fillMesh(mesh0, 0, x, y, z, fill_function, conf)

                            if clip:
                                mesh0.clip()

                            pgrid0.electrons[s,r,q] = mesh0
                            pgrid1.electrons[s,r,q] = mesh0


                            ################################################## 
                            #And another for positrons
                            mesh1 = createEmptyVelocityMesh(conf)
                            fillMesh(mesh1, 1, x, y, z, fill_function, conf)

                            if clip:
                                mesh1.clip()

                            pgrid0.positrons[s,r,q] = mesh1
                            pgrid1.positrons[s,r,q] = mesh1



