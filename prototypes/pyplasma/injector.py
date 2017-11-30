import numpy as np

import sys, os
sys.path.append('../../python')        #plasma, plasmatools
sys.path.append('../../corgi/pycorgi') #corgi mesh infrastucture

import corgi
import plasmatools as ptools
import pyplasma as plasma

from initialize import createEmptyVelocityMesh


# Thermal (Gaussian) plasma with some bulk velocity
def thermalPlasma(vx, vy, vz,
                  Tx, Ty, Tz,
                  Gx, Gy, Gz
                  ):

    #Brownian noise
    sigma = 0.1
    z1 = np.random.standard_normal() 
    z2 = np.random.standard_normal() 
    z3 = np.random.standard_normal() 

    f  = 1.0
    f *= np.exp(-(vx - Gx)**2 / Tx + sigma*z1)
    f *= np.exp(-(vy - Gy)**2 / Ty + sigma*z2)
    f *= np.exp(-(vz - Gz)**2 / Tz + sigma*z3)

    return f


# isotropic 1D (x-dir) thermal plasma
def thermalXPlasma(vx, T, G):
    return thermalPlasma(vx, 0.0, 0.0,
                         T,  T,   T,
                         G,  0.0, 0.0)



def fillMesh(mesh, pl):
    for k in range(mesh.Nblocks[2]):
        for j in range(mesh.Nblocks[1]):
            for i in range(mesh.Nblocks[0]):
                cid = mesh.getBlockID([i,j,k])
                (vx,vy,vz) = mesh.getCenter( cid )

                fval = thermalXPlasma(vx, pl["T"], pl["bulkVelo"])
                mesh[i,j,k] = [fval, fval, fval, fval]



#inject plasma into cells
def inject(n, conf):

    #loop over all *local* cells
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                mesh = createEmptyVelocityMesh(conf)

                if i == 10:
                    pl = { "T": 4.0, "bulkVelo" : 5.0, }
                else:
                    pl = { "T": 2.0, "bulkVelo" : 2.0, }


                fillMesh(mesh, pl)

                cid = n.cellId(i,j)
                c = n.getCellPtr(cid) #get cell ptr

                mesh.clip()

                c.addData(mesh)
                c.addData(mesh)




