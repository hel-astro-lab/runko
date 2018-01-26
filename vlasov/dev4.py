from __future__ import print_function

import numpy as np
from scipy.stats import multivariate_normal


import pyplasma as plasma
import pyplasmaDev as pdev



import initialize as init
import injector



def filler(xloc, uloc, ispcs, conf):

    delgam = np.sqrt(1.0)
    mux = 0.0
    muy = 0.0
    muz = 0.0

    mean = [0.0, 0.0, 1.0]
    cov  = np.zeros((3,3))
    cov[0,0] = 1.0 
    cov[1,1] = 2.0 
    cov[2,2] = 5.0

    f = multivariate_normal.pdf(uloc, mean, cov)

    return f





class Conf:

    outdir = "out"


    #-------------------------------------------------- 
    # space parameters
    Nx = 1
    Ny = 1

    dt = 1.0
    dx = 1.0
    dy = 1.0
    dz = 1.0

    NxMesh = 1
    NyMesh = 1
    NzMesh = 1

    qm = -1.0

    #-------------------------------------------------- 
    # velocity space parameters
    vxmin = -4.0
    vymin = -5.0
    vzmin = -4.0

    vxmax =  4.0
    vymax =  5.0
    vzmax =  4.0

    Nxv = 15
    Nyv = 15
    Nzv = 15

    #vmesh refinement
    refinement_level = 0
    clip = False
    clipThreshold = 1.0e-5






if __name__ == "__main__":


    conf = Conf()

    ################################################## 
    # node configuration
    node = plasma.Grid(conf.Nx, conf.Ny)
    xmin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.dy*conf.Ny*conf.NyMesh

    node.setGridLims(xmin, xmax, ymin, ymax)

    init.loadCells(node, conf)



    ################################################## 
    # load values into cells
    injector.inject(node, filler, conf)












