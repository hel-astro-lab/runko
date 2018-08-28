import pycorgi
import pyplasmabox.pic as pypic 

import numpy as np


def globalIndx(Ncoords, Mcoords, conf):

    #node coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    l, m, n = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #calculate coordinate extent
    x = 1.0*(i*NxMesh + l)
    y = 1.0*(j*NyMesh + m)
    z = 1.0*(0        + n)

    return [x, y, z]


def coord2indx(xloc, conf):
    x,y,z = xloc

    i = x/(1.0*conf.Nx*conf.NxMesh)
    j = y/(1.0*conf.Ny*conf.NyMesh)
    k = z/(1.0*conf.Nz*conf.NzMesh)

    return [i,j,k]


def spatialLoc(node, Ncoords, Mcoords, conf):

    #node coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    l, m, n = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #grid spacing
    xmin = node.getXmin()
    ymin = node.getYmin()

    dx = 1.0 #conf.dx
    dy = 1.0 #conf.dy
    dz = 1.0 #conf.dz


    #calculate coordinate extent
    x = xmin + i*(NxMesh)*dx + l*dx
    y = ymin + j*(NyMesh)*dy + m*dy
    z = 0.0                  + n*dz

    return [x, y, z]


def initialize_tile(c, i, j, n, conf):

    #initialize tile dimensions 
    #c.dt  = conf.dt
    #c.dx  = conf.dx
    c.dx  = 1.0
    c.cfl = conf.cfl
    
    ppc = conf.ppc #/ conf.Nspecies
    
    #normalization factors
    omp = conf.cfl/conf.c_omp #plasma reaction
    #gamma0 = 1.0      #relativistic dilatation
    gamma0 = np.sqrt(1.0/(1.0-conf.gamma_e**2.0)) #relativistic dilatation
    #betaN = np.sqrt(1.0 - 1.0/gamma0**2.0)
    q0 = -(gamma0*omp**2.0)/(ppc*(1.0 + np.abs(conf.me/conf.mi)) )
    #print("normalization factor: {}".format(q0))
    
    # load particle containers
    for sps in range(conf.Nspecies):
        container = pypic.ParticleBlock(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        if sps % 2 == 0:
            container.q = conf.me*q0
        else:
            container.q = conf.mi*q0
        
        #reserve memory for particles
        Nprtcls = conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc
        container.reserve(Nprtcls, 3)
        container.resizeEM(Nprtcls, 3)
    
        c.set_container( container )
    
    
    #set bounding box of the tile 
    mins = spatialLoc(n, [i,j], [0,0,0], conf)
    maxs = spatialLoc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
    c.set_tile_mins(mins[0:2])
    c.set_tile_maxs(maxs[0:2])
    
    
    # initialize analysis tiles ready for incoming simulation data
    for ip in range(conf.Nspecies):
        c.addAnalysisSpecies()



#load tiles into each node
def loadTiles(n, conf):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.getMpiGrid(i,j), ref[j,i]))

            if n.getMpiGrid(i,j) == n.rank:
                c = pypic.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                
                initialize_tile(c, i, j, n, conf)

                #add it to the node
                n.addTile(c, (i,j)) 





