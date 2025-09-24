import pycorgi
import pyrunko.pic as pypic 

import numpy as np


def globalIndx(Ncoords, Mcoords, conf):

    #grid coordinates
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


def spatialLoc(grid, Ncoords, Mcoords, conf):

    #grid coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    l, m, n = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #grid spacing
    xmin = grid.get_xmin()
    ymin = grid.get_ymin()

    dx = 1.0 #conf.dx
    dy = 1.0 #conf.dy
    dz = 1.0 #conf.dz


    #calculate coordinate extent
    x = xmin + i*(NxMesh)*dx + l*dx
    y = ymin + j*(NyMesh)*dy + m*dy
    z = 0.0                  + n*dz

    return [x, y, z]


def initialize_tile(c, indx, n, conf, D=2):

    # try 3d unpacking; if not, fallback to 2D
    try:
        i,j,k = indx
    except:
        i,j = indx

    c.cfl = conf.cfl
    
    ppc = conf.ppc #/ conf.Nspecies

    # load particle containers
    for sps in range(conf.Nspecies):
        if D == 2:
            container = pypic.twoD.ParticleContainer()
        elif D == 3:
            container = pypic.threeD.ParticleContainer()
        if sps % 2 == 0:
            container.q = -conf.qe
        else: 
            container.q = -conf.qi
        
        #reserve memory for particles
        Nprtcls = conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc
        container.reserve(Nprtcls)
    
        c.set_container( container )

    #set bounding box of the tile 
    mins = spatialLoc(n, [i,j], [0,0,0], conf)
    maxs = spatialLoc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
    c.set_tile_mins(mins[0:D])
    c.set_tile_maxs(maxs[0:D])
    
    # initialize analysis tiles ready for incoming simulation data
    #NOTE: only 2D tiles have room for analysis species
    if D <= 2:
        for ip in range(conf.Nspecies):
            c.add_analysis_species()



#load tiles into each grid
def loadTiles(n, conf, D=2):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            for k in range(n.get_Nz()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

                if D == 2:
                    if n.get_mpi_grid(i,j) == n.rank():
                        c = pypic.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                        initialize_tile(c, (i, j, 0), n, conf, D=2)
                        n.add_tile(c, (i,j)) 

                if D == 3:
                    if n.get_mpi_grid(i,j,k) == n.rank():
                        c = pypic.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                        initialize_tile(c, (i, j, k), n, conf, D=3)
                        n.add_tile(c, (i,j,k)) 


# make all tiles same type 
def initialize_virtuals(n, conf, D=2):

    for cid in n.get_virtual_tiles():
        c_orig = n.get_tile(cid)
        ind = c_orig.index

        # new prtcl tile;
        # TODO: load_metainfo *HAS* to be after add_tile

        if D == 2:
            c = pypic.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        if D == 3:
            c = pypic.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        n.add_tile(c, ind) 

        #c_orig.communication.local = False;
        c.load_metainfo(c_orig.communication)
        #print("{}: loading {} owned by {}".format(n.rank(), cid, c.communication.owner))
        
        initialize_tile(c, ind, n, conf, D=D)


