import pycorgi
import pyrunko.ffe as pyffe 

import numpy as np
from initialize_pic import globalIndx, coord2indx, spatialLoc



def initialize_tile(c, i, j, n, conf):

    #initialize tile dimensions 
    c.dx  = 1.0
    c.cfl = conf.cfl

    
    #normalization factors
    ppc = 1
    omp = conf.cfl/conf.c_omp #plasma reaction
    gamma0 = np.sqrt(1.0/(1.0-conf.gamma_e**2.0)) #relativistic dilatation
    q0 = -(gamma0*omp**2.0)/(ppc*(1.0 + np.abs(conf.me/conf.mi)) )

    #set bounding box of the tile 
    mins = spatialLoc(n, [i,j], [0,0,0], conf)
    maxs = spatialLoc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
    c.set_tile_mins(mins[0:2])
    c.set_tile_maxs(maxs[0:2])
    

#load tiles into each node
def loadTiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            if n.get_mpi_grid(i,j) == n.rank():
                c = pyffe.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                initialize_tile(c, i, j, n, conf)
                n.add_tile(c, (i,j)) 



# make all tiles same type 
def initialize_virtuals(n, conf):
    for cid in n.get_virtual_tiles():
        c_orig = n.get_tile(cid)
        (i,j) = c_orig.index

        # TODO: load_metainfo *HAS* to be after add_tile
        c = pyffe.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        n.add_tile(c, (i,j)) 

        #c_orig.communication.local = False;
        c.load_metainfo(c_orig.communication)
        initialize_tile(c, i,j,n, conf)



