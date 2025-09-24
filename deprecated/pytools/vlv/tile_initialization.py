# -*- coding: utf-8 -*-

import numpy as np

from .injector import ind2loc



def initialize_tile(c, ind, n, conf):
    i,j = ind

    #initialize tile dimensions 
    c.cfl = conf.cfl


    # initialize analysis tiles ready for incoming simulation data
    #for ip in range(conf.Nspecies):
    #    c.add_analysis_species()

    #set bounding box of the tile 
    mins = ind2loc(n, [i,j], [0,0,0], conf)
    maxs = ind2loc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
    c.set_tile_mins(mins[0:1])
    c.set_tile_maxs(maxs[0:1])
    c.threshold = conf.clipThreshold


#load tiles into each grid
def loadTiles(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

            if n.get_mpi_grid(i,j) == n.rank():
                c = pyrunko.vlv.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                initialize_tile(c, (i, j), n, conf)

                #add it to the grid
                n.add_tile(c, (i,)) 



