# -*- coding: utf-8 -*-

import pycorgi 

#all tile generator
def tiles_all(grid):
    for cid in grid.get_tile_ids():
        tile = grid.get_tile(cid)
        yield tile

#local tile generator
def tiles_local(grid):
    for cid in grid.get_local_tiles():
        tile = grid.get_tile(cid)
        yield tile

#virtual tile generator
def tiles_virtual(grid):
    for cid in grid.get_virtual_tiles():
        tile = grid.get_tile(cid)
        yield tile


# mpi boundary tiles
def tiles_boundary(grid):
    for cid in grid.get_boundary_tiles():
        tile = grid.get_tile(cid)
        yield tile




