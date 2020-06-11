# -*- coding: utf-8 -*- 

import pyrunko.ffe as pyffe


def ind2loc(gridI, tileI, conf):

    # grid coordinates
    i, j, k = gridI
    Nx = conf.Nx
    Ny = conf.Ny
    Nz = conf.Nz

    # tile coordinates
    l, m, n = tileI
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    # grid spacing; start point + step
    xmin = conf.xmin
    ymin = conf.ymin
    zmin = conf.zmin

    dx = 1.0  # conf.dx
    dy = 1.0  # conf.dy
    dz = 1.0  # conf.dz

    # calculate coordinate extent
    x = xmin + i * (NxMesh) * dx + l * dx
    y = ymin + j * (NyMesh) * dy + m * dy
    z = zmin + k * (NzMesh) * dz + n * dz

    return [x, y, z]



def initialize_tile(tile, indx, n, conf):

    # set parameters
    tile.cfl = conf.cfl
    ppc = conf.ppc  # / conf.Nspecies

    # set bounding box of the tile
    mins = ind2loc(indx, (0, 0, 0), conf)
    maxs = ind2loc(indx, (conf.NxMesh, conf.NyMesh, conf.NzMesh), conf)

    if conf.threeD:
        tile.set_tile_mins(mins[0:3])
        tile.set_tile_maxs(maxs[0:3])
    elif conf.twoD:
        tile.set_tile_mins(mins[0:2])
        tile.set_tile_maxs(maxs[0:2])

    return


# load virtual tiles
def load_virtual_tiles(n, conf):

    for cid in n.get_virtual_tiles():
        tile_orig = n.get_tile(cid)
        ind = tile_orig.index

        # new prtcl tile;
        # TODO: load_metainfo *HAS* to be after add_tile because
        # add_tile modifies tile content.

        if conf.threeD:
            i,j,k = ind
            tile = pyffe.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

            n.add_tile(tile, ind) 
            tile.load_metainfo(tile_orig.communication)
            initialize_tile(tile, (i,j,k), n, conf)

        elif conf.twoD:
            i,j = ind
            tile = pyffe.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

            n.add_tile(tile, ind) 
            tile.load_metainfo(tile_orig.communication)
            initialize_tile(tile, (i,j,0), n, conf)

    return 



# 3D loading of ffe tiles into grid
def load_tiles(n, conf):

    for k in range(n.get_Nz()):
        for j in range(n.get_Ny()):
            for i in range(n.get_Nx()):
                # print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

                if conf.threeD:
                    if n.get_mpi_grid(i, j, k) == n.rank():
                        tile = pyffe.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                        ind = (i, j, k)
                        initialize_tile(tile, (i,j,k), n, conf)
                        n.add_tile(tile, ind)

                elif conf.twoD:
                    if n.get_mpi_grid(i, j) == n.rank():
                        tile = pyffe.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                        ind = (i, j)

                        initialize_tile(tile, (i,j,k), n, conf)
                        n.add_tile(tile, ind)

    return

