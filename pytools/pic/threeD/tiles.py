import pyrunko.pic as pypic


def ind2loc(grid, gridI, tileI, conf):

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

    # grid spacing
    xmin = grid.get_xmin()
    ymin = grid.get_ymin()
    zmin = grid.get_zmin()

    dx = 1.0  # conf.dx
    dy = 1.0  # conf.dy
    dz = 1.0  # conf.dz

    # calculate coordinate extent
    x = xmin + i * (NxMesh) * dx + l * dx
    y = ymin + j * (NyMesh) * dy + m * dy
    z = zmin + k * (NzMesh) * dz + n * dz

    return [x, y, z]


def initialize_tile(tile, indx, n, conf):
    i, j, k = indx

    # set parameters
    tile.cfl = conf.cfl
    ppc = conf.ppc  # / conf.Nspecies

    # load particle containers
    for sps in range(conf.Nspecies):
        container = pypic.threeD.ParticleContainer()

        # alternate injection between - and + charged prtcls
        if sps % 2 == 0:
            container.q = -conf.qe
        else:
            container.q = -conf.qi

        # reserve memory for particles
        Nprtcls = conf.NxMesh * conf.NyMesh * conf.NzMesh * conf.ppc
        container.reserve(Nprtcls)

        tile.set_container(container)

    # set bounding box of the tile
    mins = ind2loc(n, [i, j, k], [0, 0, 0], conf)
    maxs = ind2loc(n, [i, j, k], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)

    tile.set_tile_mins(mins[0:3])
    tile.set_tile_maxs(maxs[0:3])

    # initialize analysis tiles ready for incoming simulation data
    # NOTE: only 2D tiles have room for analysis species
    # for ip in range(conf.Nspecies):
    #    c.add_analysis_species()

    return


# 3D loading of pic tiles into grid
def load_tiles(n, conf):

    for k in range(n.get_Nz()):
        for j in range(n.get_Ny()):
            for i in range(n.get_Nx()):
                # print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))

                if n.get_mpi_grid(i, j, k) == n.rank():
                    tile = pypic.threeD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                    ind = (i, j, k)
                    initialize_tile(tile, ind, n, conf)

                    n.add_tile(tile, ind)

    return
