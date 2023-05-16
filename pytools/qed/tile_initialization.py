# -*- coding: utf-8 -*- 

import pyrunko.qed as pyqed


def initialize_tile(tile, indx, n, conf):

    # set parameters
    tile.cfl = conf.cfl

    #--------------------------------------------------
    # particle container

    #ppc = conf.ppc  # / conf.Nspecies

    # load particle containers
    #for sps in range(conf.Nspecies):

    #    if conf.threeD:
    #        container = pypic.threeD.ParticleContainer()
    #    elif conf.twoD:
    #        container = pypic.twoD.ParticleContainer()

    #    # alternate injection between - and + charged prtcls
    #    if sps % 2 == 0:
    #        container.q = -conf.qe
    #    else:
    #        container.q = -conf.qi

    #    # reserve memory for particles
    #    Nprtcls = conf.NxMesh * conf.NyMesh * conf.NzMesh * conf.ppc
    #    container.reserve(Nprtcls)

    #    tile.set_container(container)

    #--------------------------------------------------
    # photon buckets

    # load particle containers
    bucket = pyrad.PhotonContainer()

    #reserve memory for particles
    Nprtcls = conf.ppt
    bucket.reserve(Nprtcls)
    tile.push_back( bucket )

    #--------------------------------------------------

    # set bounding box of the tile
    #mins = ind2loc(indx, (0, 0, 0), conf)
    #maxs = ind2loc(indx, (conf.NxMesh, conf.NyMesh, conf.NzMesh), conf)

    #if conf.threeD:
    #    tile.set_tile_mins(mins[0:3])
    #    tile.set_tile_maxs(maxs[0:3])
    #elif conf.twoD:
    #    tile.set_tile_mins(mins[0:2])
    #    tile.set_tile_maxs(maxs[0:2])


    return








