# -*- coding: utf-8 -*-

import numpy as np
import pyrunko


# small data struct that stores mesh location info
class TileInfo:

    i = 0
    j = 0
    k = 0

    q = 0
    r = 0
    s = 0

    ispcs = 0

    clip = True
    clipThreshold = 0.0


def get_mesh(f5, conf):

    # parse location
    # 
    # format is:
    #
    # GROUP "tile-9_0_0" {
    #    GROUP "loc-0_0_0" {
    #       GROUP "sp-0" {
    #          DATASET "cids" {
    # etc.

    i = conf.i
    j = conf.j
    k = conf.k

    q = conf.q
    r = conf.r
    s = conf.s

    ip = conf.ispcs

    tilename = "tile-"+str(i)+"_"+str(j)+"_"+str(k)
    locname  = "loc-" +str(q)+"_"+str(r)+"_"+str(s)
    ispname  = "sp-"  +str(ip)

    #print(tilename)
    #print(locname)
    #print(ispname)
    dset_tile = f5[tilename]
    dset_loc  = dset_tile[locname]
    dset      = dset_loc[ispname]

    #meta info about grid:
    # maximum_refinement_level
    # error_cid
    # error_index
    # top_refinement_level
    # length
    # mins
    # maxs
    # number_of_cells

    length = dset["length"][()]
    mins   = dset["mins"][()]
    maxs   = dset["maxs"][()]
    tref   = dset["top_refinement_level"][()]

    # create empty mesh with metadata
    vmesh = pyrunko.tools.AdaptiveMesh3D()
    vmesh.resize(length)
    vmesh.set_min(mins)
    vmesh.set_max(maxs)
    vmesh.top_refinement_level = tref

    # hashmap values:
    #cids
    #vals
    
    cids = dset["cids"][()]
    vals = dset["vals"][()]

    # build mesh
    for c, cid in enumerate(cids):
        rfl = vmesh.get_refinement_level(cid)
        indx = vmesh.get_indices(cid)
        uloc = vmesh.get_center(indx, rfl)

        vmesh[indx[0], indx[1], indx[2], rfl] = vals[c]

    if conf.clip:
        vmesh.clip_cells(conf.clipThreshold)

    return vmesh


