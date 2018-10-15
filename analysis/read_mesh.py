from __future__ import print_function

import numpy as np
import h5py
import os.path

import pycorgi
import pyplasmabox.tools as plasma

from visualize     import imshow
from visualize_amr import get_leaf_mesh
from visualize_amr import xSliceMid, ySliceMid, zSliceMid



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

    length = dset["length"].value.tolist()
    mins   = dset["mins"].value.tolist()
    maxs   = dset["maxs"].value.tolist()
    tref   = dset["top_refinement_level"].value

    # create empty mesh with metadata
    vmesh = plasma.AdaptiveMesh3D()
    vmesh.resize(length)
    vmesh.set_min(mins)
    vmesh.set_max(maxs)
    vmesh.top_refinement_level = tref

    # hashmap values:
    #cids
    #vals
    
    cids = dset["cids"].value
    vals = dset["vals"].value

    # build mesh
    for c, cid in enumerate(cids):
        rfl = vmesh.get_refinement_level(cid)
        indx = vmesh.get_indices(cid)
        uloc = vmesh.get_center(indx, rfl)

        vmesh[indx[0], indx[1], indx[2], rfl] = vals[c]

    if conf.clip:
        vmesh.clip_cells(conf.clipThreshold)

    return vmesh



# fast way of reading 1D meshes directly into vectors
# TODO: implement
def get_1d_mesh(f5, conf):
    #assert if Ny and Nz == 2

    return 0



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


# visualize mesh snapshot
def get_1d_meshes(prefix, lap, conf, spcs, vdir):

    # initialize array
    if vdir == "x":
        fullNvx = np.int(conf.Nvx * (2.0**conf.refinement_level))
        fullNvx = fullNvx if conf.Nvx > 1 else 2*2**conf.refinement_level
        data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, fullNvx) )
    elif vdir == "y":
        fullNvy = np.int(conf.Nvy * (2.0**conf.refinement_level))
        fullNvy = fullNvy if conf.Nvy > 1 else 2*2**conf.refinement_level
        data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, fullNvy) )
    elif vdir == "z":
        fullNvz = np.int(conf.Nvz * (2.0**conf.refinement_level))
        fullNvz = fullNvz if conf.Nvz > 1 else 2*2**conf.refinement_level
        data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, fullNvz) )

    #read from preprocessed file if it exists
    fname = prefix + "meshes_" + str(lap) + "-processed_" + str(spcs) + ".h5"
    if os.path.isfile(fname):
        print("reading from preprocesses file...")
        f5 = h5py.File(fname,'r')

        dset = f5['mesh']
        data = np.array(f5['mesh'][:,:])
        vmin = dset.attrs['vmin']
        vmax = dset.attrs['vmax']

    else:
        print("reading from raw files...")

        # initialize tile info
        tinfo = TileInfo()
        #tinfo.j = 0
        tinfo.j = 0
        tinfo.k = 0
        #tinfo.q = 0
        tinfo.r = 0
        tinfo.s = 0
        tinfo.ispcs = spcs

        rank = 0 #TODO remove hard coded rank
        fname2 = prefix + "meshes-" + str(rank) + "_" + str(lap) + ".h5"
        f5 = h5py.File(fname2,'r')
        #print(fname2)

        # TODO: parallellize
        for i in range(conf.Nx):
            tinfo.i = i
            for s in range(conf.NxMesh):
                tinfo.q = s

                vmesh = get_mesh(f5, tinfo)
                pym = get_leaf_mesh(vmesh, conf.refinement_level)
                
                if vdir == "x":
                    sl = xSliceMid(pym) #slice from middle
                    data[ i*conf.NxMesh + s, :] = sl

                    dx = pym.xx[1]-pym.xx[0]
                    vmin = pym.xx[0] -dx/2.0
                    vmax = pym.xx[-1]+dx/1.0

                elif vdir == "y":
                    sl = ySliceMid(pym) #slice from middle
                    data[ i*conf.NxMesh + s, :] = sl

                    dx = pym.yy[1]-pym.yy[0]
                    vmin = pym.yy[0] -dx/2.0
                    vmax = pym.yy[-1]+dx/1.0
                elif vdir == "z":
                    sl = zSliceMid(pym) #slice from middle

                    dx = pym.zz[1]-pym.zz[0]
                    vmin = pym.zz[0] -dx/2.0
                    vmax = pym.zz[-1]+dx/1.0

                data[ i*conf.NxMesh + s, :] = sl

        #writing the quick read file
        print("writing preprocesses file...")
        f5 = h5py.File(fname,'w')
        dset = f5.create_dataset('mesh', data=data)

        #dset = f5.create_dataset('meta')
        dset.attrs['vmin'] = vmin
        dset.attrs['vmax'] = vmax


    meshes = {
            'data': data,
            'vmin': vmin,
            'vmax': vmax,
             }


    return meshes





# testing if run as main
if __name__ == "__main__":
    #bla

    if False:
        lap = 0
        rank = 0
        fname = "../projects/tests/meshes-" + str(rank) + "_" + str(lap) + ".h5"

        f = h5py.File(fname,'r')

        tinfo = TileInfo()
        tinfo.i = 0
        tinfo.j = 0
        tinfo.k = 0
        
        tinfo.q = 0
        tinfo.r = 0
        tinfo.s = 0

        tinfo.ispcs = 0

        vm = get_mesh(f, tinfo)


    if True:
        from configSetup import Configuration

        import matplotlib.pyplot as plt

        plt.fig = plt.figure(1, figsize=(8,3))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(1, 1)
        gs.update(hspace = 0.5)

        axs = []
        axs.append( plt.subplot(gs[0]) )

        lap = 0
        ispcs = 0
        vdir = "x"
        prefix = "../projects/tests/"

        conf = Configuration(prefix + 'config-landau.ini') 

        meshes = get_1d_meshes(prefix, lap, conf, ispcs, vdir)
        data = meshes['data']
        vmin = meshes['vmin']
        vmax = meshes['vmax']


        ##################################################
        # visualize

        #print(np.max(data))
        data = data/data.max()

        xmin = 0.0
        ymin = 0.0
        xmax = conf.dx*conf.Nx*conf.NxMesh
        ymax = conf.dy*conf.Ny*conf.NyMesh

        imshow(axs[0], data,
               xmin, xmax,
               vmin, vmax,
               cmap = 'plasma_r',
               vmin =   0.0,
               vmax =   1.0,
               clip =   0.0,
               )

        if vdir == "x":
            if ispcs == 0:
                axs[0].set_ylabel(r'$v_{x,e}$')
            if ispcs == 1:
                axs[0].set_ylabel(r'$v_{x,p}$')
        if vdir == "y":
            if ispcs == 0:
                axs[0].set_ylabel(r'$v_{y,e}$')
            if ispcs == 1:
                axs[0].set_ylabel(r'$v_{y,p}$')

        slap = str(lap).rjust(4, '0')
        fname = 'mesh_{}_{}.png'.format(0, slap)
        plt.savefig(fname)






