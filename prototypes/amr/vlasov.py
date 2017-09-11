import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

from visualize import *
import vmesh



#physical "real" distribution to compare against
def physical_vel(x,y,z):

    mux = 2.0
    muy = 0.0
    muz = 0.0
    sigmax = 4.0
    sigmay = 6.0
    sigmaz = 9.0

    vx = np.exp(-(x-mux)**2 / sigmax**2 )
    vy = np.exp(-(y-muy)**2 / sigmay**2 )
    #vz = np.exp(-(z-muz)**2 / sigmaz**2 )
    vz = 1.0

    return vx*vy*vz
    #return 0.5




def populate_mesh( mesh ):

    for k in range(mesh.Nblocks[2]):
        for j in range(mesh.Nblocks[1]):
            for i in range(mesh.Nblocks[0]):
                cid = mesh.get_block_ID([i,j,k])
                (x,y,z) = mesh.get_center( cid )

                fval = physical_vel(x,y,z)
                mesh[i,j,k] = [fval, fval, fval, fval]
                #print "({},{},{}) = {}".format(i,j,k,fval)


class Params:
    mins = None
    maxs = None
    lens = None




if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )

    
    #cid = 1
    #for zi in range(10):
    #    for yi in range(10):
    #        for xi in range(10):
    #            GID = zi*(10*10) + yi*10 + xi + 1
    #            print "({},{},{}) cid: {} GID: {}".format(xi,yi,zi,cid,GID)
    #            cid += 1
    

    ################################################## 
    # set-up grid
    # xy
    params = Params()
    params.mins = [ -20.0, -20.0, -20.0 ]
    params.maxs = [  20.0,  20.0,  20.0 ]

    mesh = vmesh.vMesh()
    mesh.zFill(params.mins, params.maxs)

    populate_mesh( mesh )

    print "num of blocks:", mesh.number_of_blocks
    #cid = mesh.get_block_ID( [5,5,5] )
    #print "cid:",cid
    #print "cen:",mesh.get_center(cid)
    #print "val:",mesh[5,5,5]

    mesh.clip()
    print "num of blocks:", mesh.number_of_blocks

    #visualize_mesh(axs[0], mesh, params)
    visualize_data(axs[1], mesh, params)


    print "Memory usage: {} Mb ({} cells)".format(mesh.sizeInBytes()/1e6, mesh.number_of_blocks)
    plt.savefig("vlasov_0001.png")





