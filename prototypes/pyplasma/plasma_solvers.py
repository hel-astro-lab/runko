from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append('../../python')        #plasma, plasmatools
sys.path.append('../../corgi/pycorgi') #corgi mesh infrastucture

import corgi
import plasmatools as ptools
import pyplasma as plasma


from configSetup import Configuration
import initialize as init

from visualize import plotNode
from visualize import plotXmesh
from visualize import plotJ, plotE
from visualize import saveVisz

import injector

from visualize import getYee


def updateBoundaries(node):
    for cid in node.getCellIds():
        c = node.getCellPtr( cid )
        c.updateBoundaries(node)


def save(n, conf, lap, dset):

    #get E field
    yee = getYee(n, conf)

    dset[:, lap] = yee['ex']






if __name__ == "__main__":

    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,7))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(4, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )
    axs.append( plt.subplot(gs[3]) )



    ################################################## 
    #initialize node
    conf = Configuration('config.ini') 

    node = plasma.Grid(conf.Nx, conf.Ny)
    node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

    #node.initMpi()
    #loadMpiXStrides(node)

    init.loadCells(node, conf)


    ################################################## 
    # Path to be created 
    if node.master:
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)


    ################################################## 
    #visualize as a test
    plotNode(axs[0], node, conf)


    ################################################## 
    # test step
    #node.analyzeBoundaryCells()
    #print("{}: send queue        : {}".format(node.rank, node.send_queue))
    #print("{}: send queue address: {}".format(node.rank, node.send_queue_address))

    #node.communicateSendCells()
    #node.communicateRecvCells()
    #plot_node(axs[0], node, 1)


    ################################################## 
    print(" Node howling: {}".format(node.howl()))

    c = node.getCellPtr(1) 
    print(" Cell barking: {}".format(c.bark()))



    ################################################## 
    # initialize
    injector.inject(node, conf) #injecting plasma

    plotXmesh(axs[1], node, conf)
    saveVisz(0, node, conf)



    #setup momentum space solver
    vsol = plasma.MomentumLagrangianSolver()
    intp = ptools.BundleInterpolator4th()
    vsol.setInterpolator(intp)


    #setup spatial space solver
    ssol = plasma.SpatialLagrangianSolver2nd()
    ssol.setGrid(node)


    import h5py
    f = h5py.File("out/run.hdf5", "w")

    grp0 = f.create_group("params")
    grp0.attrs['dx']    = 1.0
    grp0.attrs['dt']    = 0.1



    grp = f.create_group("fields")
    dset = grp.create_dataset("Ex", (conf.Nx*conf.NxMesh, 1000), dtype='f')



    #simulation loop
    for lap in range(1,1000):

        #E field
        #updateBoundaries(node)
        #for cid in node.getCellIds():
        #    c = node.getCellPtr( cid )
        #    c.pushE()

        #momentum step
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                vsol.setCell(cell)
                vsol.solve()

        #spatial step
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                ssol.setTargetCell(i,j)
                ssol.solve()

        #cycle to the new fresh snapshot
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.cyclePlasma()

        #currents
        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            c.depositCurrent()


        #clip every cell
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.clip()


        #I/O
        if (lap % 1 == 0):
            print("--- lap {}".format(lap))
            #plotXmesh(axs[1], node, conf)
            #plotJ(axs[2], node, conf)
            #plotE(axs[3], node, conf)
            #saveVisz(lap, node, conf)

            #save temporarily to file
            save(node, conf, lap, dset)



    #node.finalizeMpi()
