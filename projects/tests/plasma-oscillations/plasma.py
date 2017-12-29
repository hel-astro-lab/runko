from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os

import corgi
import plasmatools as ptools
import pyplasma as plasma


from configSetup import Configuration
import initialize as init

from visualize import plotNode
from visualize import plotXmesh
from visualize import plotJ, plotE
from visualize import saveVisz
from visualize import getYee

import injector


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
    plt.fig = plt.figure(1, figsize=(8,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(5, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )
    axs.append( plt.subplot(gs[3]) )
    axs.append( plt.subplot(gs[4]) )



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
    #if node.master:
    if True:
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)



    ################################################## 
    # initialize
    injector.inject(node, conf) #injecting plasma

    plotNode(axs[0], node, conf)
    plotXmesh(axs[1], node, conf, 0)
    plotXmesh(axs[2], node, conf, 1)
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
    grp0.attrs['dx']    = 0.1
    grp0.attrs['dt']    = 0.01


    Nt = 50 

    grp = f.create_group("fields")
    dset = grp.create_dataset("Ex", (conf.Nx*conf.NxMesh, Nt), dtype='f')



    #simulation loop
    for lap in range(1,Nt):

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

            #save temporarily to file
            save(node, conf, lap, dset)


        if (lap % 1 == 0):
            plotNode(axs[0], node, conf)
            plotXmesh(axs[1], node, conf, 0) #electrons
            plotXmesh(axs[2], node, conf, 1) #positrons

            plotJ(axs[3], node, conf)
            plotE(axs[4], node, conf)
            saveVisz(lap, node, conf)


    #node.finalizeMpi()
