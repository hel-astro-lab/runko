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
from visualize import saveVisz

import injector



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



    ################################################## 
    #initialize node
    conf = Configuration('config.ini') 

    node = plasma.Grid(conf.Nx, conf.Ny)
    node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

    #node.initMpi()
    #loadMpiXStrides(node)

    init.loadCells(node)


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

    c = node.getCell(1) 
    print(" Cell barking: {}".format(c.bark()))



    ################################################## 
    # initialize

    injector.inject(node, conf) #injecting plasma

    plotXmesh(axs[1], node, conf)

    saveVisz(0, node, conf)




    #velocity space solver
    intp = ptools.BundleInterpolator4th()
    vsol = plasma.SplittedLagrangian()
    #vsol.setMesh(mesh)
    vsol.setInterpolator(intp)

    # vsol.solve()

    node.cycle()




    #node.finalizeMpi()
