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

#import initialize as init #lets redo this ourselves and add maxwell cells

from visualize import plotNode
from visualize import saveVisz

import injector




#load cells into each node
def loadCells(n):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            if n.getMpiGrid(i,j) == n.rank:
                c = plasma.PlasmaCell(i, j, n.rank, n.getNx(), n.getNy())
                #c = plasma.VlasovCell(i, j, n.rank, n.getNx(), n.getNy())
                n.addCell(c) #TODO load data to cell





if __name__ == "__main__":

    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(2, 1)
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
    loadCells(node)

    for i in range(node.getNx()):
        for j in range(node.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                print("ij ({},{})".format(i,j))
                c = node.getCellPtr(i,j)
                c.pushE()
                c.updateBoundaries(node)

                

    #c = plasma.PlasmaCell(0, 0, 0, 10, 10)
    #c.pushE()
    #c.pushHalfB()
    #c.updateBoundaries(node)




    #visualize as a test
    #plotNode(axs[0], node, conf)

    ################################################## 
    #print(" Node howling: {}".format(node.howl()))
    #c = node.getCellPtr(1) 
    #c.pushE()
    #c.pushHalfB()

    #plotXmesh(axs[1], node, conf)
    #saveVisz(0, node, conf)


    #node.finalizeMpi()
