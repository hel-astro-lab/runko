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


class Conf:

    Nx = 3
    Ny = 3

    NxMesh = 3
    NyMesh = 3

    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

    def __init__(self):
        print("initialized...")





#load cells into each node
def loadCells(n):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            if n.getMpiGrid(i,j) == n.rank:
                c = plasma.PlasmaCell(i, j, n.rank, n.getNx(), n.getNy(), 3, 3)
                #c = plasma.VlasovCell(i, j, n.rank, n.getNx(), n.getNy(), 3, 3)
                n.addCell(c) #TODO load data to cell


def wrap(ii, N):
    if ii < 0:
        ii = N-1
        return N-1
    if ii == N:
        return 0
    return ii

    #while(ii < 0):
    #    ii += N
    #while(ii >= N):
    #    ii -= N
    #return ii





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
    #conf = Configuration('config.ini') 
    conf = Conf()


    node = plasma.Grid(conf.Nx, conf.Ny)
    node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

    #node.initMpi()
    #loadMpiXStrides(node)
    loadCells(node)

    # lets put values into Yee lattice
    val = 1.0
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                c = node.getCellPtr(i,j)
                yee = c.getYee()

                for q in range(conf.NxMesh):
                    for k in range(conf.NyMesh):
                        yee.ex[q,k,0] = val
                        val += 1

    data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh))


    for cid in node.getCellIds():
        c = node.getCellPtr( cid )
        (i, j) = c.index()

        yee = c.getYee()

        for k in range(conf.NyMesh):
            for q in range(conf.NxMesh):
                data[ i*conf.NxMesh + q, j*conf.NyMesh + k ] = yee.ex[q,k,0]
    print(data)

    #update boundaries
    for cid in node.getCellIds():
        c = node.getCellPtr( cid )
        c.updateBoundaries(node)












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
