from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append('../../python')        #plasma, plasmatools
sys.path.append('../../corgi/pycorgi') #corgi mesh infrastucture

import corgi
import plasmatools as ptools
import pyplasma as plasma

from visualize import plotNode
from visualize import imshow
from visualize import saveVisz



class Conf:

    outdir = "out"

    Nx = 8
    Ny = 8
    NxMesh = 25
    NyMesh = 25

    #Nx = 1
    #Ny = 1
    #NxMesh = 200
    #NyMesh = 200

    xmin =-0.5
    xmax = 0.5

    ymin =-0.5
    ymax = 0.5

    def __init__(self):
        print("initialized...")



#load cells into each node
def loadCells(n, conf):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            if n.getMpiGrid(i,j) == n.rank:
                c = plasma.PlasmaCell(i, j, n.rank, n.getNx(), n.getNy(), conf.NxMesh, conf.NyMesh)
                #c = plasma.VlasovCell(i, j, n.rank, n.getNx(), n.getNy(), 3, 3)
                n.addCell(c) #TODO load data to cell

def wrap(ii, N):
    if ii < 0:
        ii = N-1
        return N-1
    if ii == N:
        return 0
    return ii



#Collect stuff from subgrids to global grid using funGetter
def collectMesh(n, conf, funGet):

    NxFull = conf.Nx* conf.NxMesh
    NyFull = conf.Ny* conf.NyMesh

    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh

    data = -1.0 * np.ones( (3, NxFull, NyFull) )


    for cid in n.getCellIds():
        c = n.getCellPtr( cid )
        (i, j) = c.index()
        yee = c.getYee()

        for r in range(NyMesh):
            for q in range(NxMesh):

                vec = funGet(yee, q, r)
                data[:, i*NxMesh + q, j*NyMesh + r ] = vec

    return data


#pick E field vector components from Yee lattice
def getEfield(yee, q, r):
    ex = yee.ex[q,r,0]
    ey = yee.ey[q,r,0]
    ez = yee.ez[q,r,0]

    return ex, ey, ez

#pick B field vector components from Yee lattice
def getBfield(yee, q, r):
    bx = yee.bx[q,r,0]
    by = yee.by[q,r,0]
    bz = yee.bz[q,r,0]

    return bx, by, bz


def plotQuiver(ax, data, conf):

    datax = data[0,:,:].T
    datay = data[1,:,:].T

    xx = np.linspace( node.getXmin(), node.getXmax(), conf.Nx*conf.NxMesh )
    yy = np.linspace( node.getYmin(), node.getYmax(), conf.Ny*conf.NyMesh )
    X, Y = np.meshgrid(xx, yy)

    ix = np.int(0.5*node.getNx())
    iy = np.int(0.5*node.getNy())

    #filter middle out
    #fi = 10
    #datax[ix-fi:ix+fi, ix-fi:iy+fi] = 0.0
    #datay[ix-fi:ix+fi, ix-fi:iy+fi] = 0.0
    #datax[ix-fi:ix+fi, :] = 0.0
    #datay[ix-fi:ix+fi, :] = 0.0

    #normalize
    lmax = np.max( np.sqrt(datax**2 + datay**2) )
    datax = datax / lmax 
    datay = datay / lmax

    sk = 6
    ax.quiver(X[::sk, ::sk], 
              Y[::sk, ::sk], 
              datax[::sk,::sk], 
              datay[::sk,::sk], 
              pivot='tail',
              scale_units='inches',
              scale=2.0,
              )
    


def plotEfield(ax, n, conf):
    data = collectMesh(n, conf, getEfield)

    #magnitude of vectors
    magn = np.sqrt( data[0,:,:]**2 + data[1,:,:]**2 + data[2,:,:]**2 )

    print("maximum E: {}".format(np.max(magn)))

    imshow(ax, 
            #np.log10(magn),
            magn,
            n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
            cmap = 'inferno_r',
            vmin = magn.min(),
            vmax = magn.max(),
            clip = 0.001,
            #vmin = -10.0,
            #vmax =   1.0,
            #clip = -10.0,
            )

    plotQuiver(ax, data, conf)




def plotBfield(ax, n, conf):
    data = collectMesh(n, conf, getBfield)

    #magnitude of vectors
    magn = np.sqrt( data[0,:,:]**2 + data[1,:,:]**2 + data[2,:,:]**2 )
    print("maximum B: {}".format(np.max(magn)))

    imshow(ax, 
            #np.log10(magn),
            magn,
            n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
            cmap = 'inferno_r',
            vmin = magn.min(),
            vmax = magn.max(),
            clip = 0.001,
            #vmin = -10.0,
            #vmax =   1.0,
            #clip = -10.0,
            )

    plotQuiver(ax, data, conf)



def gauss(x, xmid, sig):
    val = np.exp(-((x-xmid)/sig)**2)
    if val < 0.1:
        return 0.0
    return val


def injectRingCurrent(node, conf):

    xx = np.linspace( node.getXmin(), node.getXmax(), conf.Nx*conf.NxMesh )
    yy = np.linspace( node.getYmin(), node.getYmax(), conf.Ny*conf.NyMesh )

    for i in range(node.getNx()):
        for j in range(node.getNy()):
            #if node.getMpiGrid(i,j) == node.rank:
            if True:
                c = node.getCellPtr(i,j)
                yee1 = c.getYee()
                #yee2 = c.getNewYee()

                for q in range(conf.NxMesh):
                    for r in range(conf.NyMesh):

                        x = xx[i*conf.NxMesh + q]
                        y = yy[j*conf.NyMesh + r]

                        #ring current
                        sep = 0.1 #radius of the ring
                        if x < 0.0:
                            val = +(gauss(x,-sep, 0.05)*gauss(y, 0.0, 0.05))
                        else:
                            val = -(gauss(x, sep, 0.05)*gauss(y, 0.0, 0.05))

                        yee1.jz[q,r, 0] += val



#roll current mesh snapshot to newest updated one
def rollSnapshot(node):
    for cid in node.getCellIds():
        c = node.getCellPtr( cid )
        c.cycleYee()


def updateBoundaries(node):
    for cid in node.getCellIds():
        c = node.getCellPtr( cid )
        c.updateBoundaries(node)



if __name__ == "__main__":

    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(4,8))
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
    conf = Conf()

    node = plasma.Grid(conf.Nx, conf.Ny)
    node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

    #node.initMpi()
    #loadMpiXStrides(node)
    loadCells(node, conf)


    #insert initial current into the grid
    injectRingCurrent(node, conf)
    updateBoundaries(node)

    #plot initial condition
    plotNode(axs[0], node, conf)
    plotEfield(axs[1], node, conf)
    saveVisz(0, node, conf)


    #main loop
    for lap in range(1,51):
        print("---lap: {}".format(lap))

        #B field
        updateBoundaries(node)
        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            c.pushHalfB()
            c.pushHalfB()

        #E field
        updateBoundaries(node)
        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            c.pushE()

        #currents
        updateBoundaries(node)
        for cid in node.getCellIds():
            c = node.getCellPtr( cid )
            c.depositCurrent()

        #inject
        injectRingCurrent(node, conf)

        if (lap % 10 == 0):
            #plotNode(axs[0], node, conf)
            plotBfield(axs[0], node, conf)
            plotEfield(axs[1], node, conf)
            saveVisz(lap, node, conf)



    #node.finalizeMpi()
