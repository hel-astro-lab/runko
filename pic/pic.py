from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os

import corgi
import pyplasma as plasma
import pypic 


from configSetup import Configuration
import initialize as init

from injector import spatialLoc



from visualize import plotNode
from visualize import plotJ, plotE, plotDens
from visualize import getYee
from visualize import saveVisz


from timer import Timer





#load cells into each node
def loadCells(n, conf):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.getMpiGrid(i,j), ref[j,i]))

            if n.getMpiGrid(i,j) == n.rank:
                c = pypic.PicCell(i, j, n.rank, 
                                  n.getNx(), n.getNy(),
                                  conf.NxMesh, conf.NyMesh
                                  )

                #initialize cell dimensions 
                c.dt = conf.dt
                c.dx = conf.dx

                # use same scale for Maxwell solver
                c.yeeDt = conf.dt
                c.yeeDx = conf.dx

                #add it to the node
                n.addCell(c) 




#inject plasma into cells
def inject(node, ffunc, conf):

    #loop over all *local* cells
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = node.cellId(i,j)
                c      = node.getCellPtr(cid) #get cell ptr

                #reserve memory for particles
                Nprtcls = conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc
                c.container.reserve(Nprtcls, 3)


                for ispcs in range(conf.Nspecies):

                    #set q/m
                    if ispcs == 0:
                        qm = conf.me
                    elif ispcs == 1:
                        qm = conf.mi

                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))
                                xloc = spatialLoc(node, (i,j), (l,m,n), conf)

                                for ip in range(conf.ppc):
                                    #xloc = [0.0, 0.0, 0.0]
                                    uloc = [0.1, 0.0, 0.0]

                                    c.container.add_particle(xloc, uloc)


# visualize particle content in x-dir
def plotXmesh(ax, n, conf, spcs, vdir):

    for i in range(conf.Nx):
        cid = n.cellId(i,0)
        c = n.getCellPtr(cid)

        x = c.container.loc(0)
        y = c.container.loc(1)
        z = c.container.loc(2)

        ux = c.container.vel(0)
        uy = c.container.vel(1)
        uz = c.container.vel(2)

        ax.plot(x, ux, "k.", markersize=2, alpha=0.8)
        
    if vdir == "x":
        if spcs == 0:
            ax.set_ylabel(r'$v_{x,e}$')
        if spcs == 1:
            ax.set_ylabel(r'$v_{x,p}$')

    ax.set_xlim(n.getXmin(), n.getXmax())
    ax.set_ylim(-1.0, 1.0)



def filler(xloc, uloc, ispcs, conf):
    x = 10.0
    return x










if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(8, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )


    conf = Configuration('config.ini') 


    #node = plasma.Grid(conf.Nx, conf.Ny)
    node = corgi.Node(conf.Nx, conf.Ny)

    xmin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.dy*conf.Ny*conf.NyMesh

    node.setGridLims(xmin, xmax, ymin, ymax)

    loadCells(node, conf)


    ################################################## 
    # Path to be created 
    #if node.master:
    if True:
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)



    inject(node, filler, conf) #injecting plasma particles


    # visualize initial condition
    plotNode( axs[0], node, conf)
    plotXmesh(axs[1], node, conf, 0, "x")

    saveVisz(-1, node, conf)







