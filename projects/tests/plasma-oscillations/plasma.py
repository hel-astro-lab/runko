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

from timer import Timer


# Generic function to fill the velocity mesh
#
# Maxwellian plasma with Brownian noise
# where delgam = kT/m_i c^2
#
def filler(x, y, z, ux, uy, uz, conf, ispcs):

    #electrons
    if ispcs == 0:
        delgam  = conf.delgam * np.abs(conf.mi / conf.me) * conf.temperature_ratio

        # bulk velocities
        mux = conf.gamma_e
        muy = 0.0
        muz = 0.0


        #give electrons a nudge 
        Lx  = conf.Nx*conf.NxMesh*conf.dx
        kx = 1
        mux += conf.beta*np.sin( 2.0*np.pi*kx*x/Lx)


    #ions/positrons
    elif ispcs == 1:
        delgam  = conf.delgam

        # bulk velocities
        mux = conf.gamma_i
        muy = 0.0
        muz = 0.0

        Lx  = conf.Nx*conf.NxMesh*conf.dx
        kx = 1
        mux += conf.beta*np.sin( 2.0*np.pi*kx*x/Lx)


    #Classical Maxwellian distribution
    z1 = 0.01*np.random.standard_normal() # Brownian noise
    #z1 = 0.0

    f  = 1.0/np.sqrt(2.0*np.pi*delgam)
    f *= np.exp(-0.5*((ux - mux)**2)/delgam + z1*delgam)

    return f



# Get Yee grid components from node and save to hdf5 file
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


    # Timer for profiling
    timer = Timer(["total", "init", "step", "io"])
    timer.start("total")
    timer.start("init")


    ################################################## 
    #initialize node
    conf = Configuration('config.ini') 

    node = plasma.Grid(conf.Nx, conf.Ny)

    xmin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.dy*conf.Ny*conf.NyMesh

    node.setGridLims(xmin, xmax, ymin, ymax)


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
    injector.inject(node, filler, conf, clip=False) #injecting plasma



    # visualize initial condition
    plotNode(axs[0], node, conf)
    plotXmesh(axs[1], node, conf, 0)
    plotXmesh(axs[2], node, conf, 1)
    saveVisz(-1, node, conf)



    #setup momentum space solver
    vsol = plasma.MomentumLagrangianSolver()
    #intp = ptools.BundleInterpolator2nd()
    #intp = ptools.BundleInterpolator4th()
    intp = ptools.BundleInterpolator4PIC()
    vsol.setInterpolator(intp)


    #setup spatial space solver
    ssol = plasma.SpatialLagrangianSolver2nd()
    #ssol = plasma.SpatialLagrangianSolver4th()
    ssol.setGrid(node)
    


    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    #setup output file
    import h5py
    f = h5py.File("out/run.hdf5", "w")

    grp0 = f.create_group("params")
    grp0.attrs['dx']    = conf.dx
    grp0.attrs['dt']    = conf.interval*conf.dt
    grp = f.create_group("fields")


    #number of samples
    Nsamples    = int(conf.Nt/conf.interval) + 1
    dset = grp.create_dataset("Ex", (conf.Nx*conf.NxMesh, Nsamples), dtype='f')



    #simulation loop
    time = 0.0
    ifile = 0
    for lap in range(0, conf.Nt):

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

        timer.lap("step")

        #I/O
        if (lap % conf.interval == 0):
            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 
            timer.stats("step")


            timer.start("io")

            #save temporarily to file
            save(node, conf, ifile, dset)
            ifile += 1

            plotNode(axs[0], node, conf)
            plotXmesh(axs[1], node, conf, 0) #electrons
            plotXmesh(axs[2], node, conf, 1) #positrons

            plotJ(axs[3], node, conf)
            plotE(axs[4], node, conf)
            saveVisz(lap, node, conf)


            timer.stop("io")
            timer.stats("io")

            timer.start("step") #refresh lap counter (avoids IO profiling)




        time += conf.dt
    
    f.close()
    #node.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
