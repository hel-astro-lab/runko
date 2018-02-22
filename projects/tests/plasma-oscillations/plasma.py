from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os

import corgi
import pyplasma as plasma


from configSetup import Configuration
import initialize as init

from visualize import plotNode
from visualize_amr import plotXmesh
from visualize import plotJ, plotE, plotDens
from visualize import saveVisz
from visualize import getYee

import injector

from timer import Timer


# Generic function to fill the velocity mesh
#
# Maxwellian plasma with Brownian noise
# where delgam = kT/m_i c^2
#
def filler(xloc, uloc, ispcs, conf):

    mux_noise      = 0.0
    delgam_noise   = 0.0
    brownian_noise = 0.0

    x = xloc[0]
    y = xloc[1]
    z = xloc[2] 

    ux = uloc[0]
    uy = uloc[1]
    uz = uloc[2] 

    #print("uy={} uz={}".format(uy,uz))

    #1d filler
    if not( (uy == 0.0) and (uz == 0.0) ):
        return 0.0

    #electrons
    if ispcs == 0:
        delgam  = conf.delgam * np.abs(conf.mi / conf.me) * conf.temperature_ratio

        # bulk velocities
        mux = conf.gamma_e
        muy = 0.0
        muz = 0.0

        Lx  = conf.Nx*conf.NxMesh*conf.dx
        mux_noise += np.sum( conf.beta*np.sin( 2*np.pi*( -modes*x/Lx + random_phase)) )


    #ions/positrons
    elif ispcs == 1:
        delgam  = conf.delgam

        # bulk velocities
        mux = conf.gamma_i
        muy = 0.0
        muz = 0.0



    #Brownian noise
    #brownian_noise = 0.01*np.random.standard_normal() 
    #brownian_noise *= delgam


    #Classical Maxwellian distribution
    #f  = (1.0/(2.0*np.pi*delgam))**(3.0/2.0)
    f  = (2.0)*(1.0/(2.0*np.pi*delgam))**(1.0/2.0)

    #f *= np.exp(-0.5*((ux - mux - mux_noise)**2)/(delgam + delgam_noise) + brownian_noise)
    f *= np.exp(-0.5*( (ux - mux - mux_noise)**2 + (uy - muy)**2 + (uz - muz)**2)/(delgam))


    return f



# Get Yee grid components from node and save to hdf5 file
def save(n, conf, lap, dset):

    #get E field
    yee = getYee(n, conf)

    dset[:, lap] = yee['ex']


def insert_em(node, conf):

    for i in range(node.getNx()):
        for j in range(node.getNy()):
            c = node.getCellPtr(i,j)
            yee = c.getYee(0)

            for q in range(conf.NxMesh):
                for k in range(conf.NyMesh):
                    if q == 32:
                        yee.ex[q,k,0] =  0.5

                    #yee.ex[q,k,0] =  0.0
                    #yee.ey[q,k,0] =  0.0
                    #yee.ez[q,k,0] =  0.0

                    #yee.bx[q,k,0] =  0.0
                    #yee.by[q,k,0] =  0.0
                    #yee.bz[q,k,0] =  0.0



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
    axs.append( plt.subplot(gs[3]) )
    axs.append( plt.subplot(gs[4]) )
    axs.append( plt.subplot(gs[5]) )
    axs.append( plt.subplot(gs[6]) )
    axs.append( plt.subplot(gs[7]) )


    # Timer for profiling
    timer = Timer(["total", "init", "step", "io"])
    timer.start("total")
    timer.start("init")


    ################################################## 
    #initialize node
    conf = Configuration('config-plasmaosc.ini') 
    #conf = Configuration('config-dispersion.ini') 
    #conf = Configuration('config-twostream.ini') 

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
    Nx           = conf.Nx*conf.NxMesh
    #modes        = np.arange(Nx) 
    modes        = np.array([4])
    random_phase = np.random.rand(len(modes))

    injector.inject(node, filler, conf) #injecting plasma


    # visualize initial condition
    plotNode(axs[0], node, conf)
    plotXmesh(axs[1], node, conf, 0, "x")
    plotXmesh(axs[2], node, conf, 0, "y")
    plotXmesh(axs[3], node, conf, 1, "x")
    plotXmesh(axs[4], node, conf, 1, "y")
    plotJ(axs[5], node, conf)
    plotE(axs[6], node, conf)
    plotDens(axs[7], node, conf)
    saveVisz(-1, node, conf)



    #setup momentum space solver
    vsol = plasma.AmrMomentumLagrangianSolver()

    #setup spatial space solver
    ssol = plasma.AmrSpatialLagrangianSolver()


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


    #XXX DEBUG
    #insert_em(node, conf)


    #simulation loop
    time  = 0.0
    ifile = 0
    for lap in range(0, conf.Nt):

        #B field half update

        ##move vlasov fluid

        #update boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.updateBoundaries(node)

        #momentum step
        #for j in range(node.getNy()):
        #    for i in range(node.getNx()):
        #        cell = node.getCellPtr(i,j)
        #        vsol.solve(cell)
        plasma.stepVelocity(node)

        #cycle to the new fresh snapshot
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.cycle()

        #spatial step
        #for j in range(node.getNy()):
        #    for i in range(node.getNx()):
        #        cell = node.getCellPtr(i,j)
        #        ssol.solve(cell, node)
        plasma.stepLocation(node)

        #cycle to the new fresh snapshot
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.cycle()

        #B field second half update

        #E field (Ampere's law)
        #for cid in node.getCellIds():
        #    c = node.getCellPtr( cid )
        #    c.pushE()


        #current deposition from moving flux
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.depositCurrent()

        #clip every cell
        if conf.clip:
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    cell = node.getCellPtr(i,j)
                    cell.clip()

        # analyze
        plasma.analyze(node)


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

            plotXmesh(axs[1], node, conf, 0, "x")
            plotXmesh(axs[2], node, conf, 0, "y")
            plotXmesh(axs[3], node, conf, 1, "x")
            plotXmesh(axs[4], node, conf, 1, "y")

            plotJ(axs[5], node, conf)
            plotE(axs[6], node, conf)
            plotDens(axs[7], node, conf)

            saveVisz(lap, node, conf)


            timer.stop("io")
            timer.stats("io")

            timer.start("step") #refresh lap counter (avoids IO profiling)


        time += conf.dt
    
    f.close()
    #node.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
