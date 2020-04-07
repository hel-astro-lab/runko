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
from visualize import get_yee

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
    if not( (np.abs(uy) < 0.01) and (np.abs(uz) < 0.01) ):
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
    f  = 0.5*(1.0/(2.0*np.pi*delgam))**(0.5)

    #f *= np.exp(-0.5*((ux - mux - mux_noise)**2)/(delgam + delgam_noise) + brownian_noise)
    #f *= np.exp(-0.5*( (ux - mux - mux_noise)**2 + (uy - muy)**2 + (uz - muz)**2)/(delgam))
    f *= np.exp(-0.5*( (ux - mux - mux_noise)**2)/(delgam))


    #number density oscillations
    #Lx  = conf.Nx*conf.NxMesh*conf.dx
    #k  = 2.0*np.pi/0.568
    #f *= 1.0 + conf.beta*np.cos(k*x/conf.dx)


    return f


# Get Yee grid components from grid and save to hdf5 file
def save(n, conf, lap, f5):

    #get E field
    yee = get_yee(n, conf)

    f5['fields/Ex'  ][:,lap] = yee['ex']
    f5['fields/rho' ][:,lap] = yee['rho']
    f5['fields/ekin'][:,lap] = yee['ekin']

    return





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
    #initialize grid
    conf = Configuration('config-plasmaosc.ini') 
    #conf = Configuration('config-dispersion.ini') 

    grid = plasma.Grid(conf.Nx, conf.Ny)

    xmin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.dy*conf.Ny*conf.NyMesh

    grid.set_grid_lims(xmin, xmax, ymin, ymax)


    #grid.initMpi()
    #loadMpiXStrides(grid)

    init.loadCells(grid, conf)


    ################################################## 
    # Path to be created 
    #if grid.master:
    if True:
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)


    ################################################## 
    # initialize
    Nx           = conf.Nx*conf.NxMesh
    modes        = np.arange(Nx) 
    #modes        = np.array([2])
    random_phase = np.random.rand(len(modes))

    injector.inject(grid, filler, conf) #injecting plasma


    # visualize initial condition
    plotNode(axs[0], grid, conf)
    plotXmesh(axs[1], grid, conf, 0, "x")
    plotXmesh(axs[2], grid, conf, 0, "y")
    plotXmesh(axs[3], grid, conf, 1, "x")
    plotXmesh(axs[4], grid, conf, 1, "y")
    plotJ(axs[5], grid, conf)
    plotE(axs[6], grid, conf)
    plotDens(axs[7], grid, conf)
    saveVisz(-1, grid, conf)




    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    #setup output file
    import h5py
    f5 = h5py.File("out/run.hdf5", "w")

    grp0 = f5.create_group("params")
    grp0.attrs['dx']    = conf.dx
    grp0.attrs['dt']    = conf.dt
    grp = f5.create_group("fields")

    #number of samples (every step is saved)
    #Nsamples = int(conf.Nt/conf.interval) + 1
    Nsamples = conf.Nt
    dset  = grp.create_dataset("Ex",   (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset2 = grp.create_dataset("rho",  (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset3 = grp.create_dataset("ekin", (conf.Nx*conf.NxMesh, Nsamples), dtype='f')



    ##################################################

    #simulation loop
    time  = 0.0
    ifile = 0
    for lap in range(0, conf.Nt):

        #B field half update

        ##move vlasov fluid

        #update boundaries
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.update_boundaries(grid)

        #momentum step
        #for j in range(grid.get_Ny()):
        #    for i in range(grid.get_Nx()):
        #        cell = grid.get_tileptr(i,j)
        #        vsol.solve(cell)
        plasma.step_velocity_1d(grid)


        #cycle to the new fresh snapshot
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.cycle()

        #spatial step
        #for j in range(grid.get_Ny()):
        #    for i in range(grid.get_Nx()):
        #        cell = grid.get_tileptr(i,j)
        #        ssol.solve(cell, grid)
        plasma.step_location(grid)

        #cycle to the new fresh snapshot
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.cycle()

        #B field second half update

        #E field (Ampere's law)
        #for cid in grid.getCellIds():
        #    c = grid.get_tileptr( cid )
        #    c.push_e()


        #current deposition from moving flux
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.deposit_current()

        #clip every cell
        if conf.clip:
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    cell = grid.get_tileptr(i,j)
                    cell.clip()

        # analyze
        plasma.analyze(grid)

        #save temporarily to file
        save(grid, conf, ifile, f5)
        ifile += 1




        timer.lap("step")

        #I/O
        if (lap % conf.interval == 0):
            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 
            timer.stats("step")


            timer.start("io")


            plotNode(axs[0], grid, conf)

            plotXmesh(axs[1], grid, conf, 0, "x")
            plotXmesh(axs[2], grid, conf, 0, "y")
            plotXmesh(axs[3], grid, conf, 1, "x")
            plotXmesh(axs[4], grid, conf, 1, "y")

            plotJ(axs[5], grid, conf)
            plotE(axs[6], grid, conf)
            plotDens(axs[7], grid, conf)

            saveVisz(lap, grid, conf)


            timer.stop("io")
            timer.stats("io")

            timer.start("step") #refresh lap counter (avoids IO profiling)


        time += conf.dt
    
    f5.close()
    #grid.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
