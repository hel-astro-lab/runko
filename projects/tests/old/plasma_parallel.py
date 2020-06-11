from __future__ import print_function

import numpy as np
#import matplotlib.pyplot as plt
import h5py

import sys, os

import corgi
import pyplasma as plasma


from configSetup import Configuration
import initialize as init

#from visualize import plotNode
#from visualize_amr import plotXmesh
#from visualize import plotJ, plotE, plotDens
#from visualize import saveVisz
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

    mux = 0.0
    muy = 0.0
    muz = 0.0
    delgam = conf.delgam

    #1d filler
    if not( (np.abs(uy) < 0.01) and (np.abs(uz) < 0.01) ):
        return 0.0


    #box advection test
    #dv = (conf.vxmax - conf.vxmin)/(conf.Nvx - 1.0)
    #nn = 1.0/(dv)
    #nn = 1.0

    #if not((x >= 1.0) and (x<=1.01) ):
    #    return 0.0
    #    #if -0.005 < ux < 0.005:
    #    #    return nn
    #    #else:
    #    #    return 0.0



    #speedtest
    #if 0.0<x<0.1:
    #    if 0.98 < ux < 1.02:
    #        return 1.0


    #current test
    #if True: #positive
    #    v = 0.1
    #    dv = (conf.vxmax - conf.vxmin)/(conf.Nvx - 1.0)
    #    if (v-dv/2.0 <= ux <= v+dv/2.0):
    #        return nn
    #    else:
    #        return 0.0
    #else: #negative
    #    if -0.0995 > ux > -0.105:
    #        return nn
    #    else:
    #        return 0.0

    #    #return n0
    #return 0.0


    #electrons
    if ispcs == 0:
        delgam  = conf.delgam * np.abs(conf.mi / conf.me) * conf.temperature_ratio

        # bulk velocities
        mux = conf.gamma_e
        muy = 0.0
        muz = 0.0

        #Lx  = conf.Nx*conf.NxMesh*conf.dx
        #mux_noise += np.sum( conf.beta*np.sin( 2*np.pi*( -modes*x/Lx + random_phase)) )

    #positrons/ions/second species
    if ispcs == 1:
        delgam  = conf.delgam

        # bulk velocities
        mux = conf.gamma_i
        muy = 0.0
        muz = 0.0


    if False:
        alpha = 1.0e-2
        if ispcs < 2:
            n0 = alpha*0.5
        if ispcs >= 2:
            n0 = (1.0-alpha)*0.5

    #plasma frequency scale
    #n0 = 1.0/conf.Nspecies

    #plasma reaction
    omp = conf.cfl*conf.dx
    n0 = (omp**2.0)/conf.Nspecies

    #phase space cell volume normalization
    #dv = (conf.vxmax - conf.vxmin)/(conf.Nvx - 1.0)
    #n0 *= conf.dx/dv 
    #n0 *= 1.0/conf.dx

    #print(n0)
    #n0 = 1.0/conf.Nspecies


    #Brownian noise
    #brownian_noise = 0.01*np.random.standard_normal() 
    #brownian_noise *= delgam

    
    #velocity perturbation
    Lx = conf.Nx*conf.NxMesh*conf.dx
    kmode = conf.modes
    mux_noise = conf.beta*np.cos(2.0*np.pi*kmode*x/Lx) * (Lx/(2.0*np.pi*kmode))
    #mux_noise *= np.sqrt(delgam)/conf.dx/np.sqrt(conf.cfl) #normalize 


    #Classical Maxwellian
    f  = n0*(1.0/(2.0*np.pi*delgam))**(0.5)
    f *= np.exp(-0.5*((ux - mux - mux_noise)**2.0)/(delgam))


    #number density perturbation
    #Lx = conf.Nx*conf.NxMesh*conf.dx
    #kmode = 2.0 #mode
    #f *= 1.0 + conf.beta*np.cos(2.0*np.pi*kmode*x/Lx)


    return f



# Get Yee grid components from grid and save to hdf5 file
def save(n, conf, lap, f5):

    #get E field
    yee = get_yee(n, conf)

    f5['fields/Ex'  ][:,lap] = yee['ex']
    f5['fields/rho' ][:,lap] = yee['rho']
    #f5['fields/ekin'][:,lap] = yee['ekin']
    f5['fields/jx'  ][:,lap] = yee['jx']

    return


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(grid, conf):

    Lx  = conf.Nx*conf.NxMesh*conf.dx
    k = 2.0 #mode

    n0 = 1.0

    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            c = grid.get_tileptr(i,j)
            yee = c.get_yee(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        #get x_i+1/2 (Yee lattice so rho_i)
                        xloc0 = injector.spatialLoc(grid, (i,j), (l,  m,n), conf)
                        xloc1 = injector.spatialLoc(grid, (i,j), (l+1,m,n), conf)

                        #get x_i-1/2 (Yee lattice so rho_i)
                        #xloc0 = injector.spatialLoc(grid, (i,j), (l,  m,n), conf)
                        #xloc1 = injector.spatialLoc(grid, (i,j), (l-1,m,n), conf)

                        xmid = 0.5*(xloc0[0] + xloc1[0])
                        yee.ex[l,m,n] = n0*conf.me*conf.beta*np.sin(2.0*np.pi*k*xmid/Lx)/k

                        #yee.ex[l,m,n] = 1.0e-5


def solvePoisson(ax, grid, conf):
    yee = get_yee(n, conf)

    x   = yee['x']
    rho = yee['rho']



if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    #plt.fig = plt.figure(1, figsize=(8,9))
    #plt.rc('font', family='serif', size=12)
    #plt.rc('xtick')
    #plt.rc('ytick')
    #
    #gs = plt.GridSpec(8, 1)
    #gs.update(hspace = 0.5)
    #
    #axs = []
    #axs.append( plt.subplot(gs[0]) )
    #axs.append( plt.subplot(gs[1]) )
    #axs.append( plt.subplot(gs[2]) )
    #axs.append( plt.subplot(gs[3]) )
    #axs.append( plt.subplot(gs[4]) )
    #axs.append( plt.subplot(gs[5]) )
    #axs.append( plt.subplot(gs[6]) )
    #axs.append( plt.subplot(gs[7]) )


    # Timer for profiling
    timer = Timer(["total", "init", "step", "io"])
    timer.start("total")
    timer.start("init")


    ################################################## 
    #initialize grid
    conf = Configuration('config-landau.ini') 
    #conf = Configuration('config-twostream.ini') 
    #conf = Configuration('config-twostream-fast.ini') 
    #conf = Configuration('config-bump-on-tail.ini') 
    #conf = Configuration('config-twostream-relativistic.ini') 
    #conf = Configuration('config-plasmaosc.ini') 
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

    #insert initial electric field
    #insert_em(grid, conf)

    lap = 0
    plasma.write_yee(grid,      lap)
    plasma.write_analysis(grid, lap)
    plasma.write_mesh(grid,     lap)


    #Initial step backwards for velocity
    for j in range(grid.get_Ny()):
        for i in range(grid.get_Nx()):
            cell = grid.get_tileptr(i,j)
            cell.update_boundaries(grid)
    plasma.initial_step_1d(grid)
    for j in range(grid.get_Ny()):
        for i in range(grid.get_Nx()):
            cell = grid.get_tileptr(i,j)
            cell.cycle()


    # visualize initial condition
    #plotNode(axs[0], grid, conf)
    #plotXmesh(axs[1], grid, conf, 0, "x")
    ##plotXmesh(axs[2], grid, conf, 0, "y")
    #if conf.Nspecies == 2:
    #    plotXmesh(axs[3], grid, conf, 1, "x")
    #    #plotXmesh(axs[4], grid, conf, 1, "y")
    #plotJ(axs[5], grid, conf)
    #plotE(axs[6], grid, conf)
    #plotDens(axs[7], grid, conf)
    #saveVisz(-1, grid, conf)




    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    #setup output file
    f5 = h5py.File(conf.outdir+"/run.hdf5", "w")

    #print(conf.dt)
    #print(conf.cfl/conf.c_omp)

    grp0 = f5.create_group("params")
    #grp0.attrs['c_omp'] = conf.c_omp
    grp0.attrs['dx']    = conf.dx
    #grp0.attrs['dt']    = conf.interval*conf.dt
    grp0.attrs['dt']    = conf.dt
    grp = f5.create_group("fields")


    #number of samples (every step is saved)
    #Nsamples = int(conf.Nt/conf.interval) + 1
    Nsamples = conf.Nt
    dset  = grp.create_dataset("Ex",   (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset2 = grp.create_dataset("rho",  (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    #dset3 = grp.create_dataset("ekin", (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset4 = grp.create_dataset("jx",   (conf.Nx*conf.NxMesh, Nsamples), dtype='f')


    lap = 0
    plasma.write_yee(grid,      lap)
    plasma.write_analysis(grid, lap)
    plasma.write_mesh(grid,     lap)


    #simulation loop
    time  = 0.0
    ifile = 0
    for lap in range(0, conf.Nt):

        #xJEu loop (Umeda a la implicit FTDT)

        #configuration space push
        plasma.step_location(grid)

        #cycle to the new fresh snapshot
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.cycle()

        #current deposition from moving flux
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.deposit_current()

        #update boundaries
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.update_boundaries(grid)

        #momentum step
        plasma.step_velocity_1d(grid)

        #cycle to the new fresh snapshot
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                cell = grid.get_tileptr(i,j)
                cell.cycle()




        ##################################################
        #diagnostics

        #clip every cell
        if conf.clip:
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    cell = grid.get_tileptr(i,j)
                    cell.clip()

        # analyze
        plasma.analyze(grid)


        timer.lap("step")

        #save temporarily to file
        #save(grid, conf, ifile, f5)
        ifile += 1

        #sys.exit()

        #I/O
        if (lap % conf.interval == 0):
            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 
            timer.stats("step")

            timer.start("io")

            plasma.write_yee(grid,      lap)
            plasma.write_analysis(grid, lap)
            plasma.write_mesh(grid,     lap)


            #plotNode(axs[0], grid, conf)

            #plotXmesh(axs[1], grid, conf, 0, "x")
            ##plotXmesh(axs[2], grid, conf, 0, "y")

            #if conf.Nspecies == 2:
            #    plotXmesh(axs[3], grid, conf, 1, "x")
            #    #plotXmesh(axs[4], grid, conf, 1, "y")

            #if conf.Nspecies == 4:
            #    plotXmesh(axs[2], grid, conf, 1, "x")
            #    plotXmesh(axs[3], grid, conf, 2, "x")
            #    plotXmesh(axs[4], grid, conf, 3, "x")

            #plotJ(axs[5], grid, conf)
            #plotE(axs[6], grid, conf)
            #plotDens(axs[7], grid, conf)


            ##solve Poisson
            ##exP = solvePoisson(axs[6], grid, conf)


            #saveVisz(lap, grid, conf)


            timer.stop("io")
            timer.stats("io")

            timer.start("step") #refresh lap counter (avoids IO profiling)

        time += conf.dt
        #time += conf.cfl/conf.c_omp

    
    f5.close()
    #grid.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
