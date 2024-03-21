from __future__ import print_function

import numpy as np

import sys, os

import pycorgi
import pyrunko.vlv.oneD as plasma


from configSetup import Configuration
import initialize as init

# for on-the-fly visualization
try:
    import matplotlib.pyplot as plt
    from visualize import plotNode
    from visualize_amr import plotXmesh
    from visualize import plotJ, plotE, plotDens
    from visualize import saveVisz
except:
    pass

from visualize import get_grids
from visualize import get_analysis
import injector

from timer import Timer

import argparse


# preclipper to accelerate injector
def preclip(xloc, uloc, ispcs, conf):
    if ispcs == 0:
        delgam  = conf.delgam * np.abs(conf.mi / conf.me) * conf.temperature_ratio
        mux = conf.ub_e
    elif ispcs == 1:
        delgam  = conf.delgam
        mux = conf.ub_i
    
    if ((uloc[0] - mux)**2.0)/delgam < 5.0:
        return False
    else:
        return True


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
    #if not( (np.abs(uy) < 0.01) and (np.abs(uz) < 0.01) ):
    #    return 0.0

    #box advection test
    #dv = (conf.vxmax - conf.vxmin)/(conf.Nvx - 1.0)
    ##nn = 1.0/(dv)
    #nn = 1.0

    #if ((x >= 1.0) and (x<=1.02) ):
    #    return 1.0
    #    #if  0.016 < ux < 0.004:
    #    #    return nn
    #    #else:
    #    #    return 0.0
    #else:
    #    return 0.0

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
        #print("ub_e:", conf.ub_e)
        mux = conf.ub_e
        muy = 0.0
        muz = 0.0

        #Lx  = conf.Nx*conf.NxMesh*conf.dx
        #mux_noise += np.sum( conf.beta*np.sin( 2*np.pi*( -modes*x/Lx + random_phase)) )

    #positrons/ions/second species
    if ispcs == 1:
        delgam  = conf.delgam

        # bulk velocities
        #print("ub_i:", conf.ub_i)
        mux = conf.ub_i
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

    #phase space tile volume normalization
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
    yee = get_grids(n, conf)
    analysis1 = get_analysis(n, conf, 0)
    analysis2 = get_analysis(n, conf, 1)

    f5['fields/Ex'  ][:,lap] = yee['ex']
    f5['fields/rho' ][:,lap] = yee['rho']
    f5['fields/jx'  ][:,lap] = yee['jx']
    f5['fields/ekin'][:,lap] = analysis1['edens'] + analysis2['edens']

    return


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(grid, conf):

    Lx  = conf.Nx*conf.NxMesh*conf.dx
    k = 2.0 #mode

    n0 = 1.0

    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            c = grid.get_tile(i,j)
            yee = c.get_grids(0)

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
    yee = get_grids(n, conf)

    x   = yee['x']
    rho = yee['rho']



if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    try:
        plt.fig = plt.figure(1, figsize=(8,9))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(8, 1)
        gs.update(hspace = 0.5)
        
        axs = []
        for ai in range(8):
            axs.append( plt.subplot(gs[ai]) )
    except:
        pass

    # Timer for profiling
    timer = Timer(
            ["total", "init", "step", "io"],
            ["cycle1", "cycle2", "loc", "vel", "cur-dep", "bounds", "clip", "analyze"]
            )
    timer.start("total")
    timer.start("init")


    # parse command line arguments
    parser = argparse.ArgumentParser(description='Simple Vlasov-Maxwell simulations')
    parser.add_argument('--conf', dest='conf_filename', default=None,
                       help='Name of the configuration file (default: None)')
    args = parser.parse_args()
    if args.conf_filename == None:
        conf = Configuration('config-landau.ini') 
        #conf = Configuration('config-twostream.ini') 
        #conf = Configuration('config-twostream-fast.ini') 
        #conf = Configuration('config-bump-on-tail.ini') 
        #conf = Configuration('config-twostream-relativistic.ini') 
        #conf = Configuration('config-plasmaosc.ini') 
        #conf = Configuration('config-dispersion.ini') 
    else:
        print("Reading configuration setup from ", args.conf_filename)
        conf = Configuration(args.conf_filename)


    ################################################## 
    #initialize grid
    grid = pycorgi.oneD.Grid(conf.Nx, conf.Ny, conf.Nz)

    xmin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.dy*conf.Ny*conf.NyMesh

    grid.set_grid_lims(xmin, xmax, ymin, ymax)


    #grid.initMpi()
    #loadMpiXStrides(grid)

    init.loadTiles(grid, conf)


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
    random_phase = np.random.rand(len(modes))

    # restart
    if conf.laprestart > 0:
        lap = conf.laprestart + 1
        injector.inject(grid, injector.empty_filler, conf, empty=True) #injecting plasma

        plasma.read_grids( grid, conf.laprestart, conf.outdir)
        plasma.read_mesh(grid, conf.laprestart, conf.outdir)
    else:
        lap = 0
        injector.inject(grid, filler, conf, preclip) #injecting plasma

        #insert initial electric field
        #insert_em(grid, conf)

        #Initial step backwards for velocity
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.update_boundaries(grid)
        plasma.initial_step(grid)
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.cycle()


    # visualize initial condition
    try:
        plotNode(axs[0], grid, conf)
        plotXmesh(axs[1], grid, conf, 0, "x")
        #plotXmesh(axs[2], grid, conf, 0, "y")
        if conf.Nspecies == 2:
            plotXmesh(axs[3], grid, conf, 1, "x")
            #plotXmesh(axs[4], grid, conf, 1, "y")
        plotJ(axs[5], grid, conf)
        plotE(axs[6], grid, conf)
        plotDens(axs[7], grid, conf)
        saveVisz(-1, grid, conf)
    except:
        pass




    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    #setup output file
    import h5py
    f5 = h5py.File(conf.outdir+"/run.hdf5", "w")

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
    dset3 = grp.create_dataset("ekin", (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset4 = grp.create_dataset("jx",   (conf.Nx*conf.NxMesh, Nsamples), dtype='f')




    #simulation loop
    time  = conf.dt*lap #start time
    ifile = 0
    for lap in range(lap, conf.Nt):

        #xJEu loop (Umeda a la implicit FTDT)

        #configuration space push
        timer.start_comp("loc")
        plasma.step_location(grid)
        timer.stop_comp("loc")

        #cycle to the new fresh snapshot
        timer.start_comp("cycle1")
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.cycle()
        timer.stop_comp("cycle1")


        #current deposition from moving flux
        timer.start_comp("cur-dep")
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.deposit_current()
        timer.stop_comp("cur-dep")

        #update boundaries
        timer.start_comp("bounds")
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.update_boundaries(grid)
        timer.stop_comp("bounds")

        #momentum step
        timer.start_comp("vel")
        plasma.step_velocity(grid)
        timer.stop_comp("vel")


        #cycle to the new fresh snapshot
        timer.start_comp("cycle2")
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                tile = grid.get_tile(i,j)
                tile.cycle()
        timer.stop_comp("cycle2")




        ##################################################
        #diagnostics

        #clip every tile
        timer.start_comp("clip")
        if conf.clip:
            for j in range(grid.get_Ny()):
                for i in range(grid.get_Nx()):
                    tile = grid.get_tile(i,j)
                    #tile.clip_neighbors()
                    tile.clip()
        timer.stop_comp("clip")




        timer.lap("step")

        #save temporarily to file
        save(grid, conf, ifile, f5)
        ifile += 1

        #sys.exit()

        # analyze (this is done for every step because run.hdf5 is updated such a way)
        timer.start_comp("analyze")
        plasma.analyze(grid)
        timer.stop_comp("analyze")

        #I/O
        if (lap % conf.interval == 0):

            # make sure hdf5 does not corrupt by regularly flushing it
            f5.flush()

            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 
            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()

            timer.start("io")

            plasma.write_grids(grid,      lap, conf.outdir + "/")
            plasma.write_analysis(grid, lap, conf.outdir + "/")

            if (lap % conf.restart == 0):
                plasma.write_mesh(grid,     lap, conf.outdir + "/")

            try:
                plotNode(axs[0], grid, conf)

                plotXmesh(axs[1], grid, conf, 0, "x")
                #plotXmesh(axs[2], grid, conf, 0, "y")

                if conf.Nspecies == 2:
                    plotXmesh(axs[3], grid, conf, 1, "x")
                    #plotXmesh(axs[4], grid, conf, 1, "y")

                if conf.Nspecies == 4:
                    plotXmesh(axs[2], grid, conf, 1, "x")
                    plotXmesh(axs[3], grid, conf, 2, "x")
                    plotXmesh(axs[4], grid, conf, 3, "x")

                plotJ(axs[5], grid, conf)
                plotE(axs[6], grid, conf)
                plotDens(axs[7], grid, conf)

                saveVisz(lap, grid, conf)
            except:
                pass

            timer.stop("io")
            timer.stats("io")

            timer.start("step") #refresh lap counter (avoids IO profiling)

        time += conf.dt
        #time += conf.cfl/conf.c_omp

    
    f5.close()
    #grid.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
