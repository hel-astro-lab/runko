from __future__ import print_function
from mpi4py import MPI

import numpy as np

import sys, os
import h5py

import pycorgi.twoD as corgi
import pyrunko.tools.twoD as pytools
import pyrunko.vlv.twoD as pyvlv
import pyrunko.pic.twoD as pypic


from configSetup import Configuration
import argparse
import initialize as init
from initialize_pic import loadTiles
from initialize_pic import initialize_virtuals
from sampling import boosted_maxwellian
from initialize_pic import spatialLoc
from injector_pic import inject
from injector_pic import insert_em



np.random.seed(1)


#try:
import matplotlib.pyplot as plt
from visualize import plotNode
from visualize import plotJ, plotE, plotDens
from visualize import get_grids
from visualize import saveVisz

from visualize import getYee2D
from visualize import plot2dYee
from visualize_pic import plot2dParticles

#except:
#pass

from timer import Timer


#debug = True
debug = False


def linear_field(x, y, z):
    return 1.0*x + 2.0*y
    #return 1.0*x + 10.0*y + 100.0*z




# visualize particle content in vx direction
def plotXmesh(ax, n, conf, spcs, vdir):

    ax.clear()

    for i in range(conf.Nx):
        cid = n.id(i,0)
        c = n.get_tile(cid)

        container = c.get_container(spcs)

        x = container.loc(0)
        y = container.loc(1)
        z = container.loc(2)

        ux = container.vel(0)
        uy = container.vel(1)
        uz = container.vel(2)

        ax.plot(x, ux, "k.", markersize=2, alpha=0.8)
        
    if vdir == "x":
        if spcs == 0:
            ax.set_ylabel(r'$v_{x,e}$')
        if spcs == 1:
            ax.set_ylabel(r'$v_{x,p}$')

    ax.minorticks_on()

    ax.set_xlim(n.get_xmin(), n.get_xmax())
    ax.set_ylim(-0.2, 0.2)



def filler(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    #zz = xloc[2] + np.random.rand(1)
    zz = 0.5

    #electrons
    if ispcs == 0:
        delgam  = conf.delgam * np.abs(conf.mi / conf.me) * conf.temperature_ratio

        gamma = abs(conf.gamma_e) # bulk velocities
        direction = np.sign(conf.gamma_e)
        #Lx  = conf.Nx*conf.NxMesh*conf.dx
        #mux_noise += np.sum( conf.beta*np.sin( 2*np.pi*( -modes*x/Lx + random_phase)) )

    #positrons/ions/second species
    if ispcs == 1:
        delgam  = conf.delgam

        gamma = abs(conf.gamma_i) # bulk velocities
        direction = np.sign(conf.gamma_i)

    ux, uy, uz, uu = boosted_maxwellian(delgam, gamma, direction=direction, dims=2)

    #print("injecting into {} ({})".format(xx, xloc[0]))

    #perturb with wave
    #Lx = conf.Nx*conf.NxMesh
    #kmode = conf.modes
    #mux_noise = conf.beta*np.cos(2.0*np.pi*kmode*xx/Lx) * (Lx/(2.0*np.pi*kmode))
    #ux += vth*mux_noise

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0



# Get Yee grid components from grid and save to hdf5 file
def save(n, conf, lap, f5):

    #get E field
    yee = get_grids(n, conf)
    ex = yee['ex']
    exS = smooth(ex, 10)


    #f5['fields/Ex'  ][:,lap] = yee['ex']
    f5['fields/Ex'  ][:,lap] = exS
    f5['fields/rho' ][:,lap] = yee['rho']
    #f5['fields/ekin'][:,lap] = yee['ekin']
    f5['fields/jx'  ][:,lap] = yee['jx']

    return


def debug_print(n, msg):
    if debug:
        print("{}: {}".format(n.rank(), msg))


if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    #try:
    plt.fig = plt.figure(1, figsize=(8,10))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(4, 3)
    gs.update(hspace = 0.5)
    
    axs = []
    for ai in range(12):
        axs.append( plt.subplot(gs[ai]) )
    #except:
    #    print()
        #pass


    # Timer for profiling
    timer = Timer(["total", "init", "step", "io"])
    timer.start("total")
    timer.start("init")

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Simple PIC-Maxwell simulations')
    parser.add_argument('--conf', dest='conf_filename', default=None,
                       help='Name of the configuration file (default: None)')
    args = parser.parse_args()
    if args.conf_filename == None:
        conf = Configuration('config-landau.ini') 
    else:
        print("Reading configuration setup from ", args.conf_filename)
        conf = Configuration(args.conf_filename)


    #grid = plasma.Grid(conf.Nx, conf.Ny)
    grid = corgi.Grid(conf.Nx, conf.Ny, conf.Nz)

    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh #XXX scaled length
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh

    grid.set_grid_lims(xmin, xmax, ymin, ymax)


    #init.loadMpiRandomly(grid)
    init.loadMpiXStrides(grid)

    loadTiles(grid, conf)


    ################################################## 
    # Path to be created 
    #if grid.master:
    if True:
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)



    np.random.seed(1)
    inject(grid, filler, conf) #injecting plasma particles
    #insert_em(grid, conf, linear_field)

    #static setup; communicate neighbor info once
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    initialize_virtuals(grid, conf)


    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    # visualize initial condition
    #try:
    plotNode( axs[0], grid, conf)
    #plotXmesh(axs[1], grid, conf, 0, "x")
    saveVisz(-1, grid, conf)
    #except:
    #    print()
    #    pass


    Nsamples = conf.Nt

    fldprop  = pyfld.FDTD2()
    pusher   = pypic.BorisPusher()
    fintp    = pypic.LinearInterpolator()
    currint  = pypic.ZigZag()
    analyzer = pypic.Analyzator()
    flt     =  pytools.Filter(conf.NxMesh, conf.NyMesh)
    flt.init_gaussian_kernel(2.0, 2.0)

    #simulation loop
    time  = 0.0
    ifile = 0
    for lap in range(0, conf.Nt):

        # Tristan loop
        #advance B half
        #move particles
        #advance B half
        #advance E
        #reset current (j=0)
        #deposit
        #exchange particles
        #exchange current
        #filter
        #add current
        #inject particles from other processors
        #enlarge domain
        #reorder particles
        #pause simulation if pause file exists


        #--------------------------------------------------
        # advance Half B

        #update boundaries
        debug_print(grid,"update_boundaries 0")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            tile.update_boundaries(grid)
        #FIXME: update also virtuals (for push_b)

        #push B half
        debug_print(grid,"push_half_b 1")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            fldprop.push_half_b(tile)
        #FIXME: push also virtuals to get correct boundaries for locals

        #update boundaries
        debug_print(grid,"update_boundaries 1")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            tile.update_boundaries(grid)
        #FIXME: update also virtuals (for second push_b)

        #--------------------------------------------------
        # move particles (only locals tiles)

        #interpolate fields (can move to next asap)
        debug_print(grid,"interpolate fields")
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            fintp.solve(tile)

        #pusher 
        debug_print(grid,"push particles")
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            pusher.solve(tile)

        #--------------------------------------------------
        # advance B half

        #push B half
        debug_print(grid,"push_half_b 2")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            fldprop.push_half_b(tile)
        #FIXME: push also virtuals

        #update boundaries
        debug_print(grid,"update_boundaries 2")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            tile.update_boundaries(grid)
        #FIXME: update virtuals

        #--------------------------------------------------
        # advance E 

        #push E
        debug_print(grid,"push_e")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            fldprop.push_e(tile)
        #FIXME: push also virtuals

        #--------------------------------------------------

        #current calculation
        debug_print(grid,"current computation")
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            currint.solve(tile)
        
        # current solving is also taking place in nbor ranks
        # that is why we update virtuals here with MPI
        #
        # This is the most expensive task so we do not double it 
        # here.

        #mpi send currents
        debug_print(grid,"send 0")
        grid.send_data(0) #(indepdendent)
        debug_print(grid,"recv 0")
        grid.recv_data(0) #(indepdendent)
        debug_print(grid,"wait 0")
        grid.wait_data(0) #(indepdendent)

        #exchange currents
        debug_print(grid,"exchange_currents")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            tile.exchange_currents(grid)
        #FIXME: exchange also virtuals (to get inner boundaries right)


        ##################################################
        # particle communication (only local/boundary tiles)

        #local particle exchange (independent)
        debug_print(grid,"check_outgoing_particles")
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            tile.check_outgoing_particles()

        # global mpi exchange (independent)
        debug_print(grid,"pack_outgoing_particles")
        for cid in grid.get_boundary_tiles():
            tile = grid.get_tile(cid)
            tile.pack_outgoing_particles()

        # MPI global exchange
        # transfer primary and extra data
        debug_print(grid, "send 1")
        grid.send_data(1) #(indepdendent)
        debug_print(grid, "send 2")
        grid.send_data(2) #(indepdendent)

        debug_print(grid, "recv 1")
        grid.recv_data(1) #(indepdendent)

        debug_print(grid,"wait 1")
        grid.wait_data(1) #(indepdendent)

        debug_print(grid, "recv 2")
        grid.recv_data(2) #(indepdendent)

        debug_print(grid,"wait 2")
        grid.wait_data(2) #(indepdendent)


        # global unpacking (independent)
        debug_print(grid,"check_outgoing_particles")
        for cid in grid.get_virtual_tiles(): 
            tile = grid.get_tile(cid)
            tile.unpack_incoming_particles()
            tile.check_outgoing_particles()

        # transfer local + global
        debug_print(grid,"get_incoming_particles")
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            tile.get_incoming_particles(grid)

        # delete local transferred particles
        debug_print(grid,"delete_transferred_particles")
        for cid in grid.get_local_tiles():
            tile = grid.get_tile(cid)
            tile.delete_transferred_particles()

        debug_print(grid,"delete all virtual particles")
        for cid in grid.get_virtual_tiles(): 
            tile = grid.get_tile(cid)
            tile.delete_all_particles()

        ##################################################

        #filter
        #debug_print(grid,"filter")
        #for cid in grid.get_tile_ids():
        #    tile = grid.get_tile(cid)
        #    flt.get_padded_current(tile, grid)

        #    #flt.fft_image_forward()
        #    #flt.apply_kernel()
        #    #flt.fft_image_backward()
        #
        #    for fj in range(conf.npasses):
        #        flt.direct_convolve_3point()
        #    flt.set_current(tile)
        # FIXME: filter also virtuals
        # FIXME: or mpi communicate filtered currents

        ##cycle new and temporary currents (only if filtering)
        #debug_print(grid,"cycle currents")
        #for cid in grid.get_tile_ids():
        #    tile = grid.get_tile(cid)
        #    tile.cycle_current()
        # FIXME: cycle also virtuals

        ##################################################

        #add current to E
        debug_print(grid,"add J to E")
        for cid in grid.get_tile_ids():
            tile = grid.get_tile(cid)
            tile.deposit_current()
        #FIXME: deposit also virtuals

        ##################################################
        # data reduction and I/O

        timer.lap("step")

        #save temporarily to file
        #save(grid, conf, ifile, f5)
        #ifile += 1

        #I/O
        if (lap % conf.interval == 0):
            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 

            timer.stats("step")
            timer.start("io")

            #analyze (independent)
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                analyzer.analyze2d(tile)

            pyvlv.write_grids(grid,      lap, conf.outdir + "/")
            pyvlv.write_analysis(grid, lap, conf.outdir + "/")
            #pyvlv.write_mesh(grid,     lap, conf.outdir + "/")

            #try:
            #    plotNode( axs[0], grid, conf)
            #    plotXmesh(axs[1], grid, conf, 0, "x")
            #    plotJ(    axs[5], grid, conf)
            #    plotE(    axs[6], grid, conf)
            #    plotDebug(axs[6], grid, conf)
            #    plotDens( axs[7], grid, conf)
            #    saveVisz(lap, grid, conf)

            #--------------------------------------------------
            #2D plots
            #try:
            plotNode(axs[0], grid, conf)
            plot2dParticles(axs[1], grid, conf, downsample=0.001)

            yee = getYee2D(grid, conf)
            plot2dYee(axs[2],  yee, grid, conf, 'rho')
            plot2dYee(axs[3],  yee, grid, conf, 'jx')
            plot2dYee(axs[4],  yee, grid, conf, 'jy')
            plot2dYee(axs[5],  yee, grid, conf, 'jz')
            plot2dYee(axs[6],  yee, grid, conf, 'ex')
            plot2dYee(axs[7],  yee, grid, conf, 'ey')
            plot2dYee(axs[8],  yee, grid, conf, 'ez')
            plot2dYee(axs[9],  yee, grid, conf, 'bx')
            plot2dYee(axs[10], yee, grid, conf, 'by')
            plot2dYee(axs[11], yee, grid, conf, 'bz')
            saveVisz(lap, grid, conf)
            #except:
            #    print()
            #    pass

            timer.stop("io")
            timer.stats("io")
            timer.start("step") #refresh lap counter (avoids IO profiling)

        time += conf.cfl/conf.c_omp
    #end of loop


    timer.stop("total")
    timer.stats("total")
