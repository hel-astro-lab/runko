from __future__ import print_function

import numpy as np

import sys, os
import h5py

import pycorgi.twoD as corgi
import pyplasmabox.tools.twoD as pytools
import pyplasmabox.vlv.twoD as pyvlv
import pyplasmabox.pic.twoD as pypic


from configSetup import Configuration
import argparse
import initialize as init
from initialize_pic import loadTiles
from sampling import boosted_maxwellian
from initialize_pic import spatialLoc
from injector_pic import inject
from injector_pic import insert_em



np.random.seed(1)


try:
    import matplotlib.pyplot as plt
    from visualize import plotNode
    from visualize import plotJ, plotE, plotDens
    from visualize import getYee
    from visualize import saveVisz
    
    from visualize import plot2dYee
    from visualize_pic import plot2dParticles


except:
    pass


from timer import Timer

def linear_field(x, y, z):
    return 1.0*x + 2.0*y
    #return 1.0*x + 10.0*y + 100.0*z




# visualize particle content in vx direction
def plotXmesh(ax, n, conf, spcs, vdir):

    ax.clear()

    for i in range(conf.Nx):
        cid = n.id(i,0)
        c = n.getTile(cid)

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

    ax.set_xlim(n.getXmin(), n.getXmax())
    ax.set_ylim(-0.2, 0.2)



def filler(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    #zz = xloc[2] + np.random.rand(1)
    zz = 0.0


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



# Get Yee grid components from node and save to hdf5 file
def save(n, conf, lap, f5):

    #get E field
    yee = getYee(n, conf)
    ex = yee['ex']
    exS = smooth(ex, 10)


    #f5['fields/Ex'  ][:,lap] = yee['ex']
    f5['fields/Ex'  ][:,lap] = exS
    f5['fields/rho' ][:,lap] = yee['rho']
    #f5['fields/ekin'][:,lap] = yee['ekin']
    f5['fields/jx'  ][:,lap] = yee['jx']

    return



if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    try:
        plt.fig = plt.figure(1, figsize=(8,10))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(4, 3)
        gs.update(hspace = 0.5)
        
        axs = []
        for ai in range(12):
            axs.append( plt.subplot(gs[ai]) )
    except:
        pass


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


    #node = plasma.Grid(conf.Nx, conf.Ny)
    node = corgi.Node(conf.Nx, conf.Ny, conf.Nz)

    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh #XXX scaled length
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh

    node.setGridLims(xmin, xmax, ymin, ymax)

    loadTiles(node, conf)


    ################################################## 
    # Path to be created 
    #if node.master:
    if True:
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)



    inject(node, filler, conf) #injecting plasma particles
    #insert_em(node, conf, linear_field)


    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    # visualize initial condition
    try:
        plotNode( axs[0], node, conf)
        plotXmesh(axs[1], node, conf, 0, "x")
        saveVisz(-1, node, conf)
    except:
        pass

    
    #setup output file
    f5 = h5py.File(conf.outdir+"/run.hdf5", "w")
    grp0 = f5.create_group("params")
    grp0.attrs['dx']    = 1.0/conf.c_omp
    #grp0.attrs['dt']    = conf.interval*conf.dt
    grp0.attrs['dt']    = conf.cfl/conf.c_omp
    grp = f5.create_group("fields")


    Nsamples = conf.Nt
    dset  = grp.create_dataset("Ex",   (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset2 = grp.create_dataset("rho",  (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    #dset3 = grp.create_dataset("ekin", (conf.Nx*conf.NxMesh, Nsamples), dtype='f')
    dset4 = grp.create_dataset("jx",   (conf.Nx*conf.NxMesh, Nsamples), dtype='f')







    #TODO:

    #-DONE: field interpolator
    #-DONE: Vau/Boris vel pusher
    #   -position update
    #-DONE:deposit particles (zigzag)
    #-DONE: boundary wrapper
    #-DONE:filtering

    pusher   = pypic.BorisPusher()
    fintp    = pypic.LinearInterpolator()
    comm     = pypic.Communicator()
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
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.updateBoundaries(node)

        #push B half
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.pushHalfB()

        #update boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.updateBoundaries(node)

        #--------------------------------------------------
        # move particles

        #interpolate fields
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                fintp.solve(tile)

        #pusher
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                pusher.solve(tile)

        #--------------------------------------------------
        # advance B half

        #push B half
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.pushHalfB()

        ##update boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.updateBoundaries(node)

        #--------------------------------------------------
        # advance E 

        #push E
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.pushE()

        #--------------------------------------------------

        #deposit current
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                currint.solve(tile)

        #exchange currents
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.exchangeCurrents(node)

        ##################################################
        # particle communication 


        #update particle boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                comm.check_outgoing_particles(tile)

        #copy particles
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                comm.get_incoming_particles(tile, node)

        #delete transferred particles
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                comm.delete_transferred_particles(tile)

        # field communication


        ##################################################

        #filter
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                flt.get_padded_current(tile, node)

                #flt.fft_image_forward()
                #flt.apply_kernel()
                #flt.fft_image_backward()
        
                for fj in range(conf.npasses):
                    flt.direct_convolve_3point()
                flt.set_current(tile)


        ##cycle new and temporary currents
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.cycleCurrent()

        ##################################################
        

        #add current to E
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                tile = node.getTile(i,j)
                tile.depositCurrent()


        ##################################################
        # data reduction and I/O


        timer.lap("step")

        #save temporarily to file
        #save(node, conf, ifile, f5)
        #ifile += 1


        #I/O
        if (lap % conf.interval == 0):
            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 

            timer.stats("step")
            timer.start("io")

            #analyze
            for j in range(node.getNy()):
                for i in range(node.getNx()):
                    tile = node.getTile(i,j)
                    analyzer.analyze2d(tile)


            pyvlv.writeYee(node,      lap, conf.outdir + "/")
            pyvlv.writeAnalysis(node, lap, conf.outdir + "/")
            #pyvlv.writeMesh(node,     lap, conf.outdir + "/")

            #try:
            #    plotNode( axs[0], node, conf)
            #    plotXmesh(axs[1], node, conf, 0, "x")
            #    plotJ(    axs[5], node, conf)
            #    plotE(    axs[6], node, conf)
            #    plotDebug(axs[6], node, conf)
            #    plotDens( axs[7], node, conf)
            #    saveVisz(lap, node, conf)

            #--------------------------------------------------
            #2D plots
            #try:
            #plotNode(axs[0], node, conf)
            #plot2dParticles(axs[1], node, conf, downsample=0.001)
            #plot2dYee(axs[2], node, conf, 'rho')
            #plot2dYee(axs[3], node, conf, 'jx')
            #plot2dYee(axs[4], node, conf, 'jy')
            #plot2dYee(axs[5], node, conf, 'jz')
            #plot2dYee(axs[6], node, conf, 'ex')
            #plot2dYee(axs[7], node, conf, 'ey')
            #plot2dYee(axs[8], node, conf, 'ez')
            #plot2dYee(axs[9], node, conf, 'bx')
            #plot2dYee(axs[10],node, conf, 'by')
            #plot2dYee(axs[11],node, conf, 'bz')


            #saveVisz(lap, node, conf)

	    #except:
	    #    pass


            timer.stop("io")
            timer.stats("io")
            timer.start("step") #refresh lap counter (avoids IO profiling)

        time += conf.cfl/conf.c_omp
    #end of loop

    #node.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
