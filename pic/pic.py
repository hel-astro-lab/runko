from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os
import h5py

import corgi
import pyplasma as plasma
import pypic 


from configSetup import Configuration
import argparse
import initialize as init

#from injector import spatialLoc

from visualize import plotNode
from visualize import plotJ, plotE, plotDens
from visualize import getYee
from visualize import saveVisz


from timer import Timer



def spatialLoc(node, Ncoords, Mcoords, conf):

    #node coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    l, m, n = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #grid spacing
    xmin = node.getXmin()
    ymin = node.getYmin()

    #dx = conf.dx
    dx = 1.0 #XXX scaled length
    dy = conf.dy
    dz = conf.dz


    #calculate coordinate extent
    x = xmin + i*NxMesh*dx + l*dx
    y = ymin + j*NyMesh*dy + m*dy
    z = 0.0                + n*dz

    return [x, y, z]




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
                #c.dt  = conf.dt
                #c.dx  = conf.dx
                c.dx  = 1.0
                c.cfl = conf.cfl

                #normalization factors
                omp = conf.cfl/conf.c_omp #plasma reaction
                gamma0 = 1.0      #relativistic dilatation
                #betaN = np.sqrt(1.0 - 1.0/gamma0**2.0)
                qe = (gamma0*omp**2.0)/((conf.ppc*0.5)*(1.0 + np.abs(conf.me/conf.mi)) )
                
                #print("normalization factor: {}".format(qe))

                c.container.qe = qe
                c.container.qi = qe


                # use same scale for Maxwell solver
                #c.yeeDt = conf.dt
                #c.yeeDx = conf.dx
                c.yeeDx = 1.0

                #set bounding box of the tile
                mins = spatialLoc(n, [i,j], [0,0,0], conf)
                maxs = spatialLoc(n, [i,j], [conf.NxMesh, conf.NyMesh, conf.NzMesh], conf)
                c.set_tile_mins(mins)
                c.set_tile_maxs(maxs)

                #add it to the node
                n.addCell(c) 




#inject plasma into cells
def inject(node, ffunc, conf):

    #loop over all *local* cells
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                #print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = node.cellId(i,j)
                c      = node.getCellPtr(cid) #get cell ptr

                #reserve memory for particles
                Nprtcls = conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc
                c.container.reserve(Nprtcls, 3)
                c.container.resizeEM(Nprtcls, 3)

                #if not(i == 1 and j == 1):
                #    continue

                for ispcs in range(conf.Nspecies):

                    #set q/m
                    #if ispcs == 0:
                    #    qm = conf.me
                    #elif ispcs == 1:
                    #    qm = conf.mi

                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))
                                xloc = spatialLoc(node, (i,j), (l,m,n), conf)
                                #xloc[0] += 0.5

                                for ip in range(conf.ppc):
                                    #xloc = [0.0, 0.0, 0.0]
                                    #uloc = [0.1, 0.0, 0.0]
                                    x0, u0 = ffunc(xloc, ispcs, conf)

                                    #if not((x0[0] >= 0.0) and (x0[0]<=1.0) ):
                                    #    continue

                                    c.container.add_particle(x0, u0)

                # initialize analysis tiles ready for incoming simulation data
                for ip in range(conf.Nspecies):
                    c.addAnalysisSpecies()



# visualize particle content in vx direction
def plotXmesh(ax, n, conf, spcs, vdir):

    ax.clear()

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

    ax.minorticks_on()

    ax.set_xlim(n.getXmin(), n.getXmax())
    ax.set_ylim(-0.2, 0.2)



def filler(xloc, ispcs, conf):

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)

    ux = 0.0
    uy = 0.0
    uz = 0.0

    #if not((xx >= 200.0) and (xx<=300.0) ):
    #    x0 = [xx, 0.0, 0.0]
    #    u0 = [0.0, 0.0, 0.0]
    #    return x0, u0


    #speedtest
    #if 1.0 < xx < 2.0:
    #    x0 = [xx, 0.0, 0.0]
    #    u0 = [1.0, 0.0, 0.0]
    #    return x0, u0
    #else:
    #    x0 = [xx, 0.0, 0.0]
    #    u0 = [0.0, 0.0, 0.0]
    #    return x0, u0

    # current test
    #x0 = [xx, 0.0, 0.0]
    #u0 = [0.1, 0.0, 0.0]
    #return x0, u0


    vth = np.sqrt(conf.delgam)
    v0  = conf.gamma_e

    #Box-Muller sampling
    rr1 = np.random.rand()
    rr2 = np.random.rand()
    r1  = np.sqrt(-2.0*np.log(rr1))*np.cos(2.0*np.pi*rr2)
    r2  = np.sqrt(-2.0*np.log(rr1))*np.sin(2.0*np.pi*rr2)
    ux = r1*vth + v0

    #print("injecting into {} ({})".format(xx, xloc[0]))

    #perturb with wave
    Lx = conf.Nx*conf.NxMesh
    kmode = conf.modes
    mux_noise = conf.beta*np.cos(2.0*np.pi*kmode*xx/Lx) * (Lx/(2.0*np.pi*kmode))
    ux += vth*mux_noise


    x0 = [xx, 0.0, 0.0]
    u0 = [ux, uy, uz]
    return x0, u0


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(node, conf):

    #Lx  = conf.Nx*conf.NxMesh*conf.dx
    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    k = 2.0 #mode

    n0 = 1.0

    for i in range(node.getNx()):
        for j in range(node.getNy()):
            c = node.getCellPtr(i,j)
            yee = c.getYee(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        #get x_i+1/2 (Yee lattice so rho_i)
                        xloc0 = spatialLoc(node, (i,j), (l,  m,n), conf)
                        xloc1 = spatialLoc(node, (i,j), (l+1,m,n), conf)

                        xmid = 0.5*(xloc0[0] + xloc1[0])
                        #yee.ex[l,m,n] = n0*conf.me*conf.beta*np.sin(2.0*np.pi*k*xmid/Lx)/k

                        yee.ex[l,m,n] = 1.0e-5

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def plotDebug(ax, n, conf):

    #for i in range(conf.Nx):
    #    cid = n.cellId(i,0)
    #    c = n.getCellPtr(cid)

    #    x  = c.container.loc(0)
    #    ux = c.container.vel(0)
    #    ex = c.container.ex()

    #    ax.plot(x, ex, "k.", markersize=2, alpha=0.8)

    #smoothed ex
    yee = getYee(n, conf)
    ex = yee['ex']
    xx = yee['x']

    exS = smooth(ex, 10)
    ax.plot(xx, exS, "r-")



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
    plt.fig = plt.figure(1, figsize=(8,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(8, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    for ai in range(8):
        axs.append( plt.subplot(gs[ai]) )


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
    node = corgi.Node(conf.Nx, conf.Ny)

    xmin = 0.0
    #xmax = conf.dx*conf.Nx*conf.NxMesh 
    xmax = conf.Nx*conf.NxMesh #XXX scaled length
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

    #insert_em(node, conf)


    timer.stop("init") 
    timer.stats("init") 
    # end of initialization
    ################################################## 



    # visualize initial condition
    plotNode( axs[0], node, conf)
    plotXmesh(axs[1], node, conf, 0, "x")

    saveVisz(-1, node, conf)

    
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
    #
    #filtering

    pusher   = pypic.Pusher()
    fintp    = pypic.ParticleFieldInterpolator()
    comm     = pypic.Communicator()
    currint  = pypic.Depositer()
    analyzer = pypic.Analyzator()


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


        ##update boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.updateBoundaries(node)

        #interpolate fields
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                fintp.solve(cell)

        #pusher
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                pusher.solve(cell)

        #deposit current
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                currint.deposit(cell)

        ##update particle boundaries
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                comm.transfer(cell, node)

        #exchange currents
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.exchangeCurrents(node)

        #filter

        #add current to E
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                cell.depositCurrent()


        #analyze
        for j in range(node.getNy()):
            for i in range(node.getNx()):
                cell = node.getCellPtr(i,j)
                analyzer.analyze(cell)

        timer.lap("step")

        #save temporarily to file
        save(node, conf, ifile, f5)
        ifile += 1


        #I/O
        if (lap % conf.interval == 0):
            print("--------------------------------------------------")
            print("------ lap: {} / t: {}".format(lap, time)) 

            timer.stats("step")
            timer.start("io")


            plasma.writeYee(node,      lap, conf.outdir + "/")
            plasma.writeAnalysis(node, lap, conf.outdir + "/")
            #plasma.writeMesh(node,     lap, conf.outdir + "/")

            plotNode( axs[0], node, conf)
            plotXmesh(axs[1], node, conf, 0, "x")

            plotJ(    axs[5], node, conf)

            plotE(    axs[6], node, conf)
            plotDebug(axs[6], node, conf)

            plotDens( axs[7], node, conf)


            saveVisz(lap, node, conf)



            timer.stop("io")
            timer.stats("io")
            timer.start("step") #refresh lap counter (avoids IO profiling)

        time += conf.cfl/conf.c_omp
    #end of loop

    #node.finalizeMpi()


    timer.stop("total")
    timer.stats("total")
