from __future__ import print_function
from mpi4py import MPI

import numpy as np
import sys, os
import h5py

import pyplasmabox.ffe.twoD as pyffe
import pycorgi.twoD as corgi
import pyplasmabox.vlv.twoD as pyvlv
import pyplasmabox.fields.twoD as pyfld

from timer import Timer

#from configSetup import Configuration
from init_problem import Configuration_Reconnection as Configuration
import argparse
import initialize as init
from initialize_ffe import loadTiles, initialize_virtuals
from initialize_pic import globalIndx


from visualize import get_yee
try:
    import matplotlib.pyplot as plt
    from visualize import plotNode
    from visualize import plotJ, plotE
    from visualize import saveVisz
    
    from visualize import getYee2D
    from visualize import plot2dYee
except:
    pass

from numpy import sinh, cosh, tanh, pi, sin, cos, tan, sqrt



# Field initialization (guide field)
def insert_em(node, conf):

    kk = 0
    for cid in node.get_tile_ids():
        tile = node.get_tile(cid)
        yee = tile.get_yee(0)

        ii,jj = tile.index

        for n in range(conf.NzMesh):
            for m in range(-1, conf.NyMesh+1):
                for l in range(-1, conf.NxMesh+1):
                    
                    # FIXME
                    #print("{}: ind ({},{},{})".format(node.rank(), l,m,n))
                    # get global coordinates
                    iglob, jglob, kglob = globalIndx( (ii,jj), (l,m,n), conf)

                    yee.ex[l,m,n] = 0.0
                    yee.ey[l,m,n] = 0.0
                    yee.ez[l,m,n] = 0.0

                    yee.bx[l,m,n] = 0.0
                    yee.by[l,m,n] = 0.0
                    yee.bz[l,m,n] = conf.binit   


##################################################
# Field initialization   
def insert_em_harris_sheet(node, conf):
    dvstripe = conf.dvstripe
    dstripe  = conf.dstripe
    nstripe  = conf.nstripe

    binit  = conf.binit
    btheta = conf.btheta
    bphi   = conf.bphi

    beta = conf.beta

    mxhalf  = conf.mxhalf
    myhalf  = conf.myhalf
    mzhalf  = conf.mzhalf
    lstripe = conf.lstripe

    kk = 0
    for jj in range(node.get_Ny()):
        for ii in range(node.get_Nx()):
            if node.get_mpi_grid(ii,jj) == node.rank():
                c = node.get_tile(ii,jj)
                yee = c.get_yee(0)

                for n in range(conf.NzMesh):
                    for m in range(conf.NyMesh):
                        for l in range(conf.NxMesh):

                            # get global coordinates
                            iglob, jglob, kglob = globalIndx( (ii,jj), (l,m,n), conf)

                            triggerz = 1.0
                            #triggerz = cosh(2.0*pi*(kglob - mzhalf)*dvstripe) #3D

                            velstripe = tanh(2.*pi*(iglob-mxhalf)*dvstripe)

                            if not(conf.periodicx):
                                stripetanh = tanh(2.*pi*(iglob-mxhalf)*dstripe)
                            else:
                                stripetanh = tanh(dstripe*lstripe*sin(2.*pi*(iglob - mxhalf)/lstripe))

                            yee.bx[l,m,n] = conf.bxinit*conf.binit   
                            if conf.trigger:
                                yee.by[l,m,n] = binit*sin(bphi)*stripetanh +  binit*btheta*cos(bphi)* \
             (1.-1./(cosh(2.*pi*(jglob-myhalf)*dvstripe)*cosh(2.*pi*(iglob-mxhalf)*dstripe)*triggerz))

                                yee.bz[l,m,n] = binit*cos(bphi)*stripetanh + binit*btheta*sin(bphi)*    \
             (1.-1./(cosh(2.*pi*(jglob-myhalf)*dvstripe)*cosh(2.*pi*(iglob-mxhalf)*dstripe)*triggerz))

                                yee.ey[l,m,n] = (-beta)*velstripe*yee.bz[l,m,n] 

                                #drive to trigger reconnection in the middle of the box; 
                                #the coefficient should be the desired ExB speed
                                yee.ez[l,m,n] = (+beta)*velstripe*yee.by[l,m,n]  \
                                 + conf.trigger_field *  \
                                      abs(yee.by[l,m,n])/(cosh(2.*pi*(jglob-myhalf)*dvstripe)*
                                                          cosh(2.*pi*(iglob-mxhalf)*dstripe)*triggerz)
                            else:
                                yee.by[l,m,n] = binit*sin(bphi)*stripetanh + binit*btheta*cos(bphi)             
                                yee.bz[l,m,n] = binit*cos(bphi)*stripetanh + binit*btheta*sin(bphi)               
                                yee.ey[l,m,n] = (-beta)*velstripe*yee.bz[l,m,n] 
                                yee.ez[l,m,n] = (+beta)*velstripe*yee.by[l,m,n]
                                                       
                            yee.ex[l,m,n] = 0.0

                            #print("ijk {}/{}/{} binit: {} stripetanh: {} beta: {} velstripe {}".format(
                            #        iglob, jglob,kglob,binit, stripetanh, beta, velstripe))


                # copy values to boundary cells
                try:
                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                c.ex_ref[l,m,n] = yee.ex[l,m,n]
                                c.ey_ref[l,m,n] = yee.ey[l,m,n]
                                c.ez_ref[l,m,n] = yee.ez[l,m,n]

                                c.bx_ref[l,m,n] = yee.bx[l,m,n]
                                c.by_ref[l,m,n] = yee.by[l,m,n]
                                c.bz_ref[l,m,n] = yee.bz[l,m,n]
                except:
                    #print("cell ({},{}) is not boundary cell".format(ii,jj))
                    pass


if __name__ == "__main__":

    np.random.seed(42)

    do_plots = True
    do_print = False
    if MPI.COMM_WORLD.Get_rank() == 0:
        do_print =True

    ################################################## 
    # set up plotting and figure
    try:
        if do_plots:
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
        #print()
        pass


    # Timer for profiling
    timer = Timer()
    timer.start("total")
    timer.start("init")
    timer.do_print = do_print


    # parse command line arguments
    parser = argparse.ArgumentParser(description='Force-Free MHD simulation driver')
    parser.add_argument('--conf', dest='conf_filename', default=None,
                       help='Name of the configuration file (default: None)')
    args = parser.parse_args()
    if args.conf_filename == None:
        conf = Configuration('config-test.ini', do_print=do_print) 
    else:
        if do_print:
            print("Reading configuration setup from ", args.conf_filename)
        conf = Configuration(args.conf_filename, do_print=do_print)

    node = corgi.Node(conf.Nx, conf.Ny, conf.Nz)

    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh

    node.set_grid_lims(xmin, xmax, ymin, ymax)
    init.loadMpi2D(node)

    loadTiles(node, conf)

    ################################################## 
    # Path to be created 
    if node.master():
        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)
        if not os.path.exists( conf.outdir+"/restart" ):
            os.makedirs(conf.outdir+"/restart")
        if not os.path.exists( conf.outdir+"/full_output" ):
            os.makedirs(conf.outdir+"/full_output")


    do_initialization = True

    #check if this is the first time and we do not have any restart files
    if not os.path.exists( conf.outdir+'/restart/laps.txt'):
        conf.laprestart = -1 #to avoid next if statement

    # restart from latest file
    deep_io_switch = 0
    if conf.laprestart >= 0:
        do_initialization = False

        #switch between automatic restart and user-defined lap
        if conf.laprestart == 0:
            #get latest restart file from housekeeping file
            with open(conf.outdir+"/restart/laps.txt", "r") as lapfile:
                #lapfile.write("{},{}\n".format(lap, deep_io_switch))
                lines = lapfile.readlines()
                slap, sdeep_io_switch = lines[-1].strip().split(',')
                lap = int(slap)
                deep_io_switch = int(sdeep_io_switch)

            read_lap = deep_io_switch
            odir = conf.outdir + '/restart'

        elif conf.laprestart > 0:
            lap = conf.laprestart
            read_lap = lap
            odir = conf.outdir + '/full_output'

        if do_print:
            print("...reading Yee lattices (lap {}) from {}".format(read_lap, odir))
        pyvlv.read_yee(node, read_lap, odir)

        lap += 1 #step one step ahead


    # initialize
    if do_initialization:
        lap = 0
        np.random.seed(1)
        #insert_em(node, conf)
        insert_em_harris_sheet(node, conf)

    #static load balancing setup; communicate neighbor info once
    node.analyze_boundaries()
    node.send_tiles()
    node.recv_tiles()
    initialize_virtuals(node, conf)

    timer.stop("init") 
    timer.stats("init") 


    # end of initialization
    ################################################## 

    # visualize initial condition
    if do_plots:
        try:
            plotNode( axs[0], node, conf)
            #plotXmesh(axs[1], node, conf, 0, "x")
            saveVisz(-1, node, conf)
        except:
            pass


    # quick field snapshots
    qwriter  = pyfld.QuickWriter(conf.outdir, 
            conf.Nx, conf.NxMesh,
            conf.Ny, conf.NyMesh,
            conf.Nz, conf.NzMesh,
            conf.stride)

    node.send_data(1) 
    node.recv_data(1) 
    node.wait_data(1) 

    node.send_data(2) 
    node.recv_data(2) 
    node.wait_data(2) 

    ################################################## 
    sys.stdout.flush()

    #simulation loop
    time = lap*(conf.cfl/conf.c_omp)
    for lap in range(lap, conf.Nt+1):

        #Main loop
        #  advance B half
        #  move particles
        #  advance B half
        #  advance E
        #  reset current (j=0)
        #  deposit
        #  exchange particles
        #  exchange current
        #  filter
        #  add current
        #  inject particles from other processors
        #  enlarge domain
        #  reorder particles
        #  pause simulation if pause file exists


        ################################################## 
        # advance Half B

        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_b1")

        node.send_data(2) 
        node.recv_data(2) 
        node.wait_data(2) 

        timer.stop_comp("mpi_b1")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc0")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.update_boundaries(node)

        timer.stop_comp("upd_bc0")

        #--------------------------------------------------
        #push B half
        timer.start_comp("push_half_b1")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.push_half_b()

        timer.stop_comp("push_half_b1")

        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_b2")

        node.send_data(2) 
        node.recv_data(2) 
        node.wait_data(2) 

        timer.stop_comp("mpi_b2")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc1")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.update_boundaries(node)

        timer.stop_comp("upd_bc1")

        ##################################################
        # advance B half

        #--------------------------------------------------
        #push B half
        timer.start_comp("push_half_b2")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.push_half_b()

        timer.stop_comp("push_half_b2")


        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_e1")

        node.send_data(1) 
        node.recv_data(1) 
        node.wait_data(1) 

        timer.stop_comp("mpi_e1")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc2")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.update_boundaries(node)

        timer.stop_comp("upd_bc2")


        ##################################################
        # advance E 

        #--------------------------------------------------
        #push E
        timer.start_comp("push_e")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.push_e()

        timer.stop_comp("push_e")

        #--------------------------------------------------
        #current calculation; charge conserving current deposition
        #timer.start_comp("comp_curr")

        #for cid in node.get_local_tiles():
        #    tile = node.get_tile(cid)
        #    currint.solve(tile)

        #timer.stop_comp("comp_curr")

        #--------------------------------------------------
        timer.start_comp("clear_vir_cur")

        for cid in node.get_virtual_tiles():
            tile = node.get_tile(cid)
            tile.clear_current()

        timer.stop_comp("clear_vir_cur")


        ##################################################
        # current solving is also taking place in nbor ranks

        #mpi send currents
        timer.start_comp("mpi_cur")

        node.send_data(0)
        node.recv_data(0) 
        node.wait_data(0)

        timer.stop_comp("mpi_cur")

        #--------------------------------------------------
        #exchange currents
        timer.start_comp("cur_exchange")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.exchange_currents(node)

        timer.stop_comp("cur_exchange")


        #--------------------------------------------------
        #add current to E
        timer.start_comp("add_cur")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.deposit_current()

        timer.stop_comp("add_cur")


        ##################################################
        # data reduction and I/O

        timer.lap("step")
        if (lap % conf.interval == 0):
            if do_print:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time)) 

            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()
            
            #analyze (independent)
            timer.start("io")


            #shallow IO
            qwriter.write(node, lap) #quick field snapshots

            #deep IO
            if (conf.full_interval != -1 and (lap % conf.full_interval == 0) and (lap > 0)):
                pyvlv.write_yee(node,      lap, conf.outdir + "/full_output/" )

            #restart IO (overwrites)
            if ((lap % conf.restart == 0) and (lap > 0)):
                #flip between two sets of files
                deep_io_switch = 1 if deep_io_switch == 0 else 0

                pyvlv.write_yee(node,      deep_io_switch, conf.outdir + "/restart/" )

                #if successful adjust info file
                MPI.COMM_WORLD.barrier()
                if node.rank() == 0:
                    with open(conf.outdir+"/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(lap, deep_io_switch))

            #--------------------------------------------------
            #2D plots
            #try:
            if True:
                if do_plots:
                    plotNode(axs[0], node, conf)

                    yee = getYee2D(node, conf)
                    plot2dYee(axs[2],  yee, node, conf, 'rho')
                    plot2dYee(axs[3],  yee, node, conf, 'jx')
                    plot2dYee(axs[4],  yee, node, conf, 'jy')
                    plot2dYee(axs[5],  yee, node, conf, 'jz')
                    plot2dYee(axs[6],  yee, node, conf, 'ex')
                    plot2dYee(axs[7],  yee, node, conf, 'ey')
                    plot2dYee(axs[8],  yee, node, conf, 'ez')
                    plot2dYee(axs[9],  yee, node, conf, 'bx')
                    plot2dYee(axs[10], yee, node, conf, 'by')
                    plot2dYee(axs[11], yee, node, conf, 'bz')
                    saveVisz(lap, node, conf)
            #except:
            #    print()
            #    pass
            timer.stop("io")


            timer.stats("io")
            timer.start("step") #refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        #MPI.COMM_WORLD.barrier()
        #sleep(0.2)
        time += conf.cfl/conf.c_omp
    #end of loop


    timer.stop("total")
    timer.stats("total")

