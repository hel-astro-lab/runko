from __future__ import print_function
from mpi4py import MPI

import numpy as np
import sys, os
import h5py

import pyrunko.ffe.twoD as pyffe
import pycorgi.twoD as corgi
import pyrunko.vlv.twoD as pyvlv
import pyrunko.fields.twoD as pyfld

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

from antenna2 import Antenna


# Field initialization (guide field)
def insert_em(grid, conf):

    kk = 0
    for tile in tiles_all(grid): 
        yee = tile.get_yee(0)

        ii,jj = tile.index
        for n in range(conf.NzMesh):
            for m in range(-3, conf.NyMesh+3):
                for l in range(-3, conf.NxMesh+3):
                    
                    # FIXME
                    #print("{}: ind ({},{},{})".format(grid.rank(), l,m,n))
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
def insert_em_harris_sheet(grid, conf):
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
    for jj in range(grid.get_Ny()):
        for ii in range(grid.get_Nx()):
            if grid.get_mpi_grid(ii,jj) == grid.rank():
                c = grid.get_tile(ii,jj)
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


#all tile generator
def tiles_all(grid):
    for cid in grid.get_tile_ids():
        tile = grid.get_tile(cid)
        yield tile

#local tile generator
def tiles_local(grid):
    for cid in grid.get_local_tiles():
        tile = grid.get_tile(cid)
        yield tile



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

    grid = corgi.Grid(conf.Nx, conf.Ny, conf.Nz)

    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh

    grid.set_grid_lims(xmin, xmax, ymin, ymax)
    init.loadMpi2D(grid)

    loadTiles(grid, conf)

    ################################################## 
    # Path to be created 
    if grid.master():
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
        pyvlv.read_yee(grid, read_lap, odir)

        lap += 1 #step one step ahead


    # initialize
    if do_initialization:
        lap = 0
        np.random.seed(1)
        #insert_em(grid, conf)
        insert_em_harris_sheet(grid, conf)

    #static load balancing setup; communicate neighbor info once
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    initialize_virtuals(grid, conf)

    timer.stop("init") 
    timer.stats("init") 


    # end of initialization
    ################################################## 

    # visualize initial condition
    if do_plots:
        try:
            plotNode( axs[0], grid, conf)
            #plotXmesh(axs[1], grid, conf, 0, "x")
            saveVisz(-1, grid, conf)
        except:
            pass

    fldprop  = pyfld.FDTD2()
    driftcur = pyffe.DriftCurrent(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    # quick field snapshots
    qwriter  = pyfld.QuickWriter(conf.outdir, 
            conf.Nx, conf.NxMesh,
            conf.Ny, conf.NyMesh,
            conf.Nz, conf.NzMesh,
            conf.stride)

    ################################################## 
    # Langeving antenna
    if do_initialization:

        # direct B_{x,y} perturbation
        if False:
            antenna = Antenna(conf.min_mode, conf.max_mode, conf)
            for tile in tiles_local(grid): 
                antenna.add_driving(tile)


    grid.send_data(1) 
    grid.recv_data(1) 
    grid.wait_data(1) 

    grid.send_data(2) 
    grid.recv_data(2) 
    grid.wait_data(2) 

    ################################################## 
    sys.stdout.flush()

    #simulation loop
    time = lap*(conf.cfl/conf.c_omp)
    for lap in range(lap, conf.Nt+1):

        #FFE algorithm is:
        # Derivative
        # compute drift current
        # enforce E <= B 
        # subtract E_par

        #RK4 scheme is
        # k1 = dt * f(t         , y)
        # k2 = dt * f(t + 0.5*dt, y + 0.5*k1)
        # k3 = dt * f(t + 0.5*dt, y + 0.5*k2)
        # k4 = dt * f(t + dt    , y + k3)
        # vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6

        #Coefficients correspond to (dt, dy, y_i)
        for (dti, dyi, yi) in [
                (0.0, 0.0, 0),
                (0.5, 0.5, 1),
                (0.5, 0.5, 2),
                (1.0, 1.0, 3)]:

            #dt = 0.0 + dti
            #y = y0 + dyi*ys[yi]


            ################################################## 
            # compute current

            #--------------------------------------------------
            #update boundaries
            t1 = timer.start_comp("upd_bc0")

            for tile in tiles_all(grid): 
                tile.update_boundaries(grid)

            timer.stop_comp(t1)

            #--------------------------------------------------
            #compute drift current
            t1 = timer.start_comp("compute_current")

            for tile in tiles_all(grid): 
                driftcur.comp_drift_cur(tile)

            timer.stop_comp(t1)

            #--------------------------------------------------
            #add current to E
            t1 = timer.start_comp("add_cur")

            for tile in tiles_all(grid): 
                tile.deposit_current()

            timer.stop_comp(t1)


            #--------------------------------------------------
            #push E
            t1 = timer.start_comp("push_e")

            fldprop.dt = 1.0 #XXX time step 
            for tile in tiles_all(grid): 
                fldprop.push_e(tile)

            timer.stop_comp(t1)

            #--------------------------------------------------
            # comm E
            t1 = timer.start_comp("mpi_e0")

            grid.send_data(1) 
            grid.recv_data(1) 
            grid.wait_data(1) 

            timer.stop_comp(t1)

            #--------------------------------------------------
            #update boundaries
            t1 = timer.start_comp("upd_bc1")

            for tile in tiles_all(grid): 
                tile.update_boundaries(grid)

            timer.stop_comp(t1)


            if False:

                #--------------------------------------------------
                #push B half
                timer.start_comp("push_b0")

                for tile in tiles_all(grid): 
                    tile.cfl *= 2.0
                    fldprop.push_half_b(tile)
                    tile.cfl *= 0.5

                timer.stop_comp("push_b0")


                #--------------------------------------------------
                # comm B
                timer.start_comp("mpi_b1")

                grid.send_data(2) 
                grid.recv_data(2) 
                grid.wait_data(2) 

                timer.stop_comp("mpi_b1")

                #--------------------------------------------------
                #update boundaries
                timer.start_comp("upd_bc2")

                for tile in tiles_all(grid): 
                    tile.update_boundaries(grid)

                timer.stop_comp("upd_bc2")

                #TODO: enforce E.B=0
                #TODO: enforce E <= B

                #--------------------------------------------------
                #compute drift current
                timer.start_comp("compute_para_current")

                for tile in tiles_all(grid): 
                    driftcur.comp_parallel_cur(tile)

                timer.stop_comp("compute_para_current")

                #--------------------------------------------------
                # comm B
                timer.start_comp("mpi_e1")

                grid.send_data(1) 
                grid.recv_data(1) 
                grid.wait_data(1) 

                timer.stop_comp("mpi_e1")


                #TODO XXX --------------------------------------

                #--------------------------------------------------
                #current calculation; charge conserving current deposition
                #timer.start_comp("comp_curr")

                #for cid in grid.get_local_tiles():
                #    tile = grid.get_tile(cid)
                #    currint.solve(tile)

                #timer.stop_comp("comp_curr")

                #--------------------------------------------------
                #timer.start_comp("clear_vir_cur")

                #for cid in grid.get_virtual_tiles():
                #    tile = grid.get_tile(cid)
                #    tile.clear_current()

                #timer.stop_comp("clear_vir_cur")




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
            qwriter.write(grid, lap) #quick field snapshots

            #deep IO
            if (conf.full_interval != -1 and (lap % conf.full_interval == 0) and (lap > 0)):
                pyvlv.write_yee(grid,      lap, conf.outdir + "/full_output/" )

            #restart IO (overwrites)
            if ((lap % conf.restart == 0) and (lap > 0)):
                #flip between two sets of files
                deep_io_switch = 1 if deep_io_switch == 0 else 0

                pyvlv.write_yee(grid,      deep_io_switch, conf.outdir + "/restart/" )

                #if successful adjust info file
                MPI.COMM_WORLD.barrier()
                if grid.rank() == 0:
                    with open(conf.outdir+"/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(lap, deep_io_switch))

            #--------------------------------------------------
            #2D plots
            #try:
            if True:
                if do_plots:
                    plotNode(axs[0], grid, conf)

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

            sys.stdout.flush()

        #MPI.COMM_WORLD.barrier()
        #sleep(0.2)
        time += conf.cfl/conf.c_omp
    #end of loop


    timer.stop("total")
    timer.stats("total")

