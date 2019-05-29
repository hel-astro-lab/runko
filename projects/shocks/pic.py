from __future__ import print_function
from mpi4py import MPI

import numpy as np

import sys, os
import h5py

import pycorgi.twoD as corgi
import pyplasmabox.tools.twoD as pytools
import pyplasmabox.vlv.twoD as pyvlv
import pyplasmabox.pic.twoD as pypic
import pyplasmabox.fields.twoD as pyfld


from configSetup import Configuration
import argparse
import initialize as init
from initialize_pic import loadTiles
from initialize_pic import initialize_virtuals
from initialize_pic import globalIndx
from sampling import boosted_maxwellian
from initialize_pic import spatialLoc
from injector_pic import inject
from injector_pic import insert_em
from time import sleep
from visualize import get_yee

#global simulation seed
np.random.seed(1)

try:
    import matplotlib.pyplot as plt
    from visualize import plotNode
    from visualize import plotJ, plotE, plotDens
    from visualize import saveVisz
    
    from visualize import getYee2D
    from visualize import plot2dYee
    from visualize_pic import plot2dParticles

except:
    pass

from timer import Timer

from init_problem import Configuration_Shocks as Configuration

debug = False

def debug_print(n, msg):
    if debug:
        print("{}: {}".format(n.rank(), msg))
        sys.stdout.flush()


def filler(xloc, ispcs, conf):
    # perturb position between x0 + RUnif[0,1)

    #electrons
    if ispcs == 0:
        delgam  = conf.delgam #* np.abs(conf.mi / conf.me) * conf.temp_ratio

        xx = xloc[0] + np.random.rand(1)
        yy = xloc[1] + np.random.rand(1)
        #zz = xloc[2] + np.random.rand(1)
        zz = 0.5

    #positrons/ions/second species
    if ispcs == 1:
        delgam  = conf.delgam

        #on top of electrons
        xx = xloc[0]
        yy = xloc[1]
        zz = 0.5

    gamma = conf.gamma
    direction = -1
    ux, uy, uz, uu = boosted_maxwellian(delgam, gamma, direction=direction, dims=3)

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0


def initialize_piston(node, piston, conf):

    gam = conf.wallgamma
    beta = np.sqrt(1.-1./gam**2.)

    piston.gammawall = gam
    piston.betawall = beta
    piston.walloc = 5.0

    print("wall gamma:", piston.gammawall)
    print("wall beta:", piston.betawall)

    for cid in node.get_local_tiles():
        tile = node.get_tile(cid)

        #TODO: remove prtcls inside the wall




# Field initialization (guide field)
def insert_em(node, conf):

    #into radians
    btheta = conf.btheta/180.*np.pi
    bphi   = conf.bphi/180.*np.pi
    beta   = conf.beta

    kk = 0
    for cid in node.get_tile_ids():
        tile = node.get_tile(cid)
        yee = tile.get_yee(0)

        ii,jj = tile.index

        for n in range(conf.NzMesh):
            for m in range(-1, conf.NyMesh+1):
                for l in range(-1, conf.NxMesh+1):
                    # get global coordinates
                    iglob, jglob, kglob = globalIndx( (ii,jj), (l,m,n), conf)

                    yee.bx[l,m,n] = conf.binit*np.cos(btheta) 
                    yee.by[l,m,n] = conf.binit*np.sin(btheta)*np.sin(bphi)
                    yee.bz[l,m,n] = conf.binit*np.sin(btheta)*np.cos(bphi)   

                    yee.ex[l,m,n] = 0.0
                    yee.ey[l,m,n] =-beta*yee.bz[l,m,n]
                    yee.ez[l,m,n] = beta*yee.by[l,m,n]


if __name__ == "__main__":

    do_plots = True
    do_print = False

    if MPI.COMM_WORLD.Get_rank() == 0:
        do_print =True

    if do_print:
        print("Running with {} MPI processes.".format(MPI.COMM_WORLD.Get_size()))

    ################################################## 
    # set up plotting and figure
    try:
        if do_plots:
            plt.fig = plt.figure(1, figsize=(6,9))
            plt.rc('font', family='serif', size=8)
            plt.rc('xtick')
            plt.rc('ytick')
            
            gs = plt.GridSpec(11, 1)
            gs.update(hspace = 0.0)
            
            axs = []
            for ai in range(11):
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
    parser = argparse.ArgumentParser(description='Simple PIC-Maxwell simulations')
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
    xmax = conf.Nx*conf.NxMesh #XXX scaled length
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh
    node.set_grid_lims(xmin, xmax, ymin, ymax)

    #init.loadMpiRandomly(node)
    #init.loadMpiXStrides(node)
    debug_print(node, "load mpi 2d")
    init.loadMpi2D(node)
    debug_print(node, "load tiles")
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

        debug_print(node, "read")
        if do_print:
            print("...reading Yee lattices (lap {}) from {}".format(read_lap, odir))
        pyvlv.read_yee(node, read_lap, odir)

        if do_print:
            print("...reading particles (lap {}) from {}".format(read_lap, odir))
        pyvlv.read_particles(node, read_lap, odir)

        lap += 1 #step one step ahead

    # initialize
    if do_initialization:
        debug_print(node, "inject")
        lap = 0
        np.random.seed(1)
        inject(node, filler, conf) #injecting plasma particles
        insert_em(node, conf)

    #static load balancing setup; communicate neighbor info once
    debug_print(node, "analyze bcs")
    node.analyze_boundaries()
    debug_print(node, "send tiles")
    node.send_tiles()
    debug_print(node, "recv tiles")
    node.recv_tiles()
    MPI.COMM_WORLD.barrier()

    #sys.exit()

    debug_print(node, "init virs")
    initialize_virtuals(node, conf)


    timer.stop("init") 
    timer.stats("init") 


    # end of initialization
    ################################################## 
    debug_print(node, "solvers")


    # visualize initial condition
    if do_plots:
        try:
            plotNode( axs[0], node, conf)
            #plotXmesh(axs[1], node, conf, 0, "x")
            saveVisz(-1, node, conf)
        except:
            pass


    Nsamples = conf.Nt
    if conf.gammarad > 0:
        pusher   = pypic.BorisDragPusher()
        pusher.drag = conf.drag_amplitude
        pusher.temp = conf.radtemp
    else:
        pusher   = pypic.BorisPusher()


    fintp    = pypic.LinearInterpolator()
    currint  = pypic.ZigZag()
    analyzer = pypic.Analyzator()
    #flt     =  pytools.Filter(conf.NxMesh, conf.NyMesh)
    #flt.init_gaussian_kernel(2.0, 2.0)

    #moving walls
    piston   = pypic.Piston()
    initialize_piston(node, piston, conf)


    # quick field snapshots
    debug_print(node, "qwriter")
    qwriter  = pyfld.QuickWriter(conf.outdir, 
            conf.Nx, conf.NxMesh,
            conf.Ny, conf.NyMesh,
            conf.Nz, conf.NzMesh,
            conf.stride)

    # test particles
    debug_print(node, "tpwriter")
    tpwriter = pypic.TestPrtclWriter(
            conf.outdir, 
            conf.Nx, conf.NxMesh,
            conf.Ny, conf.NyMesh,
            conf.Nz, conf.NzMesh,
            conf.ppc, len(node.get_local_tiles()),
            conf.n_test_prtcls)




    debug_print(node, "mpi_e")
    node.send_data(1) 
    node.recv_data(1) 
    node.wait_data(1) 

    debug_print(node, "mpi_b")
    node.send_data(2) 
    node.recv_data(2) 
    node.wait_data(2) 

    ################################################## 
    sys.stdout.flush()

    #simulation loop
    time = lap*(conf.cfl/conf.c_omp)
    for lap in range(lap, conf.Nt+1):
        debug_print(node, "lap_start")

        ################################################## 
        # advance Half B

        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_b1")
        debug_print(node, "mpi_b1")

        node.send_data(2) 
        node.recv_data(2) 
        node.wait_data(2) 

        timer.stop_comp("mpi_b1")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc0")
        debug_print(node, "upd_bc0")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.update_boundaries(node)

        timer.stop_comp("upd_bc0")

        #--------------------------------------------------
        #push B half
        timer.start_comp("push_half_b1")
        debug_print(node, "push_half_b1")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.push_half_b()
            piston.field_bc(tile)

        timer.stop_comp("push_half_b1")

        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_b2")
        debug_print(node, "mpi_b2")

        node.send_data(2) 
        node.recv_data(2) 
        node.wait_data(2) 

        timer.stop_comp("mpi_b2")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc1")
        debug_print(node, "upd_bc1")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.update_boundaries(node)

        timer.stop_comp("upd_bc1")


        ################################################## 
        # move particles (only locals tiles)

        #--------------------------------------------------
        #interpolate fields (can move to next asap)
        timer.start_comp("interp_em")
        debug_print(node, "interp_em")

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            fintp.solve(tile)

        timer.stop_comp("interp_em")
        #--------------------------------------------------

        #--------------------------------------------------
        #push particles in x and u
        timer.start_comp("push")
        debug_print(node, "push")

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            pusher.solve(tile)

        timer.stop_comp("push")

        #--------------------------------------------------
        # apply moving walls
        timer.start_comp("walls")
        debug_print(node, "walls")
        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            piston.solve(tile)

        timer.stop_comp("walls")
        ##################################################
        # advance B half

        #--------------------------------------------------
        #push B half
        timer.start_comp("push_half_b2")
        debug_print(node, "push_half_b2")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.push_half_b()
            piston.field_bc(tile)

        timer.stop_comp("push_half_b2")


        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_e1")
        debug_print(node, "mpi_e1")

        node.send_data(1) 
        node.recv_data(1) 
        node.wait_data(1) 

        timer.stop_comp("mpi_e1")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc2")
        debug_print(node, "upd_bc2")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.update_boundaries(node)

        timer.stop_comp("upd_bc2")


        ##################################################
        # advance E 

        #--------------------------------------------------
        #push E
        timer.start_comp("push_e")
        debug_print(node, "push_e")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.push_e()
            piston.field_bc(tile)

        timer.stop_comp("push_e")

        #--------------------------------------------------
        #current calculation; charge conserving current deposition
        timer.start_comp("comp_curr")
        debug_print(node, "comp_cur")

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            currint.solve(tile)

        timer.stop_comp("comp_curr")

        #--------------------------------------------------
        timer.start_comp("clear_vir_cur")
        debug_print(node, "clear_vir_cur")

        for cid in node.get_virtual_tiles():
            tile = node.get_tile(cid)
            tile.clear_current()

        timer.stop_comp("clear_vir_cur")


        ##################################################
        # current solving is also taking place in nbor ranks
        # that is why we update virtuals here with MPI
        #
        # This is the most expensive task so we do not double it 
        # here.

        #--------------------------------------------------
        #mpi send currents
        timer.start_comp("mpi_cur")
        debug_print(node, "mpi_cur")

        node.send_data(0)
        node.recv_data(0) 
        node.wait_data(0)

        timer.stop_comp("mpi_cur")


        #--------------------------------------------------
        #exchange currents
        timer.start_comp("cur_exchange")
        debug_print(node, "cur_exchange")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.exchange_currents(node)

        timer.stop_comp("cur_exchange")



        ##################################################
        # particle communication (only local/boundary tiles)

        #--------------------------------------------------
        #local particle exchange (independent)
        timer.start_comp("check_outg_prtcls")
        debug_print(node, "check_outg_prtcls")

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            tile.check_outgoing_particles()

        timer.stop_comp("check_outg_prtcls")

        #--------------------------------------------------
        # global mpi exchange (independent)
        timer.start_comp("pack_outg_prtcls")
        debug_print(node, "pack_outg_prtcls")

        for cid in node.get_boundary_tiles():
            tile = node.get_tile(cid)
            tile.pack_outgoing_particles()

        timer.stop_comp("pack_outg_prtcls")

        #--------------------------------------------------
        # MPI global particle exchange
        # transfer primary and extra data
        timer.start_comp("mpi_prtcls")
        debug_print(node, "mpi_prtcls")

        debug_print(node, "mpi_prtcls: send3")
        node.send_data(3) 

        debug_print(node, "mpi_prtcls: recv3")
        node.recv_data(3) 

        debug_print(node, "mpi_prtcls: wait3")
        node.wait_data(3) 

        # orig just after send3
        debug_print(node, "mpi_prtcls: send4")
        node.send_data(4) 

        debug_print(node, "mpi_prtcls: recv4")
        node.recv_data(4) 

        debug_print(node, "mpi_prtcls: wait4")
        node.wait_data(4) 

        timer.stop_comp("mpi_prtcls")

        #--------------------------------------------------
        # global unpacking (independent)
        timer.start_comp("unpack_vir_prtcls")
        debug_print(node, "unpack_vir_prtcls")

        for cid in node.get_virtual_tiles(): 
            tile = node.get_tile(cid)
            tile.unpack_incoming_particles()
            tile.check_outgoing_particles()

        timer.stop_comp("unpack_vir_prtcls")

        #--------------------------------------------------
        # transfer local + global
        timer.start_comp("get_inc_prtcls")
        debug_print(node, "get_inc_prtcls")

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            tile.get_incoming_particles(node)

        timer.stop_comp("get_inc_prtcls")

        #--------------------------------------------------
        # delete local transferred particles
        timer.start_comp("del_trnsfrd_prtcls")
        debug_print(node, "del_trnsfrd_prtcls")

        for cid in node.get_local_tiles():
            tile = node.get_tile(cid)
            tile.delete_transferred_particles()

        timer.stop_comp("del_trnsfrd_prtcls")

        #--------------------------------------------------
        # delete all virtual particles (because new prtcls will come)
        timer.start_comp("del_vir_prtcls")
        debug_print(node, "del_vir_prtcls")

        for cid in node.get_virtual_tiles(): 
            tile = node.get_tile(cid)
            tile.delete_all_particles()

        timer.stop_comp("del_vir_prtcls")



        ##################################################
        #filter
        #timer.start_comp("filter")

        #for cid in node.get_tile_ids():
        #    tile = node.get_tile(cid)
        #    flt.get_padded_current(tile, node)

        #    #flt.fft_image_forward()
        #    #flt.apply_kernel()
        #    #flt.fft_image_backward()
        #
        #    for fj in range(conf.npasses):
        #        flt.direct_convolve_3point()
        #    flt.set_current(tile)

        ##cycle new and temporary currents (only if filttering)
        #for cid in node.get_tile_ids():
        #    tile = node.get_tile(cid)
        #    tile.cycle_current()
        # FIXME: cycle also virtuals
        #--------------------------------------------------
        #timer.stop_comp("filter")


        #--------------------------------------------------
        #add current to E
        timer.start_comp("add_cur")
        debug_print(node, "add_cur")

        for cid in node.get_tile_ids():
            tile = node.get_tile(cid)
            tile.deposit_current()

        timer.stop_comp("add_cur")


        ##################################################
        # data reduction and I/O

        timer.lap("step")
        if (lap % conf.interval == 0):
            debug_print(node, "io")
            if do_print:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time)) 

            #for cid in node.get_tile_ids():
            #    tile = node.get_tile(cid)
            #    tile.erase_temporary_arrays()

            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()
            
            #analyze (independent)
            timer.start("io")

            for cid in node.get_local_tiles():
                tile = node.get_tile(cid)
                analyzer.analyze2d(tile)

            # barrier for quick writers
            MPI.COMM_WORLD.barrier()

            debug_print(node, "qwriter")
            #shallow IO
            qwriter.write(node, lap) #quick field snapshots

            debug_print(node, "tpwriter")
            tpwriter.write(node, lap) #test particles

            #deep IO
            if (conf.full_interval != -1 and (lap % conf.full_interval == 0) and (lap > 0)):
                debug_print(node, "deep io")

                debug_print(node, "deep write_yee")
                pyvlv.write_yee(node,      lap, conf.outdir + "/full_output/" )
                debug_print(node, "deep write_analysis")
                pyvlv.write_analysis(node, lap, conf.outdir + "/full_output/" )
                debug_print(node, "deep write_prtcls")
                pyvlv.write_particles(node,lap, conf.outdir + "/full_output/" )


            #restart IO (overwrites)
            if ((lap % conf.restart == 0) and (lap > 0)):
                debug_print(node, "restart io")
                #flip between two sets of files
                deep_io_switch = 1 if deep_io_switch == 0 else 0

                print(node, "write_yee")
                pyvlv.write_yee(node,      deep_io_switch, conf.outdir + "/restart/" )
                print(node, "write_prtcls")
                pyvlv.write_particles(node,deep_io_switch, conf.outdir + "/restart/" )

                #if successful adjust info file
                MPI.COMM_WORLD.barrier() # sync everybody in case of failure before write
                if node.rank() == 0:
                    with open(conf.outdir+"/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(lap, deep_io_switch))


            #--------------------------------------------------
            #2D plots
            if do_plots:
                try:
                    plotNode(axs[0], node, conf)

                    yee = getYee2D(node, conf)
                    plot2dYee(axs[1],  yee, node, conf, 'rho', label_title=True)
                    plot2dYee(axs[2],  yee, node, conf, 'jx' , label_title=True)
                    plot2dYee(axs[3],  yee, node, conf, 'jy' , label_title=True)
                    plot2dYee(axs[4],  yee, node, conf, 'jz' , label_title=True)
                    plot2dYee(axs[5],  yee, node, conf, 'ex' , label_title=True)
                    plot2dYee(axs[6],  yee, node, conf, 'ey' , label_title=True)
                    plot2dYee(axs[7],  yee, node, conf, 'ez' , label_title=True)
                    plot2dYee(axs[8],  yee, node, conf, 'bx' , label_title=True)
                    plot2dYee(axs[9],  yee, node, conf, 'by' , label_title=True)
                    plot2dYee(axs[10], yee, node, conf, 'bz' , label_title=True)
                    saveVisz(lap, node, conf)
                except:
                    #print()
                    pass
            timer.stop("io")


            timer.stats("io")
            timer.start("step") #refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        #next step
        time += conf.cfl/conf.c_omp
    #end of loop

    timer.stop("total")
    timer.stats("total")
