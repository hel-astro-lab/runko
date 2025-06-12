# -*- coding: utf-8 -*-

from mpi4py import MPI
import numpy as np
import sys, os

# runko + auxiliary modules
import pytools  # runko python tools

# problem specific modules
from init_problem import Configuration_Turbulence as Configuration
from init_problem import velocity_profile
from init_problem import density_profile
from init_problem import weigth_profile


#--------------------------------------------------
rnd_seed_default = 1
np.random.seed(rnd_seed_default)  # global simulation seed

#--------------------------------------------------
# Field initialization (guide field)
def insert_em_fields(grid, conf, do_initialization):
    if not do_initialization:
        return

    for tile in pytools.tiles_all(grid):
        g = tile.yee_lattice()

        if not conf.use_maxwell_split: # if no static component

            g.set_E(lambda i, j, k: (0, 0, 0))
            g.set_B(lambda i, j, k: (0, 0, conf.binit))

        elif conf.use_maxwell_split:
            raise NotImplementedError()

    return

#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
if __name__ == "__main__":

    # --------------------------------------------------
    # initialize auxiliary tools
    sch  = pytools.Scheduler() # Set MPI aware task scheduler
    args = pytools.parse_args() # parse command line arguments
    tplt = pytools.TerminalPlot(25, 25) # terminal plotting tool

    # create conf object with simulation parameters based on them
    conf = Configuration(args.conf_filename, do_print=sch.is_master)
    sch.conf = conf # remember to update conf to scheduler

    if conf.mpi_task_mode: sch.switch_to_task_mode() # switch scheduler to task mode

    # --------------------------------------------------
    # Timer for profiling
    timer = pytools.Timer()
    timer.start("total")
    timer.start("init")
    timer.do_print = sch.is_example_worker
    timer.verbose = 0  # 0 normal; 1 - debug mode

    sch.timer = timer # remember to update scheduler

    # --------------------------------------------------
    # create output folders
    if sch.is_master: pytools.create_output_folders(conf)

    # --------------------------------------------------
    # load runko
    if conf.threeD:
        # 3D modules
        import pycorgi.threeD as pycorgi   # corgi ++ bindings
        import pyrunko.pic2.threeD as pypic2 # pic2 c++ bindings
        import pyrunko.emf2.threeD as pyfld2 # fld2 c++ bindings
    else:
        raise NotImplementedError()


    # --------------------------------------------------
    # setup grid
    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
    sch.grid = grid # remember to update scheduler

    # compute initial mpi ranks using Hilbert's curve partitioning
    if conf.use_injector:
        pytools.load_catepillar_track_mpi(grid, conf.mpi_track, conf)
    else:
        pytools.balance_mpi(grid, conf) #Hilbert curve optimization

    # load pic tiles into grid (functionally same as pytools.pic.load_tiles(grid, conf)
    for k in range(grid.get_Nz()):
        for j in range(grid.get_Ny()):
            for i in range(grid.get_Nx()):
                grid.add_tile(pypic2.Tile(conf), (i, j, k))

    raise NotImplementedError("Rest is not changed from original pic.py")

    # --------------------------------------------------
    # simulation restart

    # get current restart file status
    io_stat = pytools.check_for_restart(conf)


    # no restart file; initialize simulation
    if io_stat["do_initialization"]:
        if sch.is_master: print("initializing simulation...")
        lap = 0

        rseed = MPI.COMM_WORLD.Get_rank()
        np.random.seed(rseed + rnd_seed_default)  # sync rnd generator seed for different mpi ranks

        # injecting plasma particles
        if not(conf.use_injector):
            prtcl_stat = pytools.pic2.inject(grid, velocity_profile, density_profile, conf)
            if sch.is_example_worker:
                print
                print("     e- prtcls: {}".format(prtcl_stat[0]))
                print("     e+ prtcls: {}".format(prtcl_stat[1]))

        # inserting em grid
        insert_em_fields(grid, conf, do_initialization=True)

        # save a snapshot of the state to disk
        pytools.save_mpi_grid_to_disk(conf.outdir, 0, grid, conf)

    else:
        if sch.is_master:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        # insert damping reference values
        insert_em_fields(grid, conf, do_initialization=False)

        # read restart files
        pyfld2.read_grids(grid, io_stat["read_lap"], io_stat["read_dir"])
        pypic2.read_particles(grid, io_stat["read_lap"], io_stat["read_dir"])

        # set particle types
        print("NOT IMPLEMENTED: set particle types (should be done in ctor of pic2 tile)")

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    # static load balancing setup; communicate neighborhood info once

    if sch.is_master: print("load balancing grid..."); sys.stdout.flush()

    # --------------------------------------------------
    # initial driving; NOTE: turbulence specific/ needs to be here BEFORE update_boundaries call

    if False: # driven setup
        raise NotImplementedError("Driver setup not implemented for pic2.")
    else: # decaying setup
        print("NOT IMPLEMENTED: decay setup for pic2.")

    # update boundaries
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    if sch.is_master: print("loading virtual tiles..."); sys.stdout.flush()

    # load virtual mpi halo tiles
    print("NOT IMPLEMENTED: loading virtual tiles for pic2.")

    # --------------------------------------------------
    # load physics solvers

    if sch.is_master: print("loading solvers..."); sys.stdout.flush()

    print("NOT IMPLEMENTED: emf2 field propagator.")
    # sch.fldpropE = pyfld.FDTD2()
    # sch.fldpropB = pyfld.FDTD2()

    # enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    # sch.fldpropE.corr = conf.c_corr
    # sch.fldpropB.corr = conf.c_corr


    # --------------------------------------------------
    # particle puhser

    print("NOT IMPLEMENTED: pic2 particle pusher.")
    # sch.pusher = pypic2.BorisPusher()

    # Do pusher configuration in pusher ctor.

    #if conf.gammarad > 0:
    #    sch.pusher   = pypic.BorisDragPusher()
    #    sch.pusher.drag = conf.drag_amplitude
    #    sch.pusher.temp = 0.0

    # background field from external pusher components
    # if conf.use_maxwell_split:
    #     sch.pusher.bx_ext = conf.bx_ext
    #     sch.pusher.by_ext = conf.by_ext
    #     sch.pusher.bz_ext = conf.bz_ext

    #     sch.pusher.ex_ext = conf.ex_ext
    #     sch.pusher.ey_ext = conf.ey_ext
    #     sch.pusher.ez_ext = conf.ez_ext


    # --------------------------------------------------
    # particle interpolator
    print("NOT IMPLEMENTED: pic2 field interpolator.")
    # sch.fintp = pypic2.LinearInterpolator()

    # --------------------------------------------------
    # current deposit
    print("NOT IMPLEMENTED: pic2 current deposit.")
    # sch.currint = pypic2.ZigZag()

    # --------------------------------------------------
    # filter
    print("NOT IMPLEMENTED: emf2 filter.")
    # sch.flt = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    # --------------------------------------------------
    # I/O objects
    if sch.is_master: print("loading IO objects..."); sys.stdout.flush()

    # quick field snapshots
    print("NOT IMPLEMENTED: emf2 field writer.")
    print("NOT IMPLEMENTED: pic2 particle writers.")
    print("NOT IMPLEMENTED: pic2 momentum writer.")
    print("NOT IMPLEMENTED: emf2 field slice writer.")

    # --------------------------------------------------
    # --------------------------------------------------
    # --------------------------------------------------
    # end of initialization
    timer.stop("init")
    timer.stats("init")

    if sch.is_master: print('init: starting simulation...'); sys.stdout.flush()


    ##################################################
    # simulation time step loop

    # simulation loop
    time = lap * (conf.cfl / conf.c_omp)
    for lap in range(lap, conf.Nt + 1):

        # Commented out wall operations present in orginal pic.py are removed.

        # --------------------------------------------------
        # comm E and B
        sch.operate( dict(name='mpi_b0', solver='mpi', method='b', ) )
        sch.operate( dict(name='mpi_e0', solver='mpi', method='e', ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid,[1,2] ], nhood='local', ) )

        # --------------------------------------------------
        # push B half
        sch.operate( dict(name='push_half_b1', solver='fldpropB', method='push_half_b', nhood='local',) )

        # comm B
        sch.operate( dict(name='mpi_b1',  solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc ', solver='tile',method='update_boundaries',args=[grid, [2,] ], nhood='local',) )

        # --------------------------------------------------
        # move particles (only locals tiles)

        # interpolate fields and push particles in x and u
        sch.operate( dict(name='interp_em', solver='fintp',  method='solve', nhood='local', ) )
        #sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', ) ) # all particle continers
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[0]) ) # e^-
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[1]) ) # e^+

        sch.operate( dict(name='clear_cur', solver='tile',   method='clear_current', nhood='all', ) )

        # --------------------------------------------------
        # advance half B
        sch.operate( dict(name='push_half_b2', solver='fldpropB', method='push_half_b', nhood='local', ) )

        # comm B
        sch.operate( dict(name='mpi_b2', solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid, [2,] ], nhood='local',) )


        # --------------------------------------------------
        # push E
        sch.operate( dict(name='push_e',   solver='fldpropE', method='push_e',  nhood='local', ) )

        # TODO current deposit + MPI was here (never in pic2.py)

        # --------------------------------------------------
        # particle communication (only local/boundary tiles)

        # This should be refactored into a one C++ method.

        # local and global particle exchange
        sch.operate( dict(name='check_outg_prtcls',     solver='tile',  method='check_outgoing_particles',     nhood='local', ) )
        sch.operate( dict(name='pack_outg_prtcls',      solver='tile',  method='pack_outgoing_particles',      nhood='boundary', ) )

        sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p1',                           nhood='all', ) )
        sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p2',                           nhood='all', ) )

        sch.operate( dict(name='unpack_vir_prtcls',     solver='tile',  method='unpack_incoming_particles',    nhood='virtual', ) )
        sch.operate( dict(name='check_outg_vir_prtcls', solver='tile',  method='check_outgoing_particles',     nhood='virtual', ) )
        sch.operate( dict(name='get_inc_prtcls',        solver='tile',  method='get_incoming_particles',       nhood='local', args=[grid,]) )

        sch.operate( dict(name='del_trnsfrd_prtcls',    solver='tile',  method='delete_transferred_particles', nhood='local', ) )
        sch.operate( dict(name='del_vir_prtcls',        solver='tile',  method='delete_all_particles',         nhood='virtual', ) )

        # --------------------------------------------------
        # current calculation; charge conserving current deposition

        # These should be refactored into a one C++ method(?).

        # clear virtual current arrays for boundary addition after mpi, send currents, and exchange between tiles
        sch.operate( dict(name='comp_curr', solver='currint', method='solve', nhood='local', ) )
        sch.operate( dict(name='clear_vir_cur', solver='tile',method='clear_current',     nhood='virtual', ) )
        sch.operate( dict(name='mpi_cur',       solver='mpi', method='j',                 nhood='all', ) )
        sch.operate( dict(name='cur_exchange',  solver='tile',method='exchange_currents', nhood='local', args=[grid,], ) )

        # --------------------------------------------------
        # filter
        for fj in range(conf.npasses):

            # flt uses halo=2 padding so only every 3rd (0,1,2) pass needs update
            if fj % 2 == 0:
                sch.operate( dict(name='mpi_cur_flt', solver='mpi', method='j', ) )
                sch.operate( dict(name='upd_bc',      solver='tile',method='update_boundaries',args=[grid, [0,] ], nhood='local', ) )
                MPI.COMM_WORLD.barrier()
            sch.operate( dict(name='filter', solver='flt', method='solve', nhood='local', ) )


        # --------------------------------------------------
        # add antenna contribution

        # For Langevin antenna only.
        # sch.antenna.update_rnd_phases()
        # antenna.get_brms(grid) # debug tracking
        # sch.operate( dict(name='add_antenna', solver='antenna', method='add_ext_cur', nhood='local', ) )

        # --------------------------------------------------
        # add current to E
        sch.operate( dict(name='add_cur', solver='tile', method='deposit_current', nhood='local', ) )


        ##################################################
        # data reduction and I/O

        timer.lap("step")
        if lap % conf.interval == 0:
            if sch.is_master:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time))

            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()

            # io/analyze (independent)
            timer.start("io")
            # shrink particle arrays
            for tile in pytools.tiles_all(grid):
                tile.shrink_to_fit_all_particles()

            # barrier for quick writers
            MPI.COMM_WORLD.barrier()

            # shallow IO
            # NOTE: do moms before other IOs to keep rho field up-to-date
            mom_writer.write(grid, lap)  # pic distribution moments;
            fld_writer.write(grid, lap)  # quick field snapshots

            for pw in prtcl_writers:
                pw.write(grid, lap)  # particle tracking

            # pytools.save_mpi_grid_to_disk(conf.outdir, lap, grid, conf) # MPI grid

            #box peripheries
            if conf.threeD:
                slice_xy_writer.write(grid, lap)
                slice_xz_writer.write(grid, lap)
                slice_yz_writer.write(grid, lap)

            #--------------------------------------------------
            # terminal plot
            if sch.is_master:
                tplt.col_mode = False # use terminal color / use ASCII art
                if conf.threeD:
                    tplt.plot_panels( (1,1),
                    dict(axs=(0,0), data=slice_xy_writer.get_slice(8)/conf.j_norm ,   name='jz (xy)', cmap='RdBu'   ,vmin=-2, vmax=2),
                    #dict(axs=(0,0), data=slice_xy_writer.get_slice(9)/conf.p_norm   , name='ne (xy)', cmap='viridis',vmin= 0, vmax=4),
                    #dict(axs=(1,0), data=slice_xy_writer.get_slice(3)/conf.b_norm ,   name='bx (xy)', cmap='RdBu'   ,vmin=-2, vmax=2),
                    #dict(axs=(1,1), data=slice_xy_writer.get_slice(5)/conf.b_norm ,   name='bz (xy)', cmap='RdBu'   ,vmin=-2, vmax=2),
                    )

            #--------------------------------------------------
            #print statistics
            if sch.is_master:
                print('simulation time    {:7d} ({:7.1f} omp)   {:5.1f}%'.format( int(lap), time, 100.0*lap/conf.Nt))
                # For Langevin antenna only.
                # print('sim time {:7d} ({:7.1f} omp) ({:7.2f} l0/c) {:5.2f}%'.format(int(lap),
                                                                                        time,
                                                                                        sch.antenna.tcur,
                                                                                        100.0*lap/conf.Nt))

            #--------------------------------------------------
            # deep IO
            if conf.full_interval > 0 and (lap % conf.full_interval == 0) and (lap > 0):
                pyfld2.write_grids(grid, lap, conf.outdir + "/full_output/")
                pypic2.write_particles(grid, lap, conf.outdir + "/full_output/")


            # restart IO (overwrites)
            if (lap % conf.restart == 0) and (lap > 0):

                # flip between two sets of files
                io_stat["deep_io_switch"] = 1 if io_stat["deep_io_switch"] == 0 else 0

                pyfld2.write_grids(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                pypic2.write_particles(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                # if successful adjust info file
                MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write
                if grid.rank() == 0:
                    with open(conf.outdir + "/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(
                            lap,
                            io_stat["deep_io_switch"]+io_stat['restart_num']))

            MPI.COMM_WORLD.barrier() # extra barrier to synch everybody after IOs

            timer.stop("io")
            timer.stats("io")
            #--------------------------------------------------

            timer.start("step")  # refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        # next step
        time += conf.cfl / conf.c_omp
        #MPI.COMM_WORLD.barrier() # extra barrier to synch everybody
    # end of loop


    # --------------------------------------------------
    # end of simulation

    timer.stop("total")
    timer.stats("total")
