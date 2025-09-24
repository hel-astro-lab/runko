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
def insert_em_fields(grid, conf, do_initialization=True):

    for tile in pytools.tiles_all(grid):
        g = tile.get_grids(0)

        ii,jj,kk = tile.index if conf.threeD else (*tile.index, 0)


        # insert values into Yee lattices; includes halos from -3 to n+3
        if do_initialization:
            for n in range(-3, conf.NzMesh + 3):
                for m in range(-3, conf.NyMesh + 3):
                    for l in range(-3, conf.NxMesh + 3):

                        if not(conf.use_maxwell_split): # if no static component

                            # get global coordinates
                            #iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                            #r = np.sqrt(iglob ** 2 + jglob ** 2 + kglob ** 2)

                            g.ex[l,m,n] = 0.0
                            g.ey[l,m,n] = 0.0
                            g.ez[l,m,n] = 0.0

                            g.bx[l,m,n] = 0.0
                            g.by[l,m,n] = 0.0 
                            g.bz[l,m,n] = conf.binit

                        elif conf.use_maxwell_split: # static component
                            1
                            # TODO
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
        import pyrunko.pic.threeD as pypic # pic c++ bindings
        import pyrunko.emf.threeD as pyfld # fld c++ bindings
    elif conf.twoD:
        # 2D modules
        import pycorgi.twoD as pycorgi     # corgi ++ bindings
        import pyrunko.pic.twoD as pypic   # runko pic c++ bindings
        import pyrunko.emf.twoD as pyfld   # runko fld c++ bindings

    # --------------------------------------------------
    # setup grid
    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
    sch.grid = grid # remember to update scheduler


    # compute initial mpi ranks using Hilbert's curve partitioning
    if conf.use_injector:
        pytools.load_catepillar_track_mpi(grid, conf.mpi_track, conf) 
    else:
        #pytools.load_mpi_y_strides(grid, conf) #equal x stripes

        # advanced Hilbert curve optimization (leave 4 first ranks empty)
        #if grid.size() > 24: # big run
        #    pytools.balance_mpi_3D_rootmem(grid, 4) 
        #else:

        pytools.balance_mpi(grid, conf) #Hilbert curve optimization 


    # load pic tiles into grid
    pytools.pic.load_tiles(grid, conf)
    #load_damping_tiles(grid, conf)

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
            prtcl_stat = pytools.pic.inject(grid, velocity_profile, density_profile, conf)
            if sch.is_example_worker: 
                print("injected:")
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
        pyfld.read_grids(grid, io_stat["read_lap"], io_stat["read_dir"])
        pypic.read_particles(grid, io_stat["read_lap"], io_stat["read_dir"])

        # set particle types
        for tile in pytools.tiles_all(grid):
            for ispcs in range(conf.Nspecies):
                container = tile.get_container(ispcs)
                container.type = conf.prtcl_types[ispcs] # name container

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    # static load balancing setup; communicate neighborhood info once

    if sch.is_master: print("load balancing grid..."); sys.stdout.flush()

    # --------------------------------------------------
    # initial driving; NOTE: turbulence specific/ needs to be here BEFORE update_boundaries call

    if True: #driven setup
        from antenna3d_langevin import Antenna  # 3D forced Langevin antenna setup

        sch.antenna = Antenna(conf.min_mode, conf.max_mode, conf)
        if io_stat["do_initialization"]:
            for tile in pytools.tiles_local(grid):
                sch.antenna.add_driving(tile)
        else:
            # advance antenna
            for lap_tmp in range(0, lap):
                sch.antenna.update_rnd_phases()

    else: # decaying setup

        if conf.twoD:
            from antenna2d import Antenna
            antenna = Antenna(conf.min_mode, conf.max_mode, conf)
        elif conf.threeD:
            from antenna3d import Antenna
            antenna = Antenna(conf.max_mode, 2, conf)

        # add \delta B perturbations
        for tile in pytools.tiles_local(grid):
            antenna.add_driving(tile)

    # update boundaries
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    if sch.is_master: print("loading virtual tiles..."); sys.stdout.flush()

    # load virtual mpi halo tiles
    pytools.pic.load_virtual_tiles(grid, conf)


    # --------------------------------------------------
    # load physics solvers

    if sch.is_master: print("loading solvers..."); sys.stdout.flush()


    sch.fldpropE = pyfld.FDTD2()
    sch.fldpropB = pyfld.FDTD2()
    #sch.fldpropE = pyfld.FDTD4()
    #sch.fldpropB = pyfld.FDTD4()

    # enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    sch.fldpropE.corr = conf.c_corr
    sch.fldpropB.corr = conf.c_corr


    # --------------------------------------------------
    # particle puhser

    #sch.pusher = pypic.BorisPusher()
    #sch.pusher = pypic.VayPusher()
    sch.pusher = pypic.HigueraCaryPusher()
    #sch.pusher  = pypic.rGCAPusher()

    #if conf.gammarad > 0:
    #    sch.pusher   = pypic.BorisDragPusher()
    #    sch.pusher.drag = conf.drag_amplitude
    #    sch.pusher.temp = 0.0

    # background field from external pusher components
    if conf.use_maxwell_split:
        sch.pusher.bx_ext = conf.bx_ext 
        sch.pusher.by_ext = conf.by_ext
        sch.pusher.bz_ext = conf.bz_ext

        sch.pusher.ex_ext = conf.ex_ext 
        sch.pusher.ey_ext = conf.ey_ext
        sch.pusher.ez_ext = conf.ez_ext


    # --------------------------------------------------
    # particle interpolator
    sch.fintp = pypic.LinearInterpolator()
    #sch.fintp = pypic.QuadraticInterpolator() #2nd order quadratic
    #sch.fintp = pypic.CubicInterpolator()     #3rd order cubic 3d
    #sch.fintp = pypic.QuarticInterpolator()   #4th order quartic; 3d

    # --------------------------------------------------
    # current deposit
    sch.currint = pypic.ZigZag()
    #sch.currint = pypic.ZigZag_2nd()
    #sch.currint = pypic.ZigZag_3rd()
    #sch.currint = pypic.ZigZag_4th()
    #sch.currint = pypic.Esikerpov_2nd() # 3d only
    #sch.currint = pypic.Esikerpov_4th() # 3d only

    # --------------------------------------------------
    #filter
    sch.flt = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)


    # --------------------------------------------------
    # I/O objects
    if sch.is_master: print("loading IO objects..."); sys.stdout.flush()

    # quick field snapshots
    #fld_writer = pyfld.MasterFieldsWriter(
    fld_writer = pyfld.FieldsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride,
    )


    # tracked particles; only works with no injector (id's are messed up when coming out from injector)
    prtcl_writers = []

    #for ispc in [0, 1, 2]: #electrons & positrons
    for ispc in [0]: #electrons
        mpi_comm_size = grid.size() if not(conf.mpi_task_mode) else grid.size() - 1
        n_local_tiles = int(conf.Nx*conf.Ny*conf.Nz/mpi_comm_size)
        #print("average n_local_tiles:", n_local_tiles, " / ", mpi_comm_size)

        prtcl_writer = pypic.TestPrtclWriter(
                conf.outdir,
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh,
                conf.ppc,
                n_local_tiles, #len(grid.get_local_tiles()),
                conf.n_test_prtcls,)
        prtcl_writer.ispc = ispc
        prtcl_writers.append(prtcl_writer)


    # momens of particle distribution
    #mom_writer = pypic.MasterPicMomentsWriter(
    mom_writer = pypic.PicMomentsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride_mom,
    )

    # 3D box peripherals
    if conf.threeD:
        st = 1 # stride
        slice_xy_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 0, 1)
        slice_xz_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 1, 1)
        slice_yz_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 2, 1)

        # location of the slice (grid index)
        slice_xy_writer.ind = int(0.0*conf.Lz) # bottom slice
        slice_xz_writer.ind = int(0.0*conf.Ly) # side wall 1
        slice_yz_writer.ind = int(0.0*conf.Lx) # side wall 2

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

        # --------------------------------------------------
        # comm E and B
        sch.operate( dict(name='mpi_b0', solver='mpi', method='b', ) )
        sch.operate( dict(name='mpi_e0', solver='mpi', method='e', ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid,[1,2] ], nhood='local', ) )

        # --------------------------------------------------
        # push B half
        sch.operate( dict(name='push_half_b1', solver='fldpropB', method='push_half_b', nhood='local',) )
        #sch.operate( dict(name='wall_bc',      solver='lwall',    method='field_bc',    nhood='local',) )

        # comm B
        sch.operate( dict(name='mpi_b1',  solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc ', solver='tile',method='update_boundaries',args=[grid, [2,] ], nhood='local',) )

        # --------------------------------------------------
        # move particles (only locals tiles)

        # interpolate fields and push particles in x and u
        sch.operate( dict(name='interp_em', solver='fintp',  method='solve', nhood='local', ) )
        #sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', ) )
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[0]) ) # e^-
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[1]) ) # e^+


        # clear currents; need to call this before wall operations since they can deposit currents too 
        sch.operate( dict(name='clear_cur', solver='tile',   method='clear_current', nhood='all', ) )

        # apply moving/reflecting walls
        #sch.operate( dict(name='walls',     solver='lwall', method='solve', nhood='local', ) )


        # --------------------------------------------------
        # advance half B 
        sch.operate( dict(name='push_half_b2', solver='fldpropB', method='push_half_b', nhood='local', ) )
        #sch.operate( dict(name='wall_bc',      solver='lwall',    method='field_bc',    nhood='local', ) )

        # comm B
        sch.operate( dict(name='mpi_b2', solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid, [2,] ], nhood='local',) )


        # --------------------------------------------------
        # push E
        sch.operate( dict(name='push_e',   solver='fldpropE', method='push_e',  nhood='local', ) )
        #sch.operate( dict(name='wall_bc',  solver='lwall',    method='field_bc',nhood='local', ) )

        # TODO current deposit + MPI was here

        # --------------------------------------------------
        # particle communication (only local/boundary tiles)

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
        sch.antenna.update_rnd_phases()
        #antenna.get_brms(grid) # debug tracking
        sch.operate( dict(name='add_antenna', solver='antenna', method='add_ext_cur', nhood='local', ) )

        # --------------------------------------------------
        # add current to E
        sch.operate( dict(name='add_cur', solver='tile', method='deposit_current', nhood='local', ) )
        #operate( dict(name='wall_bc', solver='lwall', method='field_bc', nhood='local', ) )


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
            
            #pytools.save_mpi_grid_to_disk(conf.outdir, lap, grid, conf) # MPI grid

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
                #print('simulation time    {:7d} ({:7.1f} omp)   {:5.1f}%'.format( int(lap), time, 100.0*lap/conf.Nt))
                print('sim time {:7d} ({:7.1f} omp) ({:7.2f} l0/c) {:5.2f}%'.format(int(lap), 
                                                                                        time, 
                                                                                        sch.antenna.tcur,
                                                                                        100.0*lap/conf.Nt))

            #--------------------------------------------------
            # deep IO
            if conf.full_interval > 0 and (lap % conf.full_interval == 0) and (lap > 0):
                pyfld.write_grids(grid, lap, conf.outdir + "/full_output/")
                pypic.write_particles(grid, lap, conf.outdir + "/full_output/")


            # restart IO (overwrites)
            if (lap % conf.restart == 0) and (lap > 0):

                # flip between two sets of files
                io_stat["deep_io_switch"] = 1 if io_stat["deep_io_switch"] == 0 else 0

                pyfld.write_grids(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                pypic.write_particles(
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

