# -*- coding: utf-8 -*-

# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os

# runko + auxiliary modules
import pytools  # runko python tools

# problem specific modules
from setup_tests import Configuration_Test as Configuration

np.random.seed(1)  # global simulation seed



if __name__ == "__main__":

    # --------------------------------------------------
    # initial setup
    do_print = False
    if MPI.COMM_WORLD.Get_rank() == 0:
        do_print = True

    if do_print:
        print("Running pic.py with {} MPI processes.".format(MPI.COMM_WORLD.Get_size()))

    # --------------------------------------------------
    # Timer for profiling
    timer = pytools.Timer()

    timer.start("total")
    timer.start("init")
    timer.do_print = do_print

    # --------------------------------------------------
    # parse command line arguments
    args = pytools.parse_args()

    # create conf object with simulation parameters based on them
    conf = Configuration(args.conf_filename, do_print=do_print)

    # --------------------------------------------------
    # load runko

    if conf.threeD:
        # 3D modules
        import pycorgi.threeD as pycorgi  # corgi ++ bindings
        import pyrunko.pic.threeD as pypic  # runko pic c++ bindings
        import pyrunko.fields.threeD as pyfld  # runko fld c++ bindings

    elif conf.twoD:
        # 2D modules
        import pycorgi.twoD as pycorgi  # corgi ++ bindings
        import pyrunko.pic.twoD as pypic  # runko pic c++ bindings
        import pyrunko.fields.twoD as pyfld  # runko fld c++ bindings


    # --------------------------------------------------
    # load problem setup
    
    # Weibel
    if True:
        from weibel import velocity_profile
        from weibel import density_profile
        from weibel import insert_em_fields

    # Langmuir
    if False:
        from langmuir import velocity_profile
        from langmuir import density_profile
        from langmuir import insert_em_fields


    # --------------------------------------------------
    # setup grid
    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

    # compute initial mpi ranks using Hilbert's curve partitioning
    pytools.balance_mpi(grid, conf)

    # load pic tiles into grid
    pytools.pic.load_tiles(grid, conf)

    # --------------------------------------------------
    # simulation restart

    # create output folders
    if grid.master():
        pytools.create_output_folders(conf)

    # get current restart file status
    io_stat = pytools.check_for_restart(conf)

    # no restart file; initialize simulation
    if io_stat["do_initialization"]:
        if do_print:
            print("initializing simulation...")
        lap = 0

        np.random.seed(1)  # sync rnd generator seed for different mpi ranks

        # injecting plasma particles
        prtcl_stat = pytools.pic.inject(grid, velocity_profile, density_profile, conf)
        if do_print:
            print("injected:")
            print("    e- prtcls: {}".format(prtcl_stat[0]))
            print("    e+ prtcls: {}".format(prtcl_stat[1]))

        # inserting em grid
        insert_em_fields(grid, conf)

    else:
        if do_print:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        # read restart files
        pyfld.read_yee(grid, io_stat["read_lap"], io_stat["read_dir"])
        pypic.read_particles(grid, io_stat["read_lap"], io_stat["read_dir"])

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    # static load balancing setup; communicate neighborhood info once

    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    # load virtual mpi halo tiles
    pytools.pic.load_virtual_tiles(grid, conf)

    # --------------------------------------------------
    # load physics solvers

    # pusher   = pypic.BorisPusher()
    pusher = pypic.VayPusher()

    # fldprop  = pyfld.FDTD2()
    fldprop = pyfld.FDTD4()

    fintp = pypic.LinearInterpolator()
    currint = pypic.ZigZag()
    flt = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    # analyzer = pypic.Analyzator()

    # enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    fldprop.corr = 1.02

    # --------------------------------------------------
    # I/O objects

    # quick field snapshots
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

    # test particles
    prtcl_writer = pypic.TestPrtclWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.ppc,
        len(grid.get_local_tiles()),
        conf.n_test_prtcls,
    )

    mom_writer = pypic.PicMomentsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride,
    )


    # --------------------------------------------------
    # --------------------------------------------------
    # --------------------------------------------------
    # end of initialization

    timer.stop("init")
    timer.stats("init")
    # timer.verbose = 1  # 0 normal; 1 - debug mode

    # --------------------------------------------------
    # sync e and b fields

    # mpi e
    grid.send_data(1)
    grid.recv_data(1)
    grid.wait_data(1)

    # mpi b
    grid.send_data(2)
    grid.recv_data(2)
    grid.wait_data(2)

    for tile in pytools.tiles_all(grid):
        tile.update_boundaries(grid)

    ##################################################
    # simulation time step loop

    sys.stdout.flush()

    # simulation loop
    time = lap * (conf.cfl / conf.c_omp)
    for lap in range(lap, conf.Nt + 1):

        # --------------------------------------------------
        # push B half
        t1 = timer.start_comp("push_half_b1")
        for tile in pytools.tiles_all(grid):
            fldprop.push_half_b(tile)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # comm B
        t1 = timer.start_comp("mpi_b2")
        grid.send_data(2)
        grid.recv_data(2)
        grid.wait_data(2)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # update boundaries
        t1 = timer.start_comp("upd_bc1")
        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)
        timer.stop_comp(t1)

        ##################################################
        # move particles (only locals tiles)

        # --------------------------------------------------
        # interpolate fields
        t1 = timer.start_comp("interp_em")
        for tile in pytools.tiles_local(grid):
            fintp.solve(tile)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # push particles in x and u
        t1 = timer.start_comp("push")
        for tile in pytools.tiles_local(grid):
            pusher.solve(tile)
        timer.stop_comp(t1)

        ##################################################
        # advance B half

        # --------------------------------------------------
        # push B half
        t1 = timer.start_comp("push_half_b2")
        for tile in pytools.tiles_all(grid):
            fldprop.push_half_b(tile)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # comm B
        t1 = timer.start_comp("mpi_e1")
        grid.send_data(1)
        grid.recv_data(1)
        grid.wait_data(1)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # update boundaries
        t1 = timer.start_comp("upd_bc2")
        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)
        timer.stop_comp(t1)

        ##################################################
        # advance E

        # --------------------------------------------------
        # push E
        t1 = timer.start_comp("push_e")
        for tile in pytools.tiles_all(grid):
            fldprop.push_e(tile)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # current calculation; charge conserving current deposition
        t1 = timer.start_comp("comp_curr")
        for tile in pytools.tiles_local(grid):
            currint.solve(tile)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # clear virtual current arrays for boundary addition after mpi
        t1 = timer.start_comp("clear_vir_cur")
        for tile in pytools.tiles_virtual(grid):
            tile.clear_current()
        timer.stop_comp(t1)

        # --------------------------------------------------
        # mpi send currents
        t1 = timer.start_comp("mpi_cur")
        grid.send_data(0)
        grid.recv_data(0)
        grid.wait_data(0)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # exchange currents
        t1 = timer.start_comp("cur_exchange")
        for tile in pytools.tiles_all(grid):
            tile.exchange_currents(grid)
        timer.stop_comp(t1)

        ##################################################
        # particle communication (only local/boundary tiles)

        # --------------------------------------------------
        # local particle exchange (independent)
        t1 = timer.start_comp("check_outg_prtcls")
        for tile in pytools.tiles_local(grid):
            tile.check_outgoing_particles()
        timer.stop_comp("check_outg_prtcls")

        # --------------------------------------------------
        # global mpi exchange (independent)
        t1 = timer.start_comp("pack_outg_prtcls")
        for tile in pytools.tiles_boundary(grid):
            tile.pack_outgoing_particles()
        timer.stop_comp(t1)

        # --------------------------------------------------
        # MPI global particle exchange
        # transfer primary and extra data
        t1 = timer.start_comp("mpi_prtcls")
        grid.send_data(3)
        grid.recv_data(3)
        grid.wait_data(3)

        # orig just after send3
        grid.send_data(4)
        grid.recv_data(4)
        grid.wait_data(4)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # global unpacking (independent)
        t1 = timer.start_comp("unpack_vir_prtcls")
        for tile in pytools.tiles_virtual(grid):
            tile.unpack_incoming_particles()
            tile.check_outgoing_particles()
        timer.stop_comp(t1)

        # --------------------------------------------------
        # transfer local + global
        t1 = timer.start_comp("get_inc_prtcls")
        for tile in pytools.tiles_local(grid):
            tile.get_incoming_particles(grid)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # delete local transferred particles
        t1 = timer.start_comp("del_trnsfrd_prtcls")
        for tile in pytools.tiles_local(grid):
            tile.delete_transferred_particles()
        timer.stop_comp(t1)

        # --------------------------------------------------
        # delete all virtual particles (because new prtcls will come)
        t1 = timer.start_comp("del_vir_prtcls")
        for tile in pytools.tiles_virtual(grid):
            tile.delete_all_particles()
        timer.stop_comp(t1)

        ##################################################
        # filter
        timer.start_comp("filter")

        # sweep over npasses times
        for fj in range(conf.npasses):

            # update global neighbors (mpi)
            grid.send_data(0)
            grid.recv_data(0)
            grid.wait_data(0)

            # get halo boundaries and filter
            for tile in pytools.tiles_local(grid):
                tile.update_boundaries(grid)
            for tile in pytools.tiles_local(grid):
                flt.solve(tile)

            MPI.COMM_WORLD.barrier()  # sync everybody

        # --------------------------------------------------
        timer.stop_comp("filter")

        # --------------------------------------------------
        # add current to E
        t1 = timer.start_comp("add_cur")
        for tile in pytools.tiles_all(grid):
            tile.deposit_current()
        timer.stop_comp(t1)

        # comm E
        t1 = timer.start_comp("mpi_e2")
        grid.send_data(1)
        grid.recv_data(1)
        grid.wait_data(1)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # comm B
        t1 = timer.start_comp("mpi_b1")
        grid.send_data(2)
        grid.recv_data(2)
        grid.wait_data(2)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # update boundaries
        t1 = timer.start_comp("upd_bc0")
        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)
        timer.stop_comp(t1)

        ##################################################
        # data reduction and I/O

        timer.lap("step")
        if lap % conf.interval == 0:
            if do_print:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time))

            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()

            # analyze (independent)
            timer.start("io")

            # shrink particle arrays
            for tile in pytools.tiles_all(grid):
                tile.shrink_to_fit_all_particles()

            # analyze
            # for tile in pytools.tiles_local(grid): analyzer.analyze2d(tile)

            # barrier for quick writers
            MPI.COMM_WORLD.barrier()

            # shallow IO
            mom_writer.write(grid, lap)  # prtcl moments from vel distribution
            fld_writer.write(grid, lap)  # quick field snapshots
            prtcl_writer.write(grid, lap)  # test particles

            # deep IO
            if conf.full_interval > 0 and (lap % conf.full_interval == 0) and (lap > 0):
                pyfld.write_yee(grid, lap, conf.outdir + "/full_output/")
                pypic.write_particles(grid, lap, conf.outdir + "/full_output/")
                # pypic.write_analysis(grid, lap, conf.outdir + "/full_output/")

            # restart IO (overwrites)
            if (lap % conf.restart == 0) and (lap > 0):

                # flip between two sets of files
                io_stat["deep_io_switch"] = 1 if io_stat["deep_io_switch"] == 0 else 0

                pyfld.write_yee(
                    grid, io_stat["deep_io_switch"], conf.outdir + "/restart/"
                )
                pypic.write_particles(
                    grid, io_stat["deep_io_switch"], conf.outdir + "/restart/"
                )

                # if successful adjust info file
                MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write
                if grid.rank() == 0:
                    with open(conf.outdir + "/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(lap, io_stat["deep_io_switch"]))

            timer.stop("io")

            timer.stats("io")
            timer.start("step")  # refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        # next step
        time += conf.cfl / conf.c_omp
    # end of loop

    # --------------------------------------------------
    # end of simulation

    timer.stop("total")
    timer.stats("total")
