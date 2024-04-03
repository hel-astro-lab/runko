# -*- coding: utf-8 -*-

# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os

# runko + auxiliary modules
import pytools  # runko python tools

# problem specific modules
from problem import Configuration_Packets as Configuration

np.random.seed(1)  # global simulation seed


if __name__ == "__main__":

    # --------------------------------------------------
    # initial setup
    do_print = False
    if MPI.COMM_WORLD.Get_rank() == 0:
        do_print = True

    if do_print:
        print("Running ffe.py with {} MPI processes.".format(MPI.COMM_WORLD.Get_size()))

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
        import pyrunko.ffe.threeD as pyffe  # runkko ffe c++ bindings
        import pyrunko.emf.threeD as pyfld  # runko fld c++ bindings

    elif conf.twoD:
        # 2D modules
        import pycorgi.twoD as pycorgi  # corgi ++ bindings
        import pyrunko.ffe.twoD as pyffe  # runkko ffe c++ bindings
        import pyrunko.emf.twoD as pyfld  # runko fld c++ bindings


    # --------------------------------------------------
    # problem setup
    if True:
        from packet_setup import insert_em_fields # torsional harmonic envelope 


    # --------------------------------------------------
    # setup grid
    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

    # compute initial mpi ranks using Hilbert's curve partitioning
    pytools.balance_mpi(grid, conf)

    # load ffe tiles into grid
    pytools.ffe.load_tiles(grid, conf)

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
        insert_em_fields(grid, conf)
    else:
        if do_print:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        # read restart files
        pyfld.read_grids(grid, io_stat["read_lap"], io_stat["read_dir"])

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    # static load balancing setup; communicate neighborhood info once

    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    # load virtual mpi halo tiles
    pytools.ffe.load_virtual_tiles(grid, conf)

    # --------------------------------------------------
    # load physics solvers

    # reduced 2nd order FFE algorithm

    #algo = pyffe.rFFE2(conf.NxMesh, conf.NyMesh, conf.NzMesh)
    #algo = pyffe.rFFE2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    #algo = pyffe.FFE2(conf.NxMesh, conf.NyMesh, conf.NzMesh)
    algo = pyffe.FFE4(conf.NxMesh, conf.NyMesh, conf.NzMesh)

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

        # initialize Y^n-1 = Y
        t1 = timer.start_comp("copy_eb")
        for tile in pytools.tiles_all(grid):
            tile.copy_eb()
        timer.stop_comp(t1)

        ###################################################
        # rk substeps
        rks = 0 #substep counter

        # RK1
        # rk_coeffs = [(1.0,   0.0,   1.0,   1.0),]

        # RK2
        # rk_coeffs = [(1.0,   0.0,   1.0,   1.0),
        #             (0.5,   0.5,   0.5,   1.0),]

        # RK3
        rk_coeffs = [
            (1.0, 0.0, 1.0, 1.0),
            (0.75, 0.25, 0.25, 0.5),
            (1 / 3, 2 / 3, 2 / 3, 1.0),
        ]

        for (rk_c1, rk_c2, rk_c3, rk_dt) in rk_coeffs:
            rks += 1 # RK substep

            #--------------------------------------------------
            # comm E/B
            t1 = timer.start_comp("mpi_loop")
            grid.send_data(1)
            grid.recv_data(1)

            grid.send_data(2)
            grid.recv_data(2)

            grid.wait_data(1)
            grid.wait_data(2)
            timer.stop_comp(t1)

            ## update boundaries
            t1 = timer.start_comp("upd_bc_loop")
            for tile in pytools.tiles_local(grid):
                tile.update_boundaries(grid)
            timer.stop_comp(t1)

            #--------------------------------------------------

            # rho = div E
            t1 = timer.start_comp("comp_rho")
            for tile in pytools.tiles_local(grid):
                algo.comp_rho(tile)
            timer.stop_comp(t1)

            # dE = dt * curl B
            # dB = dt * curl E
            t1 = timer.start_comp("push_eb")
            for tile in pytools.tiles_local(grid):
                algo.push_eb(tile)
            timer.stop_comp(t1)

            # drift current j_perp
            # dE -= dt*j_perp
            t1 = timer.start_comp("add_jperp")
            for tile in pytools.tiles_local(grid):
                algo.add_jperp(tile)
            timer.stop_comp(t1)

            # parallel current j_par
            # dE -= dt*j_par
            t1 = timer.start_comp("add_jpar")
            for tile in pytools.tiles_local(grid):
                algo.add_jpar(tile)
            timer.stop_comp(t1)

            # diffusion
            if True:
                algo.eta = 1.0e-3

                # dE += eta*dt*nabla^2 E
                t1 = timer.start_comp("diffuse")
                for tile in pytools.tiles_local(grid):
                    algo.add_diffusion(tile)
                timer.stop_comp(t1)

            # update fields according to RK scheme
            # Y^n+1 = c1 * Y^n-1 + c2 * Y^n + c3 * dY
            t1 = timer.start_comp("update_eb")
            for tile in pytools.tiles_local(grid):
                tile.rk3_update(     rk_c1, rk_c2, rk_c3)
                #algo.update_eb(tile, rk_c1, rk_c2, rk_c3)
            timer.stop_comp(t1)

            # jpar
            #if True:
            #if rks == 3:
            if False:
                # comm e & b (NOTE: need mpi here because Y^n is updated above)
                t1 = timer.start_comp("mpi_jpar")
                grid.send_data(1)
                grid.recv_data(1)

                grid.send_data(2)
                grid.recv_data(2)

                grid.wait_data(1)
                grid.wait_data(2)
                timer.stop_comp(t1)

                ## update boundaries
                t1 = timer.start_comp("upd_bc_jpar")
                for tile in pytools.tiles_local(grid):
                    tile.update_boundaries(grid)
                timer.stop_comp(t1)

                # parallel current j_par
                # dE -= j_par
                t1 = timer.start_comp("remove_jpar")
                for tile in pytools.tiles_local(grid):
                    algo.remove_jpar(tile)
                timer.stop_comp(t1)

            # eGTb
            #if False:
            #if rks == 3:
            if True:
                # comm E/B comm e & b (NOTE: need mpi here because Y^n is updated above)
                t1 = timer.start_comp("mpi_egtb")
                grid.send_data(1)
                grid.recv_data(1)

                grid.send_data(2)
                grid.recv_data(2)

                grid.wait_data(1)
                grid.wait_data(2)
                timer.stop_comp(t1)

                ## update boundaries
                t1 = timer.start_comp("upd_bc_egtb")
                for tile in pytools.tiles_local(grid):
                    tile.update_boundaries(grid)
                timer.stop_comp(t1)

                # enforce E < B
                # dE = dE_lim
                t1 = timer.start_comp("limit_e")
                for tile in pytools.tiles_local(grid):
                    algo.limit_e(tile)
                timer.stop_comp(t1)

            ##################################################
            # update field halos

            # comm e & b
            t1 = timer.start_comp("mpi_eb2")
            grid.send_data(1)
            grid.recv_data(1)

            grid.send_data(2)
            grid.recv_data(2)

            grid.wait_data(1)
            grid.wait_data(2)
            timer.stop_comp(t1)

            # update boundaries
            t1 = timer.start_comp("upd_bc2")
            for tile in pytools.tiles_local(grid):
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

            # shallow IO
            fld_writer.write(grid, lap)  # quick field snapshots

            # deep IO
            if (
                conf.full_interval != -1
                and (lap % conf.full_interval == 0)
                and (lap > 0)
            ):
                pyfld.write_grids(grid, lap, conf.outdir + "/full_output/")

            # restart IO (overwrites)
            if (lap % conf.restart == 0) and (lap > 0):
                # flip between two sets of files
                io_stat["deep_io_switch"] = 1 if io_stat["deep_io_switch"] == 0 else 0

                pyfld.write_grids(
                    grid, io_stat["deep_io_switch"], conf.outdir + "/restart/"
                )

                # if successful adjust info file
                MPI.COMM_WORLD.barrier()
                if grid.rank() == 0:
                    with open(conf.outdir + "/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(lap, io_stat["deep_io_switch"]))

            timer.stop("io")

            timer.stats("io")
            timer.start("step")  # refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        # MPI.COMM_WORLD.barrier()
        time += conf.cfl / conf.c_omp
    # end of loop

    timer.stop("total")
    timer.stats("total")
