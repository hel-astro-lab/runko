# -*- coding: utf-8 -*-

# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os

# runko + auxiliary modules
import pytools  # runko python tools


# problem specific modules
from problem import Configuration_Problem as Configuration

np.random.seed(1)  # global simulation seed


# Prtcl velocity (and location modulation inside cell)
#
# NOTE: Cell extents are from xloc to xloc + 1, i.e., dx = 1 in all dims.
#       Therefore, typically we use position as x0 + RUnif[0,1).
#
# NOTE: Default injector changes odd ispcs's loc to (ispcs-1)'s prtcl loc.
#       This means that positrons/ions are on top of electrons to guarantee
#       charge conservation (because zero charge density initially).
#
def velocity_profile(xloc, ispcs, conf):

    # electrons
    if ispcs == 0:
        delgam = conf.delgam  # * np.abs(conf.mi / conf.me) * conf.temp_ratio

    # positrons/ions/second species
    elif ispcs == 1:
        delgam = conf.delgam

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    zz = xloc[2] + np.random.rand(1)

    # velocity sampling
    gamma = conf.gamma
    direction = -1
    ux, uy, uz, uu = pytools.sample_boosted_maxwellian(
        delgam, gamma, direction=direction, dims=3
    )

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0


# number of prtcls of species 'ispcs' to be added to cell at location 'xloc'
#
# NOTE: Plasma frequency is adjusted based on conf.ppc (and prtcl charge conf.qe/qi
#       and mass conf.me/mi are computed accordingly) so modulate density in respect
#       to conf.ppc only.
#
def density_profile(xloc, ispcs, conf):

    # uniform plasma with default n_0 number density
    return conf.ppc


# Field initialization
def insert_em_fields(grid, conf):

    # into radians
    btheta = conf.btheta / 180.0 * np.pi
    bphi = conf.bphi / 180.0 * np.pi
    beta = conf.beta

    for tile in pytools.tiles_all(grid):
        yee = tile.get_yee(0)

        ii, jj, kk = tile.index if conf.threeD else (*tile.index, 0)

        # insert values into Yee lattices; includes halos from -3 to n+3
        for n in range(-3, conf.NzMesh + 3):
            for m in range(-3, conf.NyMesh + 3):
                for l in range(-3, conf.NxMesh + 3):
                    # get global coordinates
                    iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                    r = np.sqrt(iglob ** 2 + jglob ** 2 + kglob ** 2)

                    yee.bx[l, m, n] = conf.binit * np.cos(bphi)
                    yee.by[l, m, n] = conf.binit * np.sin(bphi) * np.sin(btheta)
                    yee.bz[l, m, n] = conf.binit * np.sin(bphi) * np.cos(btheta)

                    yee.ex[l, m, n] = 0.0
                    yee.ey[l, m, n] = -beta * yee.bz[l, m, n]
                    yee.ez[l, m, n] = beta * yee.by[l, m, n]
    return


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

    # 3D box peripherals
    if conf.threeD:
        slice_xy_writer = pyfld.FieldSliceWriter(
            conf.outdir,
            conf.Nx,
            conf.NxMesh,
            conf.Ny,
            conf.NyMesh,
            conf.Nz,
            conf.NzMesh,
            1,
            0,
            1,
        )
        slice_xz_writer = pyfld.FieldSliceWriter(
            conf.outdir,
            conf.Nx,
            conf.NxMesh,
            conf.Ny,
            conf.NyMesh,
            conf.Nz,
            conf.NzMesh,
            1,
            1,
            1,
        )
        slice_yz_writer = pyfld.FieldSliceWriter(
            conf.outdir,
            conf.Nx,
            conf.NxMesh,
            conf.Ny,
            conf.NyMesh,
            conf.Nz,
            conf.NzMesh,
            1,
            2,
            1,
        )

    # --------------------------------------------------
    # reflecting leftmost wall
    piston = pypic.Piston()

    # set piston wall speed (for standard reflector it is non-moving so gam = 0)
    piston.walloc = 5.0  # leave 5 cell spacing between the wall for boundary conditions

    if conf.wallgamma > 0.0:
        # moving wall
        piston.gammawall = conf.wallgamma
        piston.betawall = np.sqrt(1.0 - 1.0 / conf.wallgamma ** 2.0)
    else:
        # stationary wall
        piston.betawall = 0.0
        piston.gammawall = 1.0

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
            piston.field_bc(tile)
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

        # --------------------------------------------------
        # apply moving/reflecting walls
        t1 = timer.start_comp("walls")
        for tile in pytools.tiles_local(grid):
            piston.solve(tile)
        timer.stop_comp(t1)

        ##################################################
        # advance B half

        # --------------------------------------------------
        # push B half
        t1 = timer.start_comp("push_half_b2")
        for tile in pytools.tiles_all(grid):
            fldprop.push_half_b(tile)
            piston.field_bc(tile)
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
            piston.field_bc(tile)
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

        # clean current behind piston
        if conf.npasses > 0:
            for tile in pytools.tiles_local(grid):
                piston.field_bc(tile)

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

        # moving walls
        # piston.walloc += piston.betawall*conf.cfl

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
            # NOTE: do moms before other IOs to keep rho up-to-date
            mom_writer.write(grid, lap)  # pic distribution moments;
            fld_writer.write(grid, lap)  # quick field snapshots
            prtcl_writer.write(grid, lap)  # test particles

            # box peripheries
            if conf.threeD:
                slice_xy_writer.write(grid, lap)
                slice_xz_writer.write(grid, lap)
                slice_yz_writer.write(grid, lap)

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
