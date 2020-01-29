# -*- coding: utf-8 -*-

# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os
import h5py

# runko + auxiliary modules
import pycorgi # corgi c++ bindings
import pyrunko # runko c++ bindings
import pytools # runko python tools

# problem specific modules
from problem import Configuration_Problem as Configuration


# global simulation seed
np.random.seed(1)


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
    # ux, uy, uz, uu = boosted_maxwellian(delgam, gamma, direction=direction, dims=3)
    ux, uy, uz, uu = 0.0, 0.0, 0.0, 0.0

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

    for cid in grid.get_tile_ids():
        tile = grid.get_tile(cid)
        yee = tile.get_yee(0)

        ii, jj, kk = tile.index

        # insert values into Yee lattices; includes halos from -3 to n+3
        for n in range(-3, conf.NzMesh + 3):
            for m in range(-3, conf.NyMesh + 3):
                for l in range(-3, conf.NxMesh + 3):
                    # get global coordinates
                    iglob, jglob, kglob = pytools.pic.threeD.ind2loc(
                        (ii, jj, kk), (l, m, n), conf
                    )

                    yee.bx[l, m, n] = conf.binit * np.cos(btheta)
                    yee.by[l, m, n] = conf.binit * np.sin(btheta) * np.sin(bphi)
                    yee.bz[l, m, n] = conf.binit * np.sin(btheta) * np.cos(bphi)

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
    # setup grid
    grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

    # compute initial mpi ranks using Hilbert's curve partitioning
    pytools.balance_mpi_3D(grid)

    # load pic tiles into grid
    pytools.pic.threeD.load_tiles(grid, conf)

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
        prtcl_stat = pytools.pic.threeD.inject(
            grid, velocity_profile, density_profile, conf
        )
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
        pyrunko.fields.threeD.read_yee(grid, io_stat["read_lap"], io_stat["read_dir"])
        pyrunko.pic.threeD.read_prtcls(grid, io_stat["read_lap"], io_stat["read_dir"])

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    #static load balancing setup; communicate neighborhood info once

    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    # load virtual mpi halo tiles
    pytools.pic.threeD.load_virtual_tiles(grid, conf)


    # --------------------------------------------------
    # --------------------------------------------------
    # --------------------------------------------------
    # end of initialization
    
    timer.stop("init") 
    timer.stats("init") 


    # --------------------------------------------------
    # load physics solvers


    #pusher   = pypic.BorisPusher()
    pusher   = pypic.VayPusher()

    #fldprop  = pyfld.FDTD2()
    fldprop  = pyfld.FDTD4()

    fintp    = pypic.LinearInterpolator()
    currint  = pypic.ZigZag()
    flt      = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    #analyzer = pypic.Analyzator()

    #enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    fldprop.corr = 1.02


    # --------------------------------------------------
    # I/O objects

    # quick field snapshots
    debug_print(grid, "qwriter")
    qwriter  = pyfld.QuickWriter(
            conf.outdir, 
            conf.Nx, conf.NxMesh,
            conf.Ny, conf.NyMesh,
            conf.Nz, conf.NzMesh,
            conf.stride)

    # test particles
    debug_print(grid, "tpwriter")
    tpwriter = pypic.TestPrtclWriter(
            conf.outdir, 
            conf.Nx, conf.NxMesh,
            conf.Ny, conf.NyMesh,
            conf.Nz, conf.NzMesh,
            conf.ppc, len(grid.get_local_tiles()),
            conf.n_test_prtcls)

    # --------------------------------------------------
    #reflecting leftmost wall
    piston   = pyrunko.pic.threeD.Piston()

    # set piston wall speed (for standard reflection non-moving and so 0)
    piston.gammawall = conf.wallgamma
    piston.betawall = np.sqrt(1.-1./conf.wallgamma**2.)
    piston.walloc = 5.0 #leave 5 cell spacing between the wall for boundary conditions


    # --------------------------------------------------
    # sync e and b fields

    #mpi e
    grid.send_data(1) 
    grid.recv_data(1) 
    grid.wait_data(1) 

    #mpi b
    grid.send_data(2) 
    grid.recv_data(2) 
    grid.wait_data(2) 

    ################################################## 
    # simulation time step loop

    sys.stdout.flush()

    #simulation loop
    time = lap*(conf.cfl/conf.c_omp)
    for lap in range(lap, conf.Nt+1):















    # --------------------------------------------------
    # end of simulation

    timer.stop("total")
    timer.stats("total")
