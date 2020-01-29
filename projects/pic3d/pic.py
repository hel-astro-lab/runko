# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os
import h5py

# runko + auxiliary modules
import pycorgi
import pyrunko
import pytools

# problem specific modules
from problem import Configuration_Problem as Configuration


# global simulation seed
np.random.seed(1)


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

    # create output folders
    if grid.master():
        pytools.create_output_folders(conf)

    # --------------------------------------------------
    # setup grid
    xmin = 0.0
    ymin = 0.0
    zmin = 0.0
    xmax = conf.Nx * conf.NxMesh
    ymax = conf.Ny * conf.NyMesh
    zmax = conf.Nz * conf.NzMesh

    grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(xmin, xmax, ymin, ymax, zmin, zmax)

    # compute initial mpi ranks using Hilbert's curve partitioning
    pytools.balance_mpi_3D(grid)

    # load pic tiles into grid
    pytools.pic.threeD.load_tiles(grid, conf)

    # --------------------------------------------------
    # simulation restart

    # get current restart file status
    io_stat = pytools.check_for_restart(conf)

    # no restart file; initialize simulation
    if io_stat["do_initialization"]:
        if do_print:
            print("initializing simulation...")

        lap = 0
        np.random.seed(1)
        inject(grid, filler, conf)  # injecting plasma particles
        insert_em(grid, conf)
    else:
        if do_print:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        pyrunko.fields.threeD.read_yee(grid, io_stat['read_lap'], io_stat['read_dir'])
        pyrunko.pic.threeD.read_prtcls(grid, io_stat['read_lap'], io_stat['read_dir'])

        lap = io_stat["lap"] + 1  # step one step ahead





    timer.stop("total")
    timer.stats("total")
