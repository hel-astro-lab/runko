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

    # print(pytools.__name__)
    # print(dir(pytools))
    # print(pytools.pic.__name__)
    # print(dir(pytools.pic))

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
    xmin = 0.0
    ymin = 0.0
    zmin = 0.0
    xmax = conf.Nx * conf.NxMesh
    ymax = conf.Ny * conf.NyMesh
    zmax = conf.Nz * conf.NzMesh

    grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(xmin, xmax, ymin, ymax, zmin, zmax)

    # init.loadMpi3D(grid)

    pytools.pic.threeD.load_tiles(grid, conf)

    # --------------------------------------------------

    timer.stop("total")
    timer.stats("total")
