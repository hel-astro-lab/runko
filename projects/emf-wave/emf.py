"""
Simple example of using runko module to simulate E and B fields in vacuum.
"""

import runko
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

if __name__ == "__main__":
    config = runko.Configuration(None)

    config.tile_partitioning = "hilbert_curve"
    config.Nt = 40
    config.Nx = 2
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 40
    config.NyMesh = 40
    config.NzMesh = 40
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"

    tile_grid = runko.TileGrid(config)

    # Initial conditions:

    k = 2 * np.pi * np.array((0, 0, 1)) / 10.

    def E0(x, y, z):
        r = np.array((x, y, z))
        return np.sin(np.dot(k, r)), 0, 0

    def B0(x, y, z):
        r = np.array((x, y, z))
        return 0, np.sin(np.dot(k, r)), 0

    J0 = lambda x, y, z: (0, 0, 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, config)
        tile.set_EBJ(E0, B0, J0)
        tile_grid.add_tile(tile, idx)

    # Simulation config:
    simulation = tile_grid.configure_simulation(config)
    Ecomm, Bcomm = runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B

    def lap_function(x):
        x.comm_external(Ecomm)
        x.comm_local(Ecomm)
        x.grid_push_half_b()
        x.grid_push_half_b()

        x.comm_external(Bcomm)
        x.comm_local(Bcomm)
        x.grid_push_e()

        if simulation.lap % 5 == 0:
            x.io_emf_snapshot()

        if simulation.lap % 5 == 0:
            simulation.log_timer_statistics()

    simulation.for_each_lap(lap_function)
