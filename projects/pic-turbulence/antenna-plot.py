# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Plot of current from oscillating Langevin antenna.
"""

import runko
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np


def make_config():
    config = runko.Configuration(None)

    config.io_outdir = "antenna-plot-output"
    config.n_tiles = [1, 1, 1]
    config.n_cells_per_tile = [10, 10, 10]
    config.cfl = 0.45
    config.n_laps = 1000
    config.field_propagator = "fdtd2"
    config.m0 = 1
    config.m1 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"
    config.current_filter = "binomial2"
    config.tile_partitioning = "hilbert_curve"

    return config


def make_antennas(config, N=1, sigma=10):
    modes = [(N, 0, N), (N, 0, -N), (0, N, N), (0, N, -N)]

    eddy_length = config.n_cells_per_tile[0] * config.n_tiles[0] / N
    light_crossing_time = eddy_length / config.cfl # in units of dt

    v_alfen = np.sqrt(sigma / (1 + sigma)) # in units of c
    w0 = 2 * np.pi * config.cfl * v_alfen / eddy_length # = k_z * v_A * dt

    # Note that in multi-node simulations, # it is important to give this `gen` argument
    # with numpy rng that has a defined seed # and every node uses the same seed.
    #
    # Otherwise, different parts of the modes will evolve differently.
    # Here it does not matter, because we only have one node.
    time_evolution = runko.sample_oscillating_langevin_antenna(size=config.n_laps,
                                                               characteristic_freq=0.8 * w0,
                                                               decorrelation_rate=-0.06 * w0)

    antenna_mode = runko.emf.threeD.antenna_mode
    to_antenna_modes = lambda mode: antenna_mode(A=(0, 0, 1), n=mode, lap_coeffs=time_evolution)
    return list(map(to_antenna_modes, modes))


if __name__ == "__main__":
    config = make_config()

    zero_field = lambda x, y, z: np.zeros_like(x)

    tile_grid = runko.TileGrid(config)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.emf.threeD.Tile(idx, config)

            for antenna in make_antennas(config):
                tile.register_antenna(antenna)

            tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(config)

    def zero_fields():
        # Note that this is slow way to set these to zero.
        # Normally depositi_current would implicitly zero out current
        # but we are not using it here, so we zero out the current manually.
        for tile in simulation.local_tiles():
            tile.batch_set_EBJ(zero_field, zero_field, zero_field,
                               zero_field, zero_field, zero_field,
                               zero_field, zero_field, zero_field)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set(xlabel="x", ylabel="y", zlabel="z")
    x, y, z = np.meshgrid(*map(range, config.n_cells_per_tile), indexing="ij")

    Q = ax.quiver(x, y, z, np.zeros_like(x), np.zeros_like(x), np.zeros_like(x))

    def update(frame):
        global Q
        Q.remove()

        zero_fields()
        simulation.for_one_lap(lambda x: x.grid_deposit_antenna_current())

        for tile in simulation.local_tiles():
            _, _, (Jx, Jy, Jz) = tile.get_EBJ()

            Q = ax.quiver(x, y, z, Jx, Jy, Jz)

        return Q,

    # Without init_func FuncAnimation calls update function twice with frame = 0, so frames = n_laps - 1.
    anim = animation.FuncAnimation(fig, update, frames=config.n_laps - 1, interval=50, repeat=False)
    plt.show()
