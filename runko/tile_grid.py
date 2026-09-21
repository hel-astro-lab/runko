# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import itertools
import pickle
import logging
import pathlib
import pycorgi.threeD as pycorgi
from .simulation import Simulation
from .runko_logging import runko_logger, on_main_rank
from .balance_grid import balance_mpi, load_catepillar_track_mpi
from .auto_outdir import resolve_outdir


class TileGrid:
    """
    Represents the 3D grid of runko tiles.
    Stores tiles in current locality and global information about all tiles.

    Initially tiles are distributed to different localities based on `tile_partitioning`.
    Possible values are: "hilbert_curve" and "catepillar_track"

    "catepillar_track" requires to be configured with `catepillar_track_length` (int).
    """

    def __init__(self, conf):

        required_vars = ["n_tiles",
                         "tile_partitioning"]

        for var in required_vars:
            if getattr(conf, var) is None:
                raise RuntimeError(f"Can not construct TileGrid without: {var}")

        self._initial_tiles = []

        self._Nx, self._Ny, self._Nz = conf.n_tiles
        self._NxMesh, self._NyMesh, self._NzMesh = conf.n_cells_per_tile

        valid_tile_partitions = ["hilbert_curve", "catepillar_track"]

        if conf.tile_partitioning not in valid_tile_partitions:
            raise RuntimeError(f"invalid `tile_partitioning`: {conf.tile_partitioning}")

        self._local_indices = []

        from mpi4py import MPI
        self._my_rank = MPI.COMM_WORLD.Get_rank()
        self._world_size = MPI.COMM_WORLD.Get_size()

        total_n_tiles = self._Nx * self._Ny * self._Nz

        tiles_per_rank = [total_n_tiles // self._world_size] * self._world_size

        # Add left over tiles evenly between the ranks.
        i = 0
        while sum(tiles_per_rank) != total_n_tiles:
            tiles_per_rank[i] += 1
            i = (i + 1) % self._world_size

        skipped_tiles = sum(tiles_per_rank[:self._my_rank])

        if conf.tile_partitioning == "hilbert_curve":
            is_power_of_two = lambda n: n > 0 and (n & (n - 1)) == 0
            if not is_power_of_two(self._Nx) or not is_power_of_two(self._Ny) or not is_power_of_two(self._Nz):
                raise "Hilbert curve tile partition requires n_tiles to be powers of two!"

            from .hilbert import Hilbert3D
            H = Hilbert3D(self._Nx, self._Ny, self._Nz)

            for n in range(skipped_tiles, skipped_tiles + tiles_per_rank[self._my_rank]):
                self._local_indices.append(H.inv(n))

        elif conf.tile_partitioning == "catepillar_track":
            index_space = list(itertools.product(range(self._Nx), range(self._Ny), range(self._Nz)))
            for idx in index_space[skipped_tiles:][:tiles_per_rank[self._my_rank]]:
                self._local_indices.append(idx)
        else:
            raise RuntimeError("Due to previous checking this should not happend.")

        self._logger = runko_logger("TileGrid")
        self._logger.debug(f"TileGrid constructed with configuration: {conf}")


    def add_tile(self, tile, tile_grid_idx: (int, int, int)):
        """
        Adds tile to given tile grid index.
        """

        self._initial_tiles.append(tile)


    def initialized_from_restart_file(self) -> bool:
        """
        Has the grid been initialized fully
        from the restart files specified in configuration.
        """

        return False


    def local_tile_indices(self):
        """
        Returns iterable which goes through all indices (i, j, k)
        corresponding to a local tile locations.
        """

        return self._local_indices


    def configure_simulation(self, config) -> Simulation:
        """
        Configures execution ready runko simulation based on the tile grid.
        """

        required_vars = ["n_laps"]

        for var in required_vars:
            if getattr(config, var) is None:
                raise RuntimeError(f"Can not configure simulation without: {var}")

        if config.verbose:
            self._logger.info(f"simulation configured with: {config.__dict__}")

        # Count particle species from config (q0/m0, q1/m1, ...)
        nspecies = 0
        for i in range(1000):
            if getattr(config, f"q{i}") is not None and getattr(config, f"m{i}") is not None:
                nspecies += 1
            else:
                break


        stride = 1 if not config.io_grid_stride else config.io_grid_stride
        io_config = dict(stride=stride,
                         outdir=resolve_outdir(config),
                         nspecies=nspecies,
                         n_prtcls=config.io_n_sampled_prtcls if config.io_n_sampled_prtcls else 0,
                         laps_in_timer_statistics=config.io_n_laps_in_timer_stats,
                         spectra_nbins=config.io_n_spectra_bins or 200,
                         spectra_umin=config.io_spectra_umin or 1e-4,
                         spectra_umax=config.io_spectra_umax or 1e3,
                         spectra_stride=config.io_spectra_stride or stride)

        pathlib.Path(io_config["outdir"]).mkdir(parents=True, exist_ok=True)

        if on_main_rank():
            pickled_conf_path = pathlib.Path(f"{io_config['outdir']}/config.pkl")
            config.ranks = self._world_size
            with open(pickled_conf_path, "wb") as f:
                pickle.dump(config, f)

            if config._config_path:
                import shutil
                shutil.copy2(config._config_path, io_config["outdir"])

        return Simulation(self._initial_tiles,
                          Simulation._im_not_user,
                          config=config,
                          Nt=config.n_laps,
                          io_config=io_config,
                          verbose=config.verbose)
