# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import itertools
import pickle
import logging
import pathlib
from mpi4py import MPI
import pycorgi.threeD as pycorgi
from .simulation import Simulation
from .runko_logging import runko_logger, on_main_rank
from .balance_grid import balance_mpi, load_catepillar_track_mpi
from .auto_outdir import resolve_outdir


class TileGrid:
    """
    Represents the grid of runko tiles.
    Stores tiles in current locality and global information about all tiles.

    Initially tiles are distributed to different localities based on `tile_partitioning`.
    Possible values are: "hilbert_curve" and "catepillar_track"

    "catepillar_track" requires to be configured with `catepillar_track_length` (int).
    """

    def __init__(self, conf):

        required_vars = ["n_tiles",
                         "n_cells_per_tile",
                         "tile_partitioning"]

        for var in required_vars:
            if getattr(conf, var) is None:
                raise RuntimeError(f"Can not construct TileGrid without: {var}")

        self._Nx, self._Ny, self._Nz = conf.n_tiles
        self._NxMesh, self._NyMesh, self._NzMesh = conf.n_cells_per_tile
        self._xmin, self._ymin, self._zmin = (0, 0, 0)

        valid_tile_partitions = ["hilbert_curve", "catepillar_track"]

        if conf.tile_partitioning not in valid_tile_partitions:
            raise RuntimeError(f"invalid `tile_partitioning`: {conf.tile_partitioning}")

        if conf.tile_partitioning == "catepillar_track":
            if conf.catepillar_track_length is None:
                msg = "Catepillar track requested, but no `catepillar_track_length` given!"
                raise RuntimeError(msg)


        self._corgi_grid = pycorgi.Grid(self._Nx, self._Ny, self._Nz)

        xmax = self._xmin + self._Nx * self._NxMesh
        ymax = self._ymin + self._Ny * self._NyMesh
        zmax = self._zmin + self._Nz * self._NzMesh

        self._corgi_grid.set_grid_lims(self._xmin,
                                       xmax,
                                       self._ymin,
                                       ymax,
                                       self._zmin,
                                       zmax)

        legacy_conf = type("", (), dict(oneD=False,
                                        twoD=False,
                                        threeD=True,
                                        mpi_task_mode=False))

        if conf.tile_partitioning == "hilbert_curve":
            balance_mpi(self._corgi_grid, legacy_conf, do_print=False)
        elif conf.tile_partitioning == "catepillar_track":
            load_catepillar_track_mpi(self._corgi_grid,
                                              conf.catepillar_track_length,
                                              legacy_conf,
                                              do_print=False)
        else:
            raise RuntimeError("Due to previous checking this should not happend.")

        self._logger = runko_logger("TileGrid")
        self._logger.debug(f"TileGrid constructed with configuration: {conf}")


    def add_tile(self, tile, tile_grid_idx: (int, int, int)):
        """
        Adds tile to given tile grid index.
        """

        i, j, k = tile_grid_idx
        my_rank, ijk_rank = self._corgi_grid.rank(), self._corgi_grid.get_mpi_grid(i, j, k)
        if my_rank != ijk_rank:
            msg = f"rank {my_rank}: Trying to add tile to {idx} which belongs to different MPI rank."
            raise RuntimeError(msg)

        self._corgi_grid.add_tile(tile, tile_grid_idx)

        self._logger.debug(f"Added tile {tile} at {tile_grid_idx}.")


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

        index_space = itertools.product(range(self._corgi_grid.get_Nx()),
                                        range(self._corgi_grid.get_Ny()),
                                        range(self._corgi_grid.get_Nz()))

        for i, j, k in index_space:
            if self._corgi_grid.rank() == self._corgi_grid.get_mpi_grid(i, j, k):
                yield (i, j, k)


    def configure_simulation(self, config) -> Simulation:
        """
        Configures execution ready runko simulation based on the tile grid.
        """

        required_vars = ["n_laps"]

        for var in required_vars:
            if getattr(config, var) is None:
                raise RuntimeError(f"Can not configure simulation without: {var}")

        self._corgi_grid.analyze_boundaries()
        self._corgi_grid.send_tiles()
        self._corgi_grid.recv_tiles()
        MPI.COMM_WORLD.barrier()

        for vtile_id in self._corgi_grid.get_virtual_tiles():
            vtile = self._corgi_grid.get_tile(vtile_id)

            # If virtual tile is not from pycorgi,
            # assume that it is already initialized.
            if 'pycorgi' not in type(vtile).__module__:
                continue

            tile_type_candidate = None

            for i, j, k in vtile.nhood():
                nhood_tile = self._corgi_grid.get_tile(i, j, k)

                if not nhood_tile:
                    continue

                try:
                    nhood_tile_type = type(nhood_tile).canonical_type()
                except:
                    nhood_tile_type = type(nhood_tile)

                # Treat tiles from pycorgi as non-initialized.
                if 'pycorgi' not in nhood_tile_type.__module__ :

                    if tile_type_candidate and tile_type_candidate != nhood_tile_type:
                        msg = "Can not deduce virtual tile type.\n"
                        msg += f"It is next to {tile_type_candidate} and {nhood_tile_type}."
                        raise RuntimeError(msg)

                    tile_type_candidate = nhood_tile_type

            indices = vtile.index
            new_tile = None
            try:
                # Try to create a virtual tile specialization of the tile.
                new_tile = tile_type_candidate.virtual_tile_specialization()(indices, config)
            except TypeError:
                # There is no virtual tile specialization for the tile.
                new_tile = tile_type_candidate(indices, config)
            except Exception as e:
                raise e

            self._corgi_grid.add_tile(new_tile, indices)
            new_tile.load_metainfo(vtile.communication)

            vtile_msg = f"Constructed virtual tile at {vtile.index} "
            vtile_msg += f"with deduced type {tile_type_candidate}."
            self._logger.debug(vtile_msg)

        if config.verbose:
            self._logger.info(f"simulation configured with: {config.__dict__}")

        # Count particle species from config (q0/m0, q1/m1, ...)
        nspecies = 0
        for i in range(6):
            if getattr(config, f"q{i}") is not None and getattr(config, f"m{i}") is not None:
                nspecies += 1
            else:
                break

        stride = 1 if not config.io_grid_stride else config.io_grid_stride
        io_config = dict(stride=stride,
                         outdir=resolve_outdir(config),
                         nspecies=nspecies,
                         n_prtcls=config.n_sampled_prtcls if config.n_sampled_prtcls else 0,
                         laps_in_timer_statistics=getattr(config, 'io_n_laps_in_timer_stats', None),
                         spectra_nbins=getattr(config, 'io_n_spectra_bins', 200),
                         spectra_umin=getattr(config, 'io_spectra_umin', 1e-4),
                         spectra_umax=getattr(config, 'io_spectra_umax', 1e3),
                         spectra_stride=getattr(config, 'io_spectra_stride', None) or stride)

        pathlib.Path(io_config["outdir"]).mkdir(parents=True, exist_ok=True)

        if on_main_rank():
            pickled_conf_path = pathlib.Path(f"{io_config['outdir']}/config.pkl")
            config.ranks = MPI.COMM_WORLD.size
            with open(pickled_conf_path, "wb") as f:
                pickle.dump(config, f)

            if config._config_path:
                import shutil
                shutil.copy2(config._config_path, io_config["outdir"])

        return Simulation(self,
                          Simulation._im_not_user,
                          Nt=config.n_laps,
                          io_config=io_config,
                          verbose=config.verbose)
