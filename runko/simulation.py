from .method_wrapper import MethodWrapper
from .runko_logging import runko_logger
from .runko_timer import Timer, timer_statistics

from runko_cpp_bindings.tools import _virtual_tile_sync_handshake_mode, comm_mode
from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter as FieldsWriter, _write_average_B_energy_density, _write_average_E_energy_density
from runko_cpp_bindings.emf.threeD import MpiioParticlesWriter as ParticlesWriter
from runko_cpp_bindings.emf.threeD import MpiioSpectraWriter as SpectraWriter
from runko_cpp_bindings.pic.threeD import _write_average_kinetic_energy

import logging
import time
import numpy as np
import pathlib
from mpi4py import MPI
from .ram_usage import get_rss_kB, get_gpu_mem_kB


class Simulation:
    """
    Handle to a configured runko simulation.
    """

    _im_not_user = 42

    def __init__(self, tile_grid, prevent_user_init: int, **kwargs):
        """
        Construct runko simulation from tile grid.

        `prevent_user_init` is to make sure users don't create
        objects of this calss on their own.
        """

        if prevent_user_init != Simulation._im_not_user:
            raise RuntimeError("Don't instantiate this class directly!")

        self._tile_grid = tile_grid
        self._lap = 0

        self._last_lap = kwargs['Nt']

        self._io_config = kwargs['io_config']
        self._emf_writer = None
        self._prtcl_writers = {}  # species -> ParticlesWriter, lazily constructed
        self._spectra_writer = None

        self._lap_timers = []
        self._lap_wall_times = []

        self.verbose_laps = kwargs.get('verbose', True)

        self._logger = runko_logger("Simulation")

        rank = MPI.COMM_WORLD.Get_rank()
        self._ram_file = pathlib.Path(f"{self._io_config['outdir']}/ram-usage/{rank}.csv")
        self._ram_file.parent.mkdir(exist_ok=True, parents=True)
        self._ram_file.write_text("lap,ram usage [kB]\n")

        if get_gpu_mem_kB() is not None:
            self._gpu_ram_file = pathlib.Path(f"{self._io_config['outdir']}/ram-usage/{rank}_gpu.csv")
            self._gpu_ram_file.write_text("lap,gpu mem usage [kB]\n")
        else:
            self._gpu_ram_file = None

        self._io_config['kinetic_energy_path'] = self._io_config["outdir"] + "/average_kinetic_energy.txt"
        self._io_config['average_B_energy_density_path'] = self._io_config["outdir"] + "/average_B_energy_density.txt"
        self._io_config['average_E_energy_density_path'] = self._io_config["outdir"] + "/average_E_energy_density.txt"

        # Reset txt output files.
        for name in ('kinetic_energy_path', 'average_B_energy_density_path', 'average_E_energy_density_path'):
            pathlib.Path(self._io_config[name]).unlink(missing_ok=True)

        ctor_msg = "Simulation constructed with:\n"
        ctor_msg += f"\tNt = {kwargs['Nt']}\n"
        ctor_msg += f"\tio config: {self._io_config}"
        self._logger.debug(ctor_msg)


    def _ensure_constructed_emf_writer(self):
        if self._emf_writer:
            return

        self._emf_writer = FieldsWriter(self._io_config["outdir"],
                                        self._tile_grid._Nx,
                                        self._tile_grid._NxMesh,
                                        self._tile_grid._Ny,
                                        self._tile_grid._NyMesh,
                                        self._tile_grid._Nz,
                                        self._tile_grid._NzMesh,
                                        self._io_config["stride"],
                                        self._io_config.get("nspecies", 2))

        self._logger.debug("FieldsWriter constructed.")


    def _ensure_constructed_prtcl_writers(self):
        if self._prtcl_writers:
            return

        n_prtcls = self._io_config.get("n_prtcls", 0)
        nspecies = self._io_config.get("nspecies", 2)
        outdir   = self._io_config["outdir"]

        for sp in range(nspecies):
            self._prtcl_writers[sp] = ParticlesWriter(outdir, n_prtcls, sp)

        self._logger.debug(f"ParticlesWriters constructed for {nspecies} species, n_prtcls={n_prtcls}.")


    def _ensure_constructed_spectra_writer(self):
        if self._spectra_writer:
            return

        io = self._io_config
        self._spectra_writer = SpectraWriter(
            io["outdir"],
            self._tile_grid._Nx,
            self._tile_grid._NxMesh,
            self._tile_grid._Ny,
            self._tile_grid._NyMesh,
            self._tile_grid._Nz,
            self._tile_grid._NzMesh,
            io.get("spectra_stride", io["stride"]),
            io["spectra_nbins"],
            io["spectra_umin"],
            io["spectra_umax"],
            io.get("nspecies", 2))

        self._logger.debug("SpectraWriter constructed.")


    def virtual_tiles(self):
        """
        Return iterable which goes through all virtual tiles in current rank.
        """

        for vtile_id in self._tile_grid._corgi_grid.get_virtual_tiles():
            yield self._tile_grid._corgi_grid.get_tile(vtile_id)


    def local_tiles(self):
        """
        Return iterable which goes through all local tiles in current rank.
        """

        for tile_id in self._tile_grid._corgi_grid.get_local_tiles():
            yield self._tile_grid._corgi_grid.get_tile(tile_id)


    def _boundary_tiles(self):
        """
        Return iterable which goes through all boundary tiles in current rank.
        """

        for tile_id in self._tile_grid._corgi_grid.get_boundary_tiles():
            yield self._tile_grid._corgi_grid.get_tile(tile_id)


    @property
    def lap(self):
        """Current simulation lap."""
        return self._lap


    # Original Linux-only implementation (requires /proc/{pid}/smaps):
    # def write_mem_usage(self):
    #     import subprocess, os
    #     system_mem = f"awk '/^Pss:/ {{pss+=$2}} END {{print {self.lap} \",\" pss}}' < /proc/{os.getpid()}/smaps >> {self._ram_file}"
    #     subprocess.run(system_mem, shell=True)

    def _write_ram_usage(self):
        """Append current RAM usage (RSS in kB) for this lap to the per-rank CSV."""
        with open(self._ram_file, "a") as f:
            f.write(f"{self.lap},{get_rss_kB()}\n")
        if self._gpu_ram_file is not None:
            with open(self._gpu_ram_file, "a") as f:
                f.write(f"{self.lap},{get_gpu_mem_kB()}\n")


    def _execute_lap_function(self, lap_function, disable_timing=False):
        """
        Execute a given lap function.

        FIXME: Define lap function.
        FIXME: Parallelize the loop.
        """

        lap_wall_time_begin = time.time()
        lap_timer = Timer() if not disable_timing else None

        def get_name(method, kwargs):
            if 'name' in kwargs:
                return kwargs['name']
            else:
                return method

        def pre(method, *vargs, **kwargs):
            name = get_name(method, kwargs)
            lap_timer.start(name)

        def post(method, *vargs, **kwargs):
            name = get_name(method, kwargs)
            lap_timer.stop(name)

        def action(method: str, *vargs):
            if method.startswith("prtcl_"):
                method_mapper = { "prtcl_push" : "push_particles",
                                  "prtcl_sort" : "sort_particles",
                                  "prtcl_deposit_current" : "deposit_current",
                                  "prtcl_reflect_particles" : "reflect_particles",
                                  "prtcl_advance_reflector_walls" : "advance_reflector_walls" }

                if method in method_mapper:
                    method = method_mapper[method]

                for tile in self.local_tiles():
                    getattr(tile, method)()

            elif method.startswith("grid_"):

                for tile in self.local_tiles():
                    getattr(tile, method[5:])(*vargs)

            elif method.startswith("io_"):
                match method[3:]:
                    case "emf_snapshot":
                        self._ensure_constructed_emf_writer()
                        self._emf_writer.write(self._tile_grid._corgi_grid, self.lap)
                    case "prtcl_snapshot":
                        self._ensure_constructed_prtcl_writers()
                        for writer in self._prtcl_writers.values():
                            writer.write(self._tile_grid._corgi_grid, self.lap)
                    case "average_kinetic_energy":
                        _write_average_kinetic_energy(self.lap, self._io_config["kinetic_energy_path"], self._tile_grid._corgi_grid)
                    case "average_B_energy_density":
                        _write_average_B_energy_density(self.lap,
                                                        self._io_config["average_B_energy_density_path"],
                                                        self._tile_grid._corgi_grid)
                    case "average_E_energy_density":
                        _write_average_E_energy_density(self.lap,
                                                        self._io_config["average_E_energy_density_path"],
                                                        self._tile_grid._corgi_grid)
                    case "spectra_snapshot":
                        self._ensure_constructed_spectra_writer()
                        self._spectra_writer.write(self._tile_grid._corgi_grid, self.lap)
                    case "ram_usage":
                        self._write_ram_usage()
                    case _:
                        raise AttributeError(f"{method} is not supported IO type.")

            elif method.startswith("comm_"):
                comm_modes = [*vargs]

                for mode in comm_modes:
                    if type(mode) != comm_mode:
                        msg = "Communications only accept runko.comm_mode arguments.\n"
                        msg += f"Received {mode} of type {type(mode)}."
                        raise TypeError(msg)

                    match method[5:]:
                        case "local":
                            for mode in comm_modes:
                                self._tile_grid._corgi_grid.local_communication(mode.value)

                        case "external":
                            for mode in comm_modes:
                                handshake_mode = _virtual_tile_sync_handshake_mode(mode)

                                if handshake_mode:
                                    self._logger.debug("Starting a handshake.")
                                    self._tile_grid._corgi_grid.recv_data(handshake_mode)
                                    self._tile_grid._corgi_grid.send_data(handshake_mode)
                                    self._tile_grid._corgi_grid.wait_data(handshake_mode)

                                self._logger.debug("Starting virtual tile sync.")
                                self._tile_grid._corgi_grid.recv_data(mode.value)
                                self._tile_grid._corgi_grid.send_data(mode.value)
                                self._tile_grid._corgi_grid.wait_data(mode.value)
                        case _:
                            raise AttributeError(f"{method} is not supported communication type.")
            else:
                raise RuntimeError(f"{method} is not supported!")


        pre_post = dict(pre=pre, post=post) if not disable_timing else dict()
        lap_function(MethodWrapper(action, **pre_post))

        lap_wall_time_end = time.time()
        if not disable_timing:
            self._lap_timers.append(lap_timer)
            self._lap_wall_times.append(lap_wall_time_end - lap_wall_time_begin)


    def prelude(self, lap_function):
        """
        Execute a given lap function without increasing the lap.
        """

        if self.verbose_laps:
            self._logger.info("Executing prelude lap function...")
        self._execute_lap_function(lap_function, disable_timing=True)


    def for_one_lap(self, lap_function):
        """
        Advance simulation by one lap using given lap functions.
        """

        if self.verbose_laps:
            self._logger.info(f"Executing lap: {self.lap}")

        self._execute_lap_function(lap_function)
        self._lap += 1


    def for_each_lap(self, lap_function):
        """
        Advance simulation until `Nt` laps is reached using given lap function.
        """

        while self._lap < self._last_lap:
            self.for_one_lap(lap_function)


    def get_time_statistics(self):
        return timer_statistics(self._lap_timers)


    def reset_timers(self):
        self._lap_timers = []
        self._lap_wall_times = []


    def log_timer_statistics(self, level=logging.INFO):
        stats_dict = self.get_time_statistics()
        if len(stats_dict) == 0:
            self._logger.warning("Trying to log timer non-existing statistic.")
            return
        stats = list(stats_dict.items())
        stats.sort(key=lambda x: -x[1].total)

        total_elapsed_time = 0
        for _, s in stats:
            total_elapsed_time += s.total

        nlen = len(max(stats, key=lambda x: len(x[0]))[0])

        msg = "Simulation execution time statistics:\n"
        msg += f"{'name':<{nlen}} | total [s] | % of total | average [s] | std [s] | count\n"
        for name, s in stats:
            p = 100 * s.total / total_elapsed_time
            msg += f"{name:<{nlen}} | {s.total:>9.1e} | {p:>10.4} | {s.average:>11.3e} | {s.std_dev:>7.1e} | {s.count:>5}\n"
        msg += f"Total elapsed time: {total_elapsed_time:.4}s\n"


        total_wall_time = np.sum(self._lap_wall_times)
        avg_lap_wall_time = np.mean(self._lap_wall_times)
        laps = len(self._lap_wall_times)

        msg += f"Lap wall times: {total_wall_time:.4} s / {laps} laps = {avg_lap_wall_time:.4} s / lap"
        self._logger.log(level, msg)
