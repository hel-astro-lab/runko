# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

from .method_wrapper import MethodWrapper
from .runko_logging import runko_logger, on_main_rank
from .runko_timer import Timer, timer_statistics
from .configuration import Configuration
from .tiles import ProxyTile
import runko_cpp_bindings.actions as actions

import json
import logging
import pickle
import time
import numpy as np
import pathlib
from dataclasses import dataclass
from mpi4py import MPI
from .ram_usage import get_rss_kB, get_gpu_mem_kB

@dataclass
class wait_measurement:
    label: str
    begin: float
    end: float


class Simulation:
    """
    Handle to a configured runko simulation.
    """

    _im_not_user = 42

    def __init__(self, initial_tiles, prevent_user_init: int, **kwargs):
        """
        Construct runko simulation from tile grid.

        `prevent_user_init` is to make sure users don't create
        objects of this calss on their own.
        """

        if prevent_user_init != Simulation._im_not_user:
            raise RuntimeError("Don't instantiate this class directly!")

        self._config = kwargs["config"]
        self._runtime_instance = actions.RuntimeInstance()
        self._simulation_context = actions.SimulationContext(self._config)

        self._simulation_context.add_tiles(initial_tiles)

        set_neighs = (actions.set_cartesian_neighbors, actions.current_context)
        set_comm_infos = (actions.set_cartesian_comm_infos, actions.current_context)

        self._simulation_context.eval((actions.sequence,
                                       (actions.quote, set_neighs),
                                       (actions.quote, set_comm_infos)))

        self._lap = 0

        self._last_lap = kwargs['Nt']
        self._prev_elapsed_wall_time = 0
        self._total_interval_time = 0

        self._io_config = kwargs['io_config']

        self._lap_timers = []
        self._lap_wall_times = []
        self._lap_finish_times = []

        # Each lap has their own list.
        self._wait_measurements: list[list[wait_measurement]] = []

        self.verbose_laps = kwargs.get('verbose', True)

        self._logger = runko_logger("Simulation")

        self._rank = MPI.COMM_WORLD.Get_rank()

        self._ram_file = pathlib.Path(f"{self._io_config['outdir']}/mem-usage/cpu/{self._rank}.csv")
        self._gpu_ram_file = pathlib.Path(f"{self._io_config['outdir']}/mem-usage/gpu/{self._rank}.csv")

        self._timer_stats_file = pathlib.Path(f"{self._io_config['outdir']}/timer-statistics/{self._rank}.pkl")

        self._io_config['kinetic_energy_path'] = self._io_config["outdir"] + "/average_kinetic_energy.txt"
        self._io_config['average_B_energy_density_path'] = self._io_config["outdir"] + "/average_B_energy_density.txt"
        self._io_config['average_E_energy_density_path'] = self._io_config["outdir"] + "/average_E_energy_density.txt"

        # Reset txt output files.
        for name in ('kinetic_energy_path', 'average_B_energy_density_path', 'average_E_energy_density_path'):
            pathlib.Path(self._io_config[name]).unlink(missing_ok=True)


        self._trace_file = pathlib.Path(f"{self._io_config['outdir']}/traces/{self._rank}.json")

        ctor_msg = "Simulation constructed with:\n"
        ctor_msg += f"\tNt = {kwargs['Nt']}\n"
        ctor_msg += f"\tio config: {self._io_config}"
        self._logger.debug(ctor_msg)


    def virtual_tiles(self):
        """
        Return iterable which goes through all virtual tiles in current rank.
        """

        raise NotImplementedError()


    def local_tiles(self):
        """
        Return iterable which goes through all local tiles in current rank.
        """

        ids = self._simulation_context.get_local_tile_ids()
        return [ProxyTile(x, self._simulation_context) for x in ids]


    def _boundary_tiles(self):
        """
        Return iterable which goes through all boundary tiles in current rank.
        """

        raise NotImplementedError()


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

        if not self._ram_file.exists():
            self._ram_file.parent.mkdir(exist_ok=True, parents=True)
            self._ram_file.write_text("lap,ram usage [kB]\n")

        with open(self._ram_file, "a") as f:
            f.write(f"{self.lap},{get_rss_kB()}\n")

        gpu_mem = get_gpu_mem_kB()

        if gpu_mem:
            if not self._gpu_ram_file.exists():
                self._gpu_ram_file.parent.mkdir(exist_ok=True, parents=True)
                self._gpu_ram_file.write_text("lap,gpu mem usage [kB]\n")

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
        self._wait_measurements.append([])

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
                                  "prtcl_pack_outgoing" : "pack_outgoing_particles",
                                  "prtcl_deposit_current" : "deposit_current",
                                  "prtcl_reflect_particles" : "reflect_particles",
                                  "prtcl_advance_reflector_walls" : "advance_reflector_walls" }

                if method in method_mapper:
                    method = method_mapper[method]

                if method == "pack_outgoing_particles":
                    if on_main_rank():
                        self._logger.warn("pack_outgoing_particles is depricated and implemented as no-op.")
                    return

                symbol = getattr(actions, method)
                self._simulation_context.eval((symbol, actions.current_context))

            elif method.startswith("grid_"):
                symbol = getattr(actions, method[len("grid_"):])
                self._simulation_context.eval((symbol, actions.current_context))

            elif method.startswith("io_"):
                match method[3:]:
                    case "emf_snapshot":
                        raise NotImplementedError()
                    case "prtcl_snapshot":
                        raise NotImplementedError()
                    case "average_kinetic_energy":
                        raise NotImplementedError()
                    case "average_B_energy_density":
                        raise NotImplementedError()
                    case "average_E_energy_density":
                        raise NotImplementedError()
                    case "spectra_snapshot":
                        raise NotImplementedError()
                    case "ram_usage":
                        raise NotImplementedError()
                    case _:
                        raise AttributeError(f"{method} is not supported IO type.")

            elif method.startswith("comm_"):
                if method == "comm_external":
                    symbol = actions.comm_external
                elif method == "comm_local":
                    symbol = actions.comm_local
                else:
                    raise RuntimeError(f"Unregonized communication: {method}")

                mode_to_prog = lambda m: (actions.quote, (symbol, actions.current_context, m))
                prog = (actions.sequence,)
                modes = (*vargs,)
                for mode in modes:
                    prog = prog + (mode_to_prog(mode),)
                self._simulation_context.eval(prog)
            else:
                raise RuntimeError(f"{method} is not supported!")


        pre_post = dict(pre=pre, post=post) if not disable_timing else dict()
        lap_function(MethodWrapper(action, **pre_post))

        lap_wall_time_end = time.time()
        if not disable_timing:
            self._lap_timers.append(lap_timer)
            self._lap_wall_times.append(lap_wall_time_end - lap_wall_time_begin)
            self._lap_finish_times.append(lap_wall_time_end)


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


    def get_time_statistics(self, n_latests=None):
        if n_latests:
            return timer_statistics(self._lap_timers[-n_latests:])
        else:
            return timer_statistics(self._lap_timers)


    def reset_timers(self):
        self._lap_timers = []
        self._lap_wall_times = []
        self._prev_elapsed_wall_time = self._total_interval_time


    def log_timer_statistics(self, level=logging.INFO):
        stats_dict = self.get_time_statistics(self._io_config["laps_in_timer_statistics"])
        if len(stats_dict) == 0:
            self._logger.warning("Trying to log timer non-existing statistic.")
            return
        stats = list(stats_dict.items())
        stats.sort(key=lambda x: -x[1].total)

        self._total_interval_time = 0
        for _, s in stats:
            self._total_interval_time += s.total 

        nlen = len(max(stats, key=lambda x: len(x[0]))[0])

        msg = "Simulation execution time statistics:\n"
        msg += f"{'name':<{nlen}} | total [s] | % of total | average [s] | std [s] | count\n"
        for name, s in stats:
            p = 100 * s.total / self._total_interval_time
            msg += f"{name:<{nlen}} | {s.total:>9.1e} | {p:>10.4} | {s.average:>11.3e} | {s.std_dev:>7.1e} | {s.count:>5}\n"
        msg += f"Total elapsed interval time: {self._total_interval_time:.4}s\n"

        total_wall_time = np.sum(self._lap_wall_times)
        avg_lap_wall_time = np.mean(self._lap_wall_times)
        laps = len(self._lap_wall_times)

        msg += f"Lap wall times: {total_wall_time:.4} s / {laps} laps = {avg_lap_wall_time:.4} s / lap\n"

        # add previously elapsed wall time before timer resets
        self._total_interval_time += self._prev_elapsed_wall_time 
        msg += f"Total elapsed wall time: {self._total_interval_time:.4}s\n"

        lap_per_max_lap = self._lap/self._last_lap # simulation progress in fraction 
        msg += f"Current lap: {self._lap} / {self._last_lap} ({lap_per_max_lap:.1%})"
        self._logger.log(level, msg)


    def write_trace_json(self, combine: bool = True) -> dict:
        """
        Writes trace data of each process into `<outdir>/traces/<mpi-rank>.json`
        in Trace Event Format [0] which can be viewed with e.g. Pefetto UI.
        Configuration parameter `laps_in_timer_statistics` defines number of latest laps
        from which the data is written. If it is not given, data from all laps is written.

        If `combine` is `True`, this has to be called on every rank
        and the trace jsons are combined into `<outdir>/traces/combined.json`.

        [0]: https://docs.google.com/document/d/1CvAClvFfyA5R-PhYUmn5OOQtYMH4h6I0nSsKchNAySU
        """

        obj = dict(traceEvents=[{
            "name": "process_name",
            "ph": "M",
            "pid": 0, # unused
            "tid": self._rank, # Use tid as proxy for mpi rank.
            "args": { "name": f"MPI ranks" }
        },{
            "name": "thread_name",
            "ph": "M",
            "pid": 0, # unused
            "tid": self._rank, # Use tid as proxy for mpi rank.
            "args": { "name": f"MPI rank {self._rank}" }
        }])

        n_latests = self._io_config["laps_in_timer_statistics"]
        if n_latests:
            timers = self._lap_timers[-n_latests:]
            lap_finish_times = self._lap_finish_times[-n_latests:]
            waits = self._wait_measurements[-n_latests:]
        else:
            timers = self._lap_timers
            lap_finish_times = self._lap_finish_times
            waits = self._wait_measurements

        for timer in timers:
            for name, measurement in timer.time_measurements.items():
                cat, _ = name.split("_", 1)

                obj["traceEvents"].append({
                    "name": name,
                    "cat": cat,
                    "ph": "X",
                    "pid": 0, # unused
                    "tid": self._rank, # Use tid as proxy for mpi rank.
                    "ts": 1e6 * measurement.begin,
                    "dur": 1e6 * (measurement.end - measurement.begin),
                })

        for t in lap_finish_times:
            obj["traceEvents"].append({
                "name": f"lap_finish",
                "cat": "meta",
                "ph": "i",
                "pid": 0, # unused
                "tid": self._rank, # Use tid as proxy for mpi rank.
                "ts": 1e6 * t,
            })

        # Iterate over flattened waits.
        for t in [x for y in waits for x in y]:
            obj["traceEvents"].append({
                "name": t.label,
                "cat": "wait",
                "ph": "X",
                "pid": 0, # unused
                "tid": self._rank, # Use tid as proxy for mpi rank.
                "ts": 1e6 * t.begin,
                "dur": 1e6 * (t.end - t.begin),
            })

        self._trace_file.parent.mkdir(exist_ok=True, parents=True)
        with open(self._trace_file, "w") as f:
            json.dump(obj, f)

        if not combine:
            return

        MPI.COMM_WORLD.barrier()
        if on_main_rank():
            combined = {"traceEvents": []}

            for i in range(MPI.COMM_WORLD.Get_size()):
                p = pathlib.Path(f"{self._io_config['outdir']}/traces/{i}.json")
                with open(p, "r") as f:
                    data = json.load(f)
                    combined["traceEvents"] += data["traceEvents"]

            with open(f"{self._io_config['outdir']}/traces/combined.json", "w") as f:
                json.dump(combined, f)


    def pickle_timer_statistics(self):
        """
        Writes pickled dictionary of timer statistics to `<outdir>/timer-statistics/<rank>.pkl`
        which maps component name (str) to runko.TimerStatistics.
        """

        self._timer_stats_file.parent.mkdir(exist_ok=True, parents=True)
        with open(self._timer_stats_file, "wb") as f:
            pickle.dump(self.get_time_statistics(self._io_config["laps_in_timer_statistics"]), f)
