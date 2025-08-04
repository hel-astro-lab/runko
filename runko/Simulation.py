from pytools import MethodWrapper
from pyrunko.tools import comm_mode
from pyrunko._runko_next import _virtual_tile_sync_handshake_mode
import pyrunko.emf2.threeD as emf
from .runko_logging import runko_logger
from .runko_timer import Timer, timer_statistics
import logging


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

        self._lap_timers = []

        self._logger = runko_logger("Simulation")

        ctor_msg = "Simulation constructed with:\n"
        ctor_msg += f"\tNt = {kwargs['Nt']}\n"
        ctor_msg += f"\tio config: {self._io_config}"
        self._logger.debug(ctor_msg)


    def _ensure_constructed_emf_writer(self):
        if self._emf_writer:
            return

        self._emf_writer = emf.FieldsWriter2(self._io_config["outdir"],
                                             self._tile_grid._Nx,
                                             self._tile_grid._NxMesh,
                                             self._tile_grid._Ny,
                                             self._tile_grid._NyMesh,
                                             self._tile_grid._Nz,
                                             self._tile_grid._NzMesh,
                                             self._io_config["stride"])

        self._logger.debug("FieldsWriter2 constructed.")


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



    def _execute_lap_function(self, lap_function, disable_timing=False):
        """
        Execute a given lap function.

        FIXME: Define lap function.
        FIXME: Parallelize the loop.
        """

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

        def for_each_local_tile_action(method: str, **kwargs):
            for tile in self.local_tiles():
                getattr(tile, method)()

        def communication_action(method: str, *comm_modes: comm_mode, **kwargs):
            for mode in comm_modes:
                if type(mode) != comm_mode:
                    msg = "Communications only accept runko.comm_mode arguments.\n"
                    msg += f"Received {mode} of type {type(mode)}."
                    raise TypeError(msg)

            match method:
                case "pairwise_moore":
                    for mode in comm_modes:
                        self._tile_grid._corgi_grid.pairwise_moore_communication(mode.value)

                case "virtual_tile_sync":
                    for mode in comm_modes:
                        handshake_mode = _virtual_tile_sync_handshake_mode(mode)

                        if handshake_mode:
                            self._tile_grid._corgi_grid.recv_data(handshake_mode)
                            self._tile_grid._corgi_grid.send_data(handshake_mode)
                            self._tile_grid._corgi_grid.wait_data(handshake_mode)

                        self._tile_grid._corgi_grid.recv_data(mode.value)
                        self._tile_grid._corgi_grid.send_data(mode.value)
                        self._tile_grid._corgi_grid.wait_data(mode.value)
                case _:
                    raise AttributeError(f"{method} is not supported communication type.")

        def io_action(method: str, **kwargs):
            match method:
                case "emf_snapshot":
                    self._ensure_constructed_emf_writer()
                    self._emf_writer.write(self._tile_grid._corgi_grid, self.lap)
                case _:
                    raise AttributeError(f"{method} is not supported IO type.")

        pre_post = dict(pre=pre, post=post) if not disable_timing else dict()
        for_each_local_tile = MethodWrapper(for_each_local_tile_action, **pre_post)
        communications = MethodWrapper(communication_action, **pre_post)
        io = MethodWrapper(io_action, **pre_post)

        lap_function(for_each_local_tile, communications, io)

        if not disable_timing:
            self._lap_timers.append(lap_timer)


    def prelude(self, lap_function):
        """
        Execute a given lap function without increasing the lap.
        """

        self._logger.info("Executing prelude lap function...")
        self._execute_lap_function(lap_function, disable_timing=True)


    def for_one_lap(self, lap_function):
        """
        Advance simulation by one lap using given lap functions.
        """

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


    def log_timer_statistics(self, level=logging.INFO):
        stats = list(self.get_time_statistics().items())
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
        msg += f"Total elapsed time [s]: {total_elapsed_time}"
        self._logger.log(level, msg)
