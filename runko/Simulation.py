from pytools import MethodWrapper
from pyrunko.tools import comm_mode
from pyrunko._runko_next import _virtual_tile_sync_handshake_mode
import pyrunko.emf2.threeD as emf
from .runko_logging import runko_logger


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



    def _execute_lap_function(self, lap_function):
        """
        Execute a given lap function.

        FIXME: Define lap function.
        FIXME: Parallelize the loop.
        """

        def for_each_local_tile_action(name: str):
            for tile in self.local_tiles():
                getattr(tile, name)()

        def communication_action(name: str, *comm_modes: comm_mode):
            for mode in comm_modes:
                if type(mode) != comm_mode:
                    msg = "Communications only accept runko.comm_mode arguments.\n"
                    msg += f"Received {mode} of type {type(mode)}."
                    raise TypeError(msg)

            match name:
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
                    raise AttributeError(f"{name} is not supported communication type.")

        def io_action(name: str):
            match name:
                case "emf_snapshot":
                    self._ensure_constructed_emf_writer()
                    self._emf_writer.write(self._tile_grid._corgi_grid, self.lap)
                case _:
                    raise AttributeError(f"{name} is not supported IO type.")

        for_each_local_tile = MethodWrapper(for_each_local_tile_action)
        communications = MethodWrapper(communication_action)
        io = MethodWrapper(io_action)

        lap_function(for_each_local_tile, communications, io)

        self._logger.debug(f"Executed lap function at lap: {self.lap}")


    def prelude(self, lap_function):
        """
        Execute a given lap function without increasing the lap.
        """

        self._logger.info("Executing prelude lap function...")
        self._execute_lap_function(lap_function)


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
