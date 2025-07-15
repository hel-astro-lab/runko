from pytools import MethodWrapper
from pyrunko.tools import comm_mode

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
                        self._tile_grid._corgi_grid.recv_data(mode.value)
                        self._tile_grid._corgi_grid.send_data(mode.value)
                        self._tile_grid._corgi_grid.wait_data(mode.value)
                case _:
                    raise AttributeError(f"{name} is not supported communication type.")

        for_each_local_tile = MethodWrapper(for_each_local_tile_action)
        communications = MethodWrapper(communication_action)

        lap_function(for_each_local_tile, communications)


    def prelude(self, lap_function):
        """
        Execute a given lap function without increasing the lap.
        """

        self._execute_lap_function(lap_function)


    def for_one_lap(self, lap_function):
        """
        Advance simulation by one lap using given lap functions.
        """

        self._execute_lap_function(lap_function)
        self._lap += 1


    def for_each_lap(self, lap_function):
        """
        Advance simulation until `Nt` laps is reached using given lap function.
        """

        while self._lap < self._last_lap:
            self.for_one_lap(lap_function)
