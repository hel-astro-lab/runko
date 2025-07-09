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


    @property
    def lap(self):
        """Current simulation lap."""
        return self._lap


    def for_one_lap(self, lap_function):
        """Advance simulation by one lap using given lap functions."""

        self._lap += 1


    def for_each_lap(self, lap_function):
        """
        Advance simulation until `Nt` laps is reached using given lap function.
        """

        while self._lap < self._last_lap:
            self.for_one_lap(lap_function)
