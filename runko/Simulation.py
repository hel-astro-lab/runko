class Simulation:
    """
    Handle to a configured runko simulation.
    """

    _im_not_user = 42

    def __init__(self, tile_grid, prevent_user_init: int):
        """
        Construct runko simulation from tile grid.

        `prevent_user_init` is to make sure users don't create
        objects of this calss on their own.
        """

        if prevent_user_init != Simulation._im_not_user:
            raise RuntimeError("Don't instantiate this class directly!")

        self._tile_grid = tile_grid


    def virtual_tiles(self):
        """
        Return iterable which goes through all virtual tiles in current rank.
        """

        for vtile_id in self._tile_grid._corgi_grid.get_virtual_tiles():
            yield self._tile_grid._corgi_grid.get_tile(vtile_id)
