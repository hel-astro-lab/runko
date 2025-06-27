import pytools
import pycorgi.threeD as pycorgi

class TileGrid:
    """
    Represents the grid of runko tiles.
    Stores tiles in current locality and global information about all tiles.

    Initially tiles are distributed to different localities based on `tile_partitioning`.
    Possible values are: "hilbert_curve" and "catepillar_track"

    "catepillar_track" requires to be configured with `catepillar_track_length` (int).
    """

    def __init__(self, conf):

        required_vars = ["Nx",
                         "Ny",
                         "Nz",
                         "NxMesh",
                         "NyMesh",
                         "NzMesh",
                         "xmin",
                         "ymin",
                         "zmin",
                         "tile_partitioning"]

        for var in required_vars:
            if getattr(conf, var) is None:
                raise RuntimeError(f"Can not construct TileGrid without: {var}")

        valid_tile_partitions = ["hilbert_curve", "catepillar_track"]

        if conf.tile_partitioning not in valid_tile_partitions:
            raise RuntimeError(f"invalid `tile_partitioning`: {conf.tile_partitioning}")

        if conf.tile_partitioning == "catepillar_track":
            if conf.catepillar_track_length is None:
                msg = "Catepillar track requested, but no `catepillar_track_length` given!"
                raise RuntimeError(msg)


        self._corgi_grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)

        xmax = conf.xmin + conf.Nx * conf.NxMesh
        ymax = conf.ymin + conf.Ny * conf.NyMesh
        zmax = conf.zmin + conf.Nz * conf.NzMesh

        self._corgi_grid.set_grid_lims(conf.xmin,
                                       xmax,
                                       conf.ymin,
                                       ymax,
                                       conf.zmin,
                                       zmax)

        legacy_conf = type("", (), dict(oneD=False,
                                        twoD=False,
                                        threeD=False,
                                        mpi_task_mode=False))

        if conf.tile_partitioning == "hilbert_curve":
            pytools.balance_mpi(self._corgi_grid, legacy_conf, do_print=False)
        elif conf.tile_partitioning == "catepillar_track":
            pytools.load_catepillar_track_mpi(self._corgi_grid,
                                              conf.catepillar_track_length,
                                              legacy_conf,
                                              do_print=False)
        else:
            raise RuntimeError("Due to previous checking this should not happend.")

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
