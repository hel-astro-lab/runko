import unittest
import runko
import numpy as np

def create_test_grid():
    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.Nt = 4
    config.Nx = 1
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 7
    config.NyMesh = 7
    config.NzMesh = 7
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"

    return config, runko.TileGrid(config)


class simulation(unittest.TestCase):
    def test_prelude_does_not_increase_lap_count(self):
        conf, tile_grid = create_test_grid()

        one_field = lambda x, y, z: (1, 1, 1)
        zero_field = lambda x, y, z: (0, 0, 0)

        for idx in tile_grid.local_tile_indices():
            tile = runko.emf.Tile(idx, conf)
            tile.set_EBJ(zero_field, one_field, zero_field)
            tile_grid.add_tile(tile, idx)

        simulation = tile_grid.configure_simulation(conf)

        for tile in simulation.local_tiles():
            _, (hBx, hBy, hBz), _ = tile.get_EBJ_with_halo()

            self.assertFalse(np.all(hBx == 1))

        prelude_ran = [False] # Has to be mutable so func below can modify it.
        def prelude_func(_, comm, *__):
            comm.virtual_tile_sync(runko.comm_mode.emf_B)
            comm.pairwise_moore(runko.comm_mode.emf_B)
            prelude_ran[0] = True

        self.assertEqual(simulation.lap, 0)
        self.assertFalse(prelude_ran[0])

        simulation.prelude(prelude_func)

        self.assertEqual(simulation.lap, 0)
        self.assertTrue(prelude_ran[0])

        for tile in simulation.local_tiles():
            _, (hBx, hBy, hBz), _ = tile.get_EBJ_with_halo()

            self.assertTrue(np.all(hBx == 1))

        simulation.for_one_lap(lambda *_: None)
        self.assertEqual(simulation.lap, 1)


if __name__ == "__main__":
    unittest.main()
