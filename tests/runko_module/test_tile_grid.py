import unittest
import runko

class tile_grid(unittest.TestCase):
    def setUp(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "hilbert_curve"
        conf.Nx = 2
        conf.Ny = 2
        conf.Nz = 2
        conf.xmin = 0
        conf.xmax = 1
        conf.ymin = 0
        conf.ymax = 2
        conf.zmin = 0
        conf.zmax = 3

        self.test_grid = runko.TileGrid(conf)

    def test_mandatory_configuration_variables(self):
        empty_conf = runko.Configuration(None)

        with self.assertRaises(RuntimeError):
            runko.TileGrid(empty_conf)

    def test_mins_cant_be_larger_than_maxs(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "hilbert_curve"
        conf.Nx = 2
        conf.Ny = 2
        conf.Nz = 2
        conf.xmin = 0
        conf.xmax = 1
        conf.ymin = 0
        conf.ymax = 0
        conf.zmin = 0
        conf.zmax = 3

        with self.assertRaisesRegex(RuntimeError, r"ymin >= ymax"):
            runko.TileGrid(conf)

        conf.ymax = 2
        conf.zmax = -1

        with self.assertRaisesRegex(RuntimeError, r"zmin >= zmax"):
            runko.TileGrid(conf)

    def test_nonsense_tile_partitioning(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "nonsense"
        conf.Nx = 2
        conf.Ny = 2
        conf.Nz = 2
        conf.xmin = 0
        conf.xmax = 1
        conf.ymin = 0
        conf.ymax = 1
        conf.zmin = 0
        conf.zmax = 1

        with self.assertRaisesRegex(RuntimeError, r"tile_partitioning"):
            runko.TileGrid(conf)

    def test_catepillar_track_requires_length(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "catepillar_track"
        conf.Nx = 2
        conf.Ny = 2
        conf.Nz = 2
        conf.xmin = 0
        conf.xmax = 1
        conf.ymin = 0
        conf.ymax = 1
        conf.zmin = 0
        conf.zmax = 1

        with self.assertRaisesRegex(RuntimeError, r"catepillar_track_length"):
            runko.TileGrid(conf)

if __name__ == "__main__":
    unittest.main()
