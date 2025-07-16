import unittest
import runko

def make_valid_emf2_config():
    config = runko.Configuration(None)
    config.Nx = 4
    config.Ny = 4
    config.Nz = 4
    config.NxMesh = 10
    config.NyMesh = 11
    config.NzMesh = 13
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"

    return config


class pic2_tile(unittest.TestCase):

    def test_pic_tiles_are_initially_empty(self):

        config = make_valid_emf2_config()
        config.qe = 1
        config.me = 1
        config.delgam = 1.0e-5
        config.temperature_ratio = 1.0
        config.sigma = 40
        config.c_omp = 1
        config.ppc = 1

        tile_grid_idx = (0, 1, 2)
        tile = runko.pic.Tile(tile_grid_idx, config)

        (pos0_x, pos0_y, pos0_z) = tile.get_positions(runko.particle.electron)
        (vel0_x, vel0_y, vel0_z) = tile.get_velocities(runko.particle.electron)

        self.assertEqual(0, len(pos0_x))
        self.assertEqual(0, len(pos0_y))
        self.assertEqual(0, len(pos0_z))
        self.assertEqual(0, len(vel0_x))
        self.assertEqual(0, len(vel0_y))
        self.assertEqual(0, len(vel0_z))

        weights0 = tile.get_weights(runko.particle.electron)
        self.assertEqual(0, len(weights0))


if __name__ == "__main__":
    unittest.main()
