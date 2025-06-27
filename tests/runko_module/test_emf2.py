import unittest
import re
import itertools
import numpy as np

import runko

class emf2(unittest.TestCase):

    def test_unregonized_config_value_type_fails_gracefully(self):

        config = runko.Configuration(None)
        config.Nx = 2
        config.Ny = 3
        config.Nz = 4
        config.NxMesh = 10
        config.NyMesh = 12
        config.NzMesh = 14
        config.xmin = -3.2
        config.ymin = -2.3
        config.zmin = 1

        class foo:
            def __inti__(self):
                self.bar = 42

        config.cfl = foo()

        tile_grid_idx = (0, 1, 2)

        regex = re.compile("cfl.*unsupported type.*foo", re.IGNORECASE)
        with self.assertRaisesRegex(RuntimeError, regex):
            runko.emf.Tile(tile_grid_idx, config)

    def test_field_set_and_get_roundtrip(self):

        config = runko.Configuration(None)
        config.Nx = 2
        config.Ny = 3
        config.Nz = 4
        config.NxMesh = 10
        config.NyMesh = 12
        config.NzMesh = 14
        config.xmin = -3.2
        config.ymin = -2.3
        config.zmin = 1
        config.cfl = 1

        tile_grid_idx = (0, 1, 2)
        tile = runko.emf.Tile(tile_grid_idx, config)

        (E0x, E0y, E0z), (B0x, B0y, B0z), (J0x, J0y, J0z) = tile.get_EBJ()

        self.assertTrue(np.all(E0x == 0))
        self.assertTrue(np.all(E0y == 0))
        self.assertTrue(np.all(E0z == 0))
        self.assertTrue(np.all(B0x == 0))
        self.assertTrue(np.all(B0y == 0))
        self.assertTrue(np.all(B0z == 0))
        self.assertTrue(np.all(J0x == 0))
        self.assertTrue(np.all(J0y == 0))
        self.assertTrue(np.all(J0z == 0))

        # Remember that E and B are stored in a Yee Lattice.
        # Grid indices (i, j, k) corresponds to location:
        #
        #     X(i, j, k) = i * dx * e_x + j * dy * e_y + k * dz * e_z
        #
        # where e_x, e_y and e_z are unit vectors.
        #
        # E and B components corrseponding to grid indices (i, j, k)
        # are located at different positions:
        #
        # Component E_l is located at: X(i, j, k) + 0.5 * dl * e_l
        # Component B_l is located at: X(i + 0.5, j + 0.5, k + 0.5) - 0.5 * dl * e_l
        # Components of current J are located at same places as E.

        Einit = lambda x, y, z: (x, y, z)
        Binit = lambda x, y, z: (2 * x, 2 * y, 2 * z)
        Jinit = lambda x, y, z: (3 * x, 3 * y, 3 * z)

        tile.set_EBJ(Einit, Binit, Jinit)

        (Ex, Ey, Ez), (Bx, By, Bz), (Jx, Jy, Jz) = tile.get_EBJ()

        index_space = itertools.product(range(config.NxMesh),
                                        range(config.NyMesh),
                                        range(config.NzMesh))
        for i, j, k in index_space:

            ii = config.xmin + tile_grid_idx[0] * config.NxMesh + i
            jj = config.ymin + tile_grid_idx[1] * config.NyMesh + j
            kk = config.zmin + tile_grid_idx[2] * config.NzMesh + k

            self.assertAlmostEqual(Ex[i, j, k], Einit(ii + 0.5, jj, kk)[0], places=5)
            self.assertAlmostEqual(Ey[i, j, k], Einit(ii, jj + 0.5, kk)[1], places=5)
            self.assertAlmostEqual(Ez[i, j, k], Einit(ii, jj, kk + 0.5)[2], places=5)

            self.assertAlmostEqual(Bx[i, j, k], Binit(ii, jj + 0.5, kk + 0.5)[0], places=5)
            self.assertAlmostEqual(By[i, j, k], Binit(ii + 0.5, jj, kk + 0.5)[1], places=5)
            self.assertAlmostEqual(Bz[i, j, k], Binit(ii + 0.5, jj + 0.5, kk)[2], places=5)

            self.assertAlmostEqual(Jx[i, j, k], Jinit(ii + 0.5, jj, kk)[0], places=5)
            self.assertAlmostEqual(Jy[i, j, k], Jinit(ii, jj + 0.5, kk)[1], places=5)
            self.assertAlmostEqual(Jz[i, j, k], Jinit(ii, jj, kk + 0.5)[2], places=5)

if __name__ == "__main__":
    unittest.main()
