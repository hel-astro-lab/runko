import unittest
import itertools
import numpy as np

import runko

class emf2_current_filter_binomial2(unittest.TestCase):

    def setUp(self):
        self.config = runko.Configuration(None)
        self.config.Nx = 4
        self.config.Ny = 4
        self.config.Nz = 4
        self.config.NxMesh = 10
        self.config.NyMesh = 11
        self.config.NzMesh = 13
        self.config.xmin = 0
        self.config.ymin = 0
        self.config.zmin = 0
        self.config.cfl = 1
        self.config.field_propagator = "FDTD2"
        self.config.current_filter = "binomial2"

        def f():
            tile_grid_idx = (0, 0, 0)
            return runko.emf.Tile(tile_grid_idx, self.config)

        self.make_test_tile = f

        # Halo regions can not easily be updated in this unittest,
        # se we only test inner region of the mesh.
        def index_space():
            return itertools.product(range(1, -1 + self.config.NxMesh),
                                     range(1, -1 + self.config.NyMesh),
                                     range(1, -1 + self.config.NzMesh))

        self.get_interior_index_space = index_space

        n = self.config.NxMesh - 2
        m = self.config.NyMesh - 2
        l = self.config.NzMesh - 2

        c = 0
        for _ in self.get_interior_index_space():
            c += 1

        self.assertEqual(c, n * m * l)


    def test_constant_J_is_unchanged(self):

        tile = self.make_test_tile()

        E = lambda x, y, z: (0, 0, 0)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (1, 2, 3)

        tile.set_EBJ(E, B, J)

        def assertConstantJ():
            _, _, (Jx, Jy, Jz) = tile.get_EBJ()

            for i, j, k in self.get_interior_index_space():
                self.assertAlmostEqual(Jx[i, j, k], 1, places=5)
                self.assertAlmostEqual(Jy[i, j, k], 2, places=5)
                self.assertAlmostEqual(Jz[i, j, k], 3, places=5)

        assertConstantJ()
        tile.filter_current()
        assertConstantJ()


    def test_alternating_J_smooths_out(self):

        tile = self.make_test_tile()

        E = lambda x, y, z: (0, 0, 0)
        B = lambda x, y, z: (0, 0, 0)
        def J(x, y, z):
            c = int(x + y + z) % 2
            return (c, c, c)

        tile.set_EBJ(E, B, J)

        _, _, (J0x, J0y, J0z) = tile.get_EBJ()

        def almost_equal(a, b, epsilon=0.00001):
            return epsilon - abs(a - b) > 0

        for i, j, k in self.get_interior_index_space():
            xy_same = almost_equal(J0x[i, j, k], J0y[i, j, k])
            yz_same = almost_equal(J0y[i, j, k], J0z[i, j, k])
            all_same = xy_same and yz_same
            self.assertTrue(all_same)

            zero_or_one = J0x[i, j, k] == 0 or J0x[i, j, k] == 1
            self.assertTrue(zero_or_one)

        tile.filter_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        for i, j, k in self.get_interior_index_space():
            xy_same = almost_equal(Jx[i, j, k], Jy[i, j, k])
            yz_same = almost_equal(Jy[i, j, k], Jz[i, j, k])
            all_same = xy_same and yz_same
            self.assertTrue(all_same)

            between_zero_and_one = Jx[i, j, k] > 0 and Jx[i, j, k] < 1
            self.assertTrue(between_zero_and_one)


if __name__ == "__main__":
    unittest.main()
