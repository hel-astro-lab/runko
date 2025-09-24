import unittest
import itertools
import numpy as np

import runko

class emf_FDTD2(unittest.TestCase):

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

    def test_push_half_b_x(self):
        """
        E = (0, z, -y) => curl(E) = (-2, 0, 0) => push_half_b => B = (b, 0, 0)
        """

        tile = self.make_test_tile()

        E = lambda x, y, z: (0, z, -y)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = Bx[1:-1, 1:-1, 1:-1].flat[0]
        self.assertTrue(b > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], b, places=5)
            self.assertAlmostEqual(By[i, j, k], 0, places=5)
            self.assertAlmostEqual(Bz[i, j, k], 0, places=5)

    def test_push_half_b_y(self):
        """
        E = (z, 0, -x) => curl(E) = (0, 2, 0) => push_half_b => B = (0, -b, 0)
        """

        tile = self.make_test_tile()

        E = lambda x, y, z: (z, 0, -x)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = By[1:-1, 1:-1, 1:-1].flat[0]
        self.assertTrue(b < 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], 0, places=5)
            self.assertAlmostEqual(By[i, j, k], b, places=5)
            self.assertAlmostEqual(Bz[i, j, k], 0, places=5)

    def test_push_half_b_z(self):
        """
        E = (y, -x, 0) => curl(E) = (0, 0, -2) => push_half_b => B = (0, 0, b)
        """

        tile = self.make_test_tile()

        E = lambda x, y, z: (y, -x, 0)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = Bz[1:-1, 1:-1, 1:-1].flat[0]
        self.assertTrue(b > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], 0, places=5)
            self.assertAlmostEqual(By[i, j, k], 0, places=5)
            self.assertAlmostEqual(Bz[i, j, k], b, places=5)

    def test_push_e_x(self):
        """
        B = (0, z, -y) => curl(B) = (-2, 0, 0) => push_e => E = (-e, 0, 0)
        """

        tile = self.make_test_tile()

        E = lambda x, y, z: (0, 0, 0)
        B = lambda x, y, z: (0, z, -y)
        J = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, B, J)
        tile.push_e()

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        e = Ex[1:-1, 1:-1, 1:-1].flat[0]
        self.assertTrue(e < 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Ex[i, j, k], e, places=5)
            self.assertAlmostEqual(Ey[i, j, k], 0, places=5)
            self.assertAlmostEqual(Ez[i, j, k], 0, places=5)

    def test_push_e_y(self):
        """
        B = (z, 0, -x) => curl(B) = (0, 2, 0) => push_e => E = (0, e, 0)
        """

        tile = self.make_test_tile()

        E = lambda x, y, z: (0, 0, 0)
        B = lambda x, y, z: (z, 0, -x)
        J = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, B, J)
        tile.push_e()

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        e = Ey[1:-1, 1:-1, 1:-1].flat[0]
        self.assertTrue(e > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Ex[i, j, k], 0, places=5)
            self.assertAlmostEqual(Ey[i, j, k], e, places=5)
            self.assertAlmostEqual(Ez[i, j, k], 0, places=5)

    def test_push_e_z(self):
        """
        B = (y, -x, 0) => curl(B) = (0, 0, -2) => push_e => E = (0, 0, -e)
        """

        tile = self.make_test_tile()

        E = lambda x, y, z: (0, 0, 0)
        B = lambda x, y, z: (y, -x, 0)
        J = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, B, J)
        tile.push_e()

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        e = Ez[1:-1, 1:-1, 1:-1].flat[0]
        self.assertTrue(e < 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Ex[i, j, k], 0, places=5)
            self.assertAlmostEqual(Ey[i, j, k], 0, places=5)
            self.assertAlmostEqual(Ez[i, j, k], e, places=5)

if __name__ == "__main__":
    unittest.main()
