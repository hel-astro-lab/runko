import unittest
import itertools
import numpy as np

import runko

class emf_stencil(unittest.TestCase):

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
        self.config.field_propagator = "stencil"
        # All zeros -> alpha=1 -> reduces to FDTD2

        def f():
            tile_grid_idx = (0, 0, 0)
            return runko.emf.threeD.Tile(tile_grid_idx, self.config)

        self.make_test_tile = f

        # Stencil reaches further into halo than FDTD2,
        # so test interior must be 3 cells in from each edge.
        def index_space():
            return itertools.product(range(3, -3 + self.config.NxMesh),
                                     range(3, -3 + self.config.NyMesh),
                                     range(3, -3 + self.config.NzMesh))

        self.get_interior_index_space = index_space

    def test_zero_coeffs_push_half_b_x(self):
        """
        With all stencil coefficients = 0 (alpha=1), push_half_b should
        match FDTD2 exactly.

        E = (0, z, -y) => curl(E) = (-2, 0, 0) => push_half_b => Bx > 0
        """

        # Stencil tile with zero coefficients
        tile_stencil = self.make_test_tile()
        E = lambda x, y, z: (0, z, -y)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)
        tile_stencil.set_EBJ(E, B, J)
        tile_stencil.push_half_b()

        # FDTD2 tile for comparison
        config2 = runko.Configuration(None)
        for attr in ['Nx', 'Ny', 'Nz', 'NxMesh', 'NyMesh', 'NzMesh',
                      'xmin', 'ymin', 'zmin', 'cfl']:
            setattr(config2, attr, getattr(self.config, attr))
        config2.field_propagator = "FDTD2"
        tile_fdtd2 = runko.emf.threeD.Tile((0, 0, 0), config2)
        tile_fdtd2.set_EBJ(E, B, J)
        tile_fdtd2.push_half_b()

        _, (Bx_s, By_s, Bz_s), _ = tile_stencil.get_EBJ()
        _, (Bx_f, By_f, Bz_f), _ = tile_fdtd2.get_EBJ()

        np.testing.assert_allclose(Bx_s, Bx_f, atol=1e-6)
        np.testing.assert_allclose(By_s, By_f, atol=1e-6)
        np.testing.assert_allclose(Bz_s, Bz_f, atol=1e-6)

    def test_push_half_b_x(self):
        """
        E = (0, z, -y) => curl(E) = (-2, 0, 0) => push_half_b => Bx > 0
        """
        tile = self.make_test_tile()
        E = lambda x, y, z: (0, z, -y)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = Bx[3:-3, 3:-3, 3:-3].flat[0]
        self.assertTrue(b > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], b, places=5)
            self.assertAlmostEqual(By[i, j, k], 0, places=5)
            self.assertAlmostEqual(Bz[i, j, k], 0, places=5)

    def test_push_half_b_y(self):
        """
        E = (z, 0, -x) => curl(E) = (0, 2, 0) => push_half_b => By < 0
        """
        tile = self.make_test_tile()
        E = lambda x, y, z: (z, 0, -x)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = By[3:-3, 3:-3, 3:-3].flat[0]
        self.assertTrue(b < 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], 0, places=5)
            self.assertAlmostEqual(By[i, j, k], b, places=5)
            self.assertAlmostEqual(Bz[i, j, k], 0, places=5)

    def test_push_half_b_z(self):
        """
        E = (y, -x, 0) => curl(E) = (0, 0, -2) => push_half_b => Bz > 0
        """
        tile = self.make_test_tile()
        E = lambda x, y, z: (y, -x, 0)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = Bz[3:-3, 3:-3, 3:-3].flat[0]
        self.assertTrue(b > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], 0, places=5)
            self.assertAlmostEqual(By[i, j, k], 0, places=5)
            self.assertAlmostEqual(Bz[i, j, k], b, places=5)

    def test_push_e_uses_fdtd2(self):
        """
        With stencil propagator, push_e should use standard FDTD2 and give
        the same result as a pure FDTD2 tile.
        """
        tile_stencil = self.make_test_tile()
        E = lambda x, y, z: (0, 0, 0)
        B = lambda x, y, z: (0, z, -y)
        J = lambda x, y, z: (0, 0, 0)
        tile_stencil.set_EBJ(E, B, J)
        tile_stencil.push_e()

        config2 = runko.Configuration(None)
        for attr in ['Nx', 'Ny', 'Nz', 'NxMesh', 'NyMesh', 'NzMesh',
                      'xmin', 'ymin', 'zmin', 'cfl']:
            setattr(config2, attr, getattr(self.config, attr))
        config2.field_propagator = "FDTD2"
        tile_fdtd2 = runko.emf.threeD.Tile((0, 0, 0), config2)
        tile_fdtd2.set_EBJ(E, B, J)
        tile_fdtd2.push_e()

        (Ex_s, Ey_s, Ez_s), _, _ = tile_stencil.get_EBJ()
        (Ex_f, Ey_f, Ez_f), _, _ = tile_fdtd2.get_EBJ()

        np.testing.assert_allclose(Ex_s, Ex_f, atol=1e-6)
        np.testing.assert_allclose(Ey_s, Ey_f, atol=1e-6)
        np.testing.assert_allclose(Ez_s, Ez_f, atol=1e-6)

    def test_nonzero_delta_linear_field(self):
        """
        With nonzero delta, a linear E field still has constant curl.
        The extended stencil of a linear field reduces to the same
        result as FDTD2 (since the higher-order finite differences
        of a linear function are zero).
        """
        self.config.stencil_delta = -0.065
        tile = self.make_test_tile()

        E = lambda x, y, z: (0, z, -y)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = Bx[3:-3, 3:-3, 3:-3].flat[0]
        self.assertTrue(b > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], b, places=5)
            self.assertAlmostEqual(By[i, j, k], 0, places=5)
            self.assertAlmostEqual(Bz[i, j, k], 0, places=5)

    def test_nonzero_beta_linear_field(self):
        """
        With nonzero beta transverse coefficients, a linear E field
        still gives uniform B in the interior (transverse averages of
        a constant field are still constant).
        """
        self.config.stencil_delta = -0.065
        self.config.stencil_beta_p1 = -0.065
        self.config.stencil_beta_p2 = -0.065
        tile = self.make_test_tile()

        E = lambda x, y, z: (0, z, -y)
        B = lambda x, y, z: (0, 0, 0)
        J = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(E, B, J)
        tile.push_half_b()

        _, (Bx, By, Bz), _ = tile.get_EBJ()

        b = Bx[3:-3, 3:-3, 3:-3].flat[0]
        self.assertTrue(b > 0)

        for i, j, k in self.get_interior_index_space():
            self.assertAlmostEqual(Bx[i, j, k], b, places=5)
            self.assertAlmostEqual(By[i, j, k], 0, places=5)
            self.assertAlmostEqual(Bz[i, j, k], 0, places=5)


if __name__ == "__main__":
    unittest.main()
