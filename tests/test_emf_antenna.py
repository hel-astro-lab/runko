import unittest
import itertools
import numpy as np

import runko

class emf_antenna(unittest.TestCase):

    def setUp(self):
        self.config = runko.Configuration(None)
        self.config.Nx = 4
        self.config.Ny = 4
        self.config.Nz = 4
        self.config.NxMesh = 40
        self.config.NyMesh = 40
        self.config.NzMesh = 40
        self.config.xmin = 0
        self.config.ymin = 0
        self.config.zmin = 0
        self.config.cfl = 0.45
        self.config.field_propagator = "FDTD2"

        self.zero_field = lambda x, y, z: np.zeros_like(x)

        def f():
            tile_grid_idx = (0, 0, 0)
            tile = runko.emf.threeD.Tile(tile_grid_idx, self.config)
            tile.batch_set_EBJ(self.zero_field, self.zero_field, self.zero_field,
                               self.zero_field, self.zero_field, self.zero_field,
                               self.zero_field, self.zero_field, self.zero_field)
            return tile


        self.make_test_tile = f

        # This has to be function, as itertools.product seems to be singel pass.
        self.get_index_space = lambda : itertools.product(range(self.config.NxMesh),
                                                          range(self.config.NyMesh),
                                                          range(self.config.NzMesh))

    def test_antenna_mode_wave_vector_x(self):
        """
        A = (1, 2, 3) * exp(i * x * k)
        B = curl(A) = (0, -3, 2) * i * k * exp(i * x * k)
        J = cfl * curl(B) = cfl * (0, 2, 3) * k^2 * Re(exp(i * x * k))
        """

        tile = self.make_test_tile()

        K = 0.06
        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), k=(K, 0, 0))

        tile.register_antenna(antenna)
        tile.deposit_antenna_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        cfl = self.config.cfl
        gcmap = tile.global_coordinate_map()

        for i, j, k in self.get_index_space():
            # xi, _, _ = gcmap((i + 0.5, j, k))
            xj, _, _ = gcmap((i, j + 0.5, k))
            xk, _, _ = gcmap((i, j, k + 0.5))

            self.assertAlmostEqual(Jx[i, j, k], 0)
            self.assertAlmostEqual(Jy[i, j, k], 2 * cfl * K**2 * np.cos(K * xj), places=5, msg=f"{(i, j, k)}")
            self.assertAlmostEqual(Jz[i, j, k], 3 * cfl * K**2 * np.cos(K * xk), places=5, msg=f"{(i, j, k)}")


    def test_antenna_mode_wave_vector_y(self):
        """
        A = (1, 2, 3) * exp(i * y * k)
        B = curl(A) = (3, 0, -1) * i * k * exp(i * y * k)
        J = cfl * curl(B) = cfl * (1, 0, 3) * k^2 * Re(exp(i * y * k))
        """

        tile = self.make_test_tile()

        K = 0.05
        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), k=(0, K, 0))

        tile.register_antenna(antenna)
        tile.deposit_antenna_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        cfl = self.config.cfl
        gcmap = tile.global_coordinate_map()

        for i, j, k in self.get_index_space():
            _, yi, _ = gcmap((i + 0.5, j, k))
            _, yj, _ = gcmap((i, j + 0.5, k))
            _, yk, _ = gcmap((i, j, k + 0.5))

            self.assertAlmostEqual(Jx[i, j, k], cfl * K**2 * np.cos(K * yi), places=5, msg=f"{(i, j, k)}")
            self.assertAlmostEqual(Jy[i, j, k], 0)
            self.assertAlmostEqual(Jz[i, j, k], 3 * cfl * K**2 * np.cos(K * yk), places=5, msg=f"{(i, j, k)}")


    def test_antenna_mode_wave_vector_z(self):
        """
        A = (1, 2, 3) * exp(i * z * k)
        B = curl(A) = (-2, 1, 0) * i * k * exp(i * y * k)
        J = cfl * curl(B) = cfl * (1, 2, 0) * k^2 * Re(exp(i * y * k))
        """

        tile = self.make_test_tile()

        K = 0.07
        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), k=(0, 0, K))

        tile.register_antenna(antenna)
        tile.deposit_antenna_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        cfl = self.config.cfl
        gcmap = tile.global_coordinate_map()

        for i, j, k in self.get_index_space():
            _, _, zi = gcmap((i + 0.5, j, k))
            _, _, zj = gcmap((i, j + 0.5, k))

            self.assertAlmostEqual(Jx[i, j, k], cfl * K**2 * np.cos(K * zi), places=5, msg=f"{(i, j, k)}")
            self.assertAlmostEqual(Jy[i, j, k], 2 * cfl * K**2 * np.cos(K * zj), places=5, msg=f"{(i, j, k)}")
            self.assertAlmostEqual(Jz[i, j, k], 0)


    def test_antenna_mode_wave_vector_xyz(self):
        """
        k = (kx, ky, kz)
        A = (1, 0, 0) * exp(i * dot(k, x))
        B = curl(A) = (0, kz, -ky) * i * exp(i * dot(k, x))
        J = cfl * curl(B) = cfl * (ky^2 + kz^2, -kx * ky, -kx * kz) * Re(exp(i * dot(k, x)))
        """

        tile = self.make_test_tile()

        kx, ky, kz = 0.13, 0.1, 0.08
        antenna = runko.emf.threeD.antenna_mode(A=(1, 0, 0), k=(kx, ky, kz))

        tile.register_antenna(antenna)
        tile.deposit_antenna_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        cfl = self.config.cfl
        gcmap = tile.global_coordinate_map()

        for i, j, k in self.get_index_space():
            xi, yi, zi = gcmap((i + 0.5, j, k))
            xj, yj, zj = gcmap((i, j + 0.5, k))
            xk, yk, zk = gcmap((i, j, k + 0.5))

            phi_i = xi * kx + yi * ky + zi * kz
            phi_j = xj * kx + yj * ky + zj * kz
            phi_k = xk * kx + yk * ky + zk * kz

            self.assertAlmostEqual(Jx[i, j, k], cfl * (ky**2 + kz**2) * np.cos(phi_i), places=2, msg=f"{(i, j, k)}")
            self.assertAlmostEqual(Jy[i, j, k], cfl * -kx * ky * np.cos(phi_j), places=2, msg=f"{(i, j, k)}")
            self.assertAlmostEqual(Jz[i, j, k], cfl * -kx * kz * np.cos(phi_k), places=2, msg=f"{(i, j, k)}")


if __name__ == "__main__":
    unittest.main()
