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
        self.config.xmin = -2.1
        self.config.ymin = 4.2
        self.config.zmin = 1.2
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

    def test_antenna_mode_trivial_time_evolution(self):
        """
        phi = 1 # One for the trivial time evolution
        A = (1, 2, 3) * exp(i * x * k) * exp(i * phi)
        B = curl(A) = (0, -3, 2) * i * k * exp(i * x * k) * exp(i * phi)
        J = cfl * curl(B) = cfl * (0, 2, 3) * k^2 * Re(exp(i * x * k) * exp(i * phi))
        """

        tile = self.make_test_tile()

        K = 0.06
        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), k=(K, 0, 0), lap_coeffs=[1, 1, 1])
        tile.register_antenna(antenna)
        z = lambda x, y, z: np.zeros_like(x)

        for _ in range(3):
            tile.batch_set_EBJ(z, z, z, z, z, z, z, z, z)
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


    def test_antenna_mode_flipping_time_evolution(self):
        """
        phi = exp(i * n * pi)
        A = (1, 2, 3) * exp(i * x * k) * exp(i * phi)
        B = curl(A) = (0, -3, 2) * i * k * exp(i * x * k) * exp(i * phi)
        J = cfl * curl(B) = cfl * (0, 2, 3) * k^2 * Re(exp(i * x * k) * exp(i * phi))
        """

        tile = self.make_test_tile()

        K = 0.06
        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), k=(K, 0, 0), lap_coeffs=[1, -2, 3])
        tile.register_antenna(antenna)
        z = lambda x, y, z: np.zeros_like(x)

        sign = [1, -2, 3]

        for M in range(3):
            tile.batch_set_EBJ(z, z, z, z, z, z, z, z, z)
            tile.deposit_antenna_current()

            _, _, (Jx, Jy, Jz) = tile.get_EBJ()

            cfl = self.config.cfl
            gcmap = tile.global_coordinate_map()

            for i, j, k in self.get_index_space():
                # xi, _, _ = gcmap((i + 0.5, j, k))
                xj, _, _ = gcmap((i, j + 0.5, k))
                xk, _, _ = gcmap((i, j, k + 0.5))

                self.assertAlmostEqual(Jx[i, j, k], 0)
                self.assertAlmostEqual(Jy[i, j, k], sign[M] * 2 * cfl * K**2 * np.cos(K * xj), places=5, msg=f"{(i, j, k)}")
                self.assertAlmostEqual(Jz[i, j, k], sign[M] * 3 * cfl * K**2 * np.cos(K * xk), places=5, msg=f"{(i, j, k)}")


    def test_antenna_mode_rotating_time_evolution(self):
        """
        A = (1, 2, 3) * exp(i * y * k) * exp(i * phi)
        B = curl(A) = (3, 0, -1) * i * k * exp(i * y * k) * exp(i * phi)
        J = cfl * curl(B) = cfl * (1, 0, 3) * k^2 * Re(exp(i * y * k) * exp(i * phi))
        """

        tile = self.make_test_tile()

        N = 3
        phi = 0.23
        p = 1.12
        q = p * complex(np.cos(phi), np.sin(phi))

        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), n=(0, N, 0), lap_coeffs=[q, q**2, q**3])
        tile.register_antenna(antenna)
        z = lambda x, y, z: np.zeros_like(x)

        for M in range(1, 4):
            tile.batch_set_EBJ(z, z, z, z, z, z, z, z, z)
            tile.deposit_antenna_current()

            _, _, (Jx, Jy, Jz) = tile.get_EBJ()

            cfl = self.config.cfl
            gcmap = tile.global_coordinate_map()


            K = 2 * np.pi * N / (self.config.Ny * self.config.NyMesh)

            for i, j, k in self.get_index_space():
                _, yi, _ = gcmap((i + 0.5, j, k))
                _, yj, _ = gcmap((i, j + 0.5, k))
                _, yk, _ = gcmap((i, j, k + 0.5))

                self.assertAlmostEqual(Jx[i, j, k], p**M * cfl * K**2 * np.cos(K * yi + M * phi), places=4, msg=f"{(i, j, k)}")
                self.assertAlmostEqual(Jy[i, j, k], 0)
                self.assertAlmostEqual(Jz[i, j, k], p**M * 3 * cfl * K**2 * np.cos(K * yk + M * phi), places=4, msg=f"{(i, j, k)}")


    def test_antenna_mode_throws_when_lap_coeffs_run_out(self):

        tile = self.make_test_tile()

        antenna = runko.emf.threeD.antenna_mode(A=(1, 2, 3), n=(0, 1, 0), lap_coeffs=[1, 1, 1])
        tile.register_antenna(antenna)

        with self.assertRaises(Exception):
            for _ in range(4):
                tile.deposit_antenna_current()


if __name__ == "__main__":
    unittest.main()
