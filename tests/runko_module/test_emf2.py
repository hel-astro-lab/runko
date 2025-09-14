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
        config.field_propagator = "FDTD2"

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
        config.field_propagator = "FDTD2"

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

        # Make sure numpy mdarrays own the data.
        self.assertIsNone(Ex.base)
        self.assertIsNone(Ey.base)
        self.assertIsNone(Ez.base)
        self.assertIsNone(Bx.base)
        self.assertIsNone(By.base)
        self.assertIsNone(Bz.base)
        self.assertIsNone(Jx.base)
        self.assertIsNone(Jy.base)
        self.assertIsNone(Jz.base)

        s = (config.NxMesh, config.NyMesh, config.NzMesh)
        self.assertEqual(Ex.shape, s)
        self.assertEqual(Ey.shape, s)
        self.assertEqual(Ez.shape, s)
        self.assertEqual(Bx.shape, s)
        self.assertEqual(By.shape, s)
        self.assertEqual(Bz.shape, s)
        self.assertEqual(Jx.shape, s)
        self.assertEqual(Jy.shape, s)
        self.assertEqual(Jz.shape, s)

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


    def test_field_batch_set(self):

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
        config.field_propagator = "FDTD2"

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


        Exinit = lambda x, y, z: x
        Eyinit = lambda x, y, z: y
        Ezinit = lambda x, y, z: z
        Bxinit = lambda x, y, z: 2 * x
        Byinit = lambda x, y, z: 2 * y
        Bzinit = lambda x, y, z: 2 * z
        Jxinit = lambda x, y, z: 3 * x
        Jyinit = lambda x, y, z: 3 * y
        Jzinit = lambda x, y, z: 3 * z

        tile.batch_set_EBJ(Exinit, Eyinit, Ezinit,
                           Bxinit, Byinit, Bzinit,
                           Jxinit, Jyinit, Jzinit)

        (Ex, Ey, Ez), (Bx, By, Bz), (Jx, Jy, Jz) = tile.get_EBJ()

        index_space = itertools.product(range(config.NxMesh),
                                        range(config.NyMesh),
                                        range(config.NzMesh))
        for i, j, k in index_space:

            ii = config.xmin + tile_grid_idx[0] * config.NxMesh + i
            jj = config.ymin + tile_grid_idx[1] * config.NyMesh + j
            kk = config.zmin + tile_grid_idx[2] * config.NzMesh + k

            self.assertAlmostEqual(Ex[i, j, k], Exinit(ii + 0.5, jj, kk), places=5)
            self.assertAlmostEqual(Ey[i, j, k], Eyinit(ii, jj + 0.5, kk), places=5)
            self.assertAlmostEqual(Ez[i, j, k], Ezinit(ii, jj, kk + 0.5), places=5)

            self.assertAlmostEqual(Bx[i, j, k], Bxinit(ii, jj + 0.5, kk + 0.5), places=5)
            self.assertAlmostEqual(By[i, j, k], Byinit(ii + 0.5, jj, kk + 0.5), places=5)
            self.assertAlmostEqual(Bz[i, j, k], Bzinit(ii + 0.5, jj + 0.5, kk), places=5)

            self.assertAlmostEqual(Jx[i, j, k], Jxinit(ii + 0.5, jj, kk), places=5)
            self.assertAlmostEqual(Jy[i, j, k], Jyinit(ii, jj + 0.5, kk), places=5)
            self.assertAlmostEqual(Jz[i, j, k], Jzinit(ii, jj, kk + 0.5), places=5)


    def test_batch_setting_fields_with_incorrect_shape_throws(self):

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
        config.field_propagator = "FDTD2"

        tile_grid_idx = (0, 1, 2)
        tile = runko.emf.Tile(tile_grid_idx, config)

        f0 = lambda x, y, z: np.ones((1))
        f1 = lambda x, y, z: np.ones((1, 1))
        f2 = lambda x, y, z: np.ones((1, 1, 1))
        f3 = lambda x, y, z: np.ones((10, 12, 15))
        f4 = lambda x, y, z: np.ones((10, 11, 14))
        f5 = lambda x, y, z: np.ones((11, 12, 14))
        f6 = lambda x, y, z: np.ones((1, 1, 1, 1))
        f7 = lambda x, y, z: np.ones((1, 1, 1, 1, 1))
        f7 = lambda x, y, z: np.ones(())
        fc = lambda x, y, z: np.ones((10, 12, 14))

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(f0, fc, fc,
                               fc, fc, fc,
                               fc, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, f1, fc,
                               fc, fc, fc,
                               fc, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, f2,
                               fc, fc, fc,
                               fc, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, fc,
                               f3, fc, fc,
                               fc, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, fc,
                               fc, f4, fc,
                               fc, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, fc,
                               fc, fc, f5,
                               fc, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, fc,
                               fc, fc, fc,
                               f6, fc, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, fc,
                               fc, fc, fc,
                               fc, f7, fc)

        with self.assertRaises(Exception):
            tile.batch_set_EBJ(fc, fc, fc,
                               fc, fc, fc,
                               fc, fc, f8)

        tile.batch_set_EBJ(fc, fc, fc,
                           fc, fc, fc,
                           fc, fc, fc)

        (Ex, Ey, Ez), (Bx, By, Bz), (Jx, Jy, Jz) = tile.get_EBJ()
        self.assertTrue(np.all(Ex == 1))
        self.assertTrue(np.all(Ey == 1))
        self.assertTrue(np.all(Ez == 1))
        self.assertTrue(np.all(Bx == 1))
        self.assertTrue(np.all(By == 1))
        self.assertTrue(np.all(Bz == 1))
        self.assertTrue(np.all(Jx == 1))
        self.assertTrue(np.all(Jy == 1))
        self.assertTrue(np.all(Jz == 1))


    def test_subtract_J_from_E(self):

        config = runko.Configuration(None)
        config.Nx = 2
        config.Ny = 3
        config.Nz = 4
        config.NxMesh = 10
        config.NyMesh = 12
        config.NzMesh = 14
        config.xmin = 0
        config.ymin = 0
        config.zmin = 0
        config.cfl = 1
        config.field_propagator = "FDTD2"

        tile_grid_idx = (0, 0, 0)
        tile = runko.emf.Tile(tile_grid_idx, config)

        Einit = lambda x, y, z: (0, 0, 0)
        Binit = lambda x, y, z: (0, 0, 0)
        Jinit = lambda x, y, z: (1 + y, 2 + z, 3 + x)

        tile.set_EBJ(Einit, Binit, Jinit)

        (E0x, E0y, E0z), _, _ = tile.get_EBJ()

        self.assertTrue(np.all(E0x == 0))
        self.assertTrue(np.all(E0y == 0))
        self.assertTrue(np.all(E0z == 0))

        tile.subtract_J_from_E()

        # Now E = -A * J, for some positive scalar A.
        # We can calculate A from each component and test that they are eqal.
        # Technically due to the units used in runko A should be 1.
        # However, I don't want to assume that in this test.

        (Ex, Ey, Ez), _, (Jx, Jy, Jz) = tile.get_EBJ()

        Ax_arr, Ay_arr, Az_arr = Ex / Jx, Ey / Jy, Ez / Jz
        Ax, Ay, Az = Ax_arr.flat[0], Ay_arr.flat[0], Az_arr.flat[0]

        index_space = itertools.product(range(config.NxMesh),
                                        range(config.NyMesh),
                                        range(config.NzMesh))

        for i, j, k in index_space:
            self.assertAlmostEqual(Ax, Ax_arr[i, j, k], places=5)
            self.assertAlmostEqual(Ay, Ay_arr[i, j, k], places=5)
            self.assertAlmostEqual(Az, Az_arr[i, j, k], places=5)
            self.assertTrue(Ax_arr[i, j, k] < 0)
            self.assertTrue(Ay_arr[i, j, k] < 0)
            self.assertTrue(Az_arr[i, j, k] < 0)


    def test_tile_index_is_set(self):
        config = runko.Configuration(None)
        config.Nx = 2
        config.Ny = 3
        config.Nz = 4
        config.NxMesh = 10
        config.NyMesh = 12
        config.NzMesh = 14
        config.xmin = 0
        config.ymin = 0
        config.zmin = 0
        config.cfl = 1
        config.field_propagator = "FDTD2"

        tile_grid_idx = [1, 2, 3]
        tile = runko.emf.Tile(tile_grid_idx, config)

        """
        This is not 100% neccesseary, as corgi grids add_tile method
        will set the index too. However, many tile methods assume that
        it is set, so this prevents bugs when tiles are used outside
        of corgi grid.
        """

        self.assertEqual(tile.index, tile_grid_idx)


    def test_bogus_tile_index_raises_exception(self):
        config = runko.Configuration(None)
        config.Nx = 2
        config.Ny = 3
        config.Nz = 4
        config.NxMesh = 10
        config.NyMesh = 12
        config.NzMesh = 14
        config.xmin = 0
        config.ymin = 0
        config.zmin = 0
        config.cfl = 1
        config.field_propagator = "FDTD2"

        with self.assertRaises(Exception):
            runko.emf.Tile((2, 2, 3), config)

        with self.assertRaises(Exception):
            runko.emf.Tile((1, 3, 3), config)

        with self.assertRaises(Exception):
            runko.emf.Tile((1, 2, 4), config)

        with self.assertRaises(Exception):
            runko.emf.Tile((-1, 0, 0), config)

        with self.assertRaises(Exception):
            runko.emf.Tile((0, -1, 0), config)

        with self.assertRaises(Exception):
            runko.emf.Tile((0, 0, -1), config)

if __name__ == "__main__":
    unittest.main()
