# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import itertools
import numpy as np
import runko


def make_test_tile():
    config = runko.Configuration(None)
    config.n_tiles = [4, 4, 4]
    config.n_cells_per_tile = [10, 11, 13]
    config.cfl = 1
    config.field_propagator = "fdtd2"
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.q2 = 2
    config.m2 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"

    return config, runko.pic.threeD.Tile((0, 0, 0), config)


def p_func(v):
    def moving_particle(x, y, z):
        base = np.array((x + 0.5, y + 0.5, z + 0.5))
        dloc = 0.1 * np.array((np.sin(x), np.cos(y), np.sin(x + z)))

        return [runko.pic.threeD.ParticleState(pos=base + dloc, vel=v)]

    return moving_particle


def interior_index_space(config):
    I = range(1, config.n_cells_per_tile[0] - 1)
    J = range(1, config.n_cells_per_tile[1] - 1)
    K = range(1, config.n_cells_per_tile[2] - 1)
    return itertools.product(I, J, K)


class pic_currend_depositer_zigzag_1st(unittest.TestCase):

    def test_nonexisting_particles_do_not_create_current(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.all(Jx == 0))
        self.assertTrue(np.all(Jy == 0))
        self.assertTrue(np.all(Jz == 0))


    def test_unmoving_particles_do_not_create_current(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        tile.inject_to_each_cell(0, p_func((0, 0, 0)))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.all(Jx == 0))
        self.assertTrue(np.all(Jy == 0))
        self.assertTrue(np.all(Jz == 0))


    def test_Jx_only(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        tile.inject_to_each_cell(0, p_func((0.1, 0, 0)))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.all(Jx < 0))
        self.assertTrue(np.all(Jy == 0))
        self.assertTrue(np.all(Jz == 0))


    def test_Jy_only(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        tile.inject_to_each_cell(0, p_func((0, 0.1, 0)))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.all(Jx == 0))
        self.assertTrue(np.all(Jy < 0))
        self.assertTrue(np.all(Jz == 0))


    def test_Jz_only(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        tile.inject_to_each_cell(0, p_func((0, 0, 0.1)))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.all(Jx == 0))
        self.assertTrue(np.all(Jy == 0))
        self.assertTrue(np.all(Jz < 0))


    def test_particles_moving_in_opposite_directions_cancel_each_other(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        v = np.array((0.2, 0.1, -0.09))
        tile.inject_to_each_cell(1, p_func(-v))
        tile.inject_to_each_cell(1, p_func(v))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        # Without virtual tiles the boundaries can not be tested.
        for idx in interior_index_space(config):
            self.assertAlmostEqual(Jx[*idx], 0, places=2)
            self.assertAlmostEqual(Jy[*idx], 0, places=2)
            self.assertAlmostEqual(Jz[*idx], 0, places=2)


    def test_particles_opposite_charges_cancel_each_other(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        v = np.array((-0.2, 0.1, 0.09))
        tile.inject_to_each_cell(0, p_func(v))
        tile.inject_to_each_cell(1, p_func(v))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        # Without virtual tiles the boundaries can not be tested.
        for idx in interior_index_space(config):
            self.assertAlmostEqual(Jx[*idx], 0, places=2)
            self.assertAlmostEqual(Jy[*idx], 0, places=2)
            self.assertAlmostEqual(Jz[*idx], 0, places=2)


    def test_three_particles_cancel_each_other(self):
        config, tile = make_test_tile()

        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, zero, zero)

        v = np.array((0.1, -0.2, -0.03))
        tile.inject_to_each_cell(0, p_func(v))
        tile.inject_to_each_cell(1, p_func(-v))
        tile.inject_to_each_cell(2, p_func(v))

        _, _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        self.assertTrue(np.all(Jx0 == 0))
        self.assertTrue(np.all(Jy0 == 0))
        self.assertTrue(np.all(Jz0 == 0))

        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        # Without virtual tiles the boundaries can not be tested.
        for idx in interior_index_space(config):
            self.assertAlmostEqual(Jx[*idx], 0, places=2)
            self.assertAlmostEqual(Jy[*idx], 0, places=2)
            self.assertAlmostEqual(Jz[*idx], 0, places=2)


if __name__ == "__main__":
    unittest.main()
