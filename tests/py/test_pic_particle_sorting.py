# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import numpy as np
import runko

def make_test_tile():
    config = runko.Configuration(None)
    config.n_tiles = [4, 4, 4]
    config.n_cells_per_tile = [6, 7, 6]
    config.cfl = 1
    config.field_propagator = "fdtd2"
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"

    return config, runko.pic.threeD.Tile((2, 1, 3), config)


class pic_particle_sorting(unittest.TestCase):
    def test_empty_tile_can_be_sorted(self):
        conf, tile = make_test_tile()

        def assertNoParticles(ptype):
            pos_x, pos_y, pos_z = tile.get_positions(ptype)
            vel_x, vel_y, vel_z = tile.get_velocities(ptype)

            self.assertEqual(0, len(pos_x))
            self.assertEqual(0, len(pos_y))
            self.assertEqual(0, len(pos_z))
            self.assertEqual(0, len(vel_x))
            self.assertEqual(0, len(vel_y))
            self.assertEqual(0, len(vel_z))
            self.assertEqual(0, len(tile.get_ids(ptype)))

        assertNoParticles(0)
        assertNoParticles(1)
        tile.sort_particles()
        assertNoParticles(0)
        assertNoParticles(1)


    def test_sorted_particles_stay_same(self):
        conf, tile = make_test_tile()

        def pgen(x, y, z):
            return runko.pic.threeD.ParticleStateBatch(pos=(x, y, z), vel=(x/2, y/2, z/2))

        def pgen1(x, y, z):
            return runko.pic.threeD.ParticleStateBatch(pos=(x+0.1, y+0.1, z+0.1), vel=(x+0.2, y+0.2, z+0.2))

        tile.batch_inject_to_cells(0, pgen)
        tile.batch_inject_to_cells(0, pgen)
        tile.batch_inject_to_cells(0, pgen)
        tile.batch_inject_to_cells(1, pgen1)
        tile.batch_inject_to_cells(1, pgen1)

        def get_hashes():
            pos0 = tile.get_positions(0)
            vel0 = tile.get_velocities(0)
            ids0 = tile.get_ids(0)
            pos1 = tile.get_positions(1)
            vel1 = tile.get_velocities(1)
            ids1 = tile.get_ids(1)

            hashes = []
            for x in [*pos0, *vel0, ids0, *pos1, *vel1, ids1]:
                x.sort()
                hashes.append(hash(tuple(x)))

            return hashes

        hashes0 = get_hashes();


        self.assertEqual(hashes0, get_hashes())
        tile.sort_particles()
        self.assertEqual(hashes0, get_hashes())
        tile.sort_particles()
        tile.sort_particles()
        self.assertEqual(hashes0, get_hashes())


class pic_payload(unittest.TestCase):
    """species 2 carries a payload of two components; it follows each particle"""

    def make_tile(self):
        conf, _ = make_test_tile()
        conf.q2, conf.m2, conf.n_payload2 = -1, 1, 2
        return runko.pic.threeD.Tile((2, 1, 3), conf)

    def test_sort_keeps_the_payload_with_its_particle(self):
        tile = self.make_tile()
        rng = np.random.default_rng(1)
        n = 300
        lo = np.array([12.0, 7.0, 18.0])  # tile (2, 1, 3) of 6 x 7 x 6 cells
        pos = lo + rng.uniform(0.0, 1.0, (n, 3)) * [6, 7, 6]
        P = runko.pic.threeD.ParticleState
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)  # the sort needs the lattice
        tile.inject(2, [P(pos=tuple(p), vel=(0.0, 0.0, 0.0)) for p in pos], payload=[pos[:, 0], -pos[:, 1]])
        ids = np.asarray(tile.get_ids(2))
        tile.sort_particles()
        self.assertFalse(np.array_equal(np.asarray(tile.get_ids(2)), ids))  # reordered
        x, y, _ = tile.get_positions(2)
        np.testing.assert_array_equal(tile.get_payload(2), x)
        np.testing.assert_array_equal(tile.get_payload(2, 1), -np.asarray(y))

    def test_generators_carry_the_payload(self):
        tile = self.make_tile()
        def pgen(x, y, z):
            vel = [np.zeros(len(x))] * 3
            return runko.pic.threeD.ParticleStateBatch(pos=[x + 0.5, y, z], vel=vel, payload=[x + 0.5, 2 * y])
        # entries dropped outside the stripe take their payload with them
        tile.batch_inject_in_x_stripe(2, pgen, 13.0, 15.0)
        x, y, _ = tile.get_positions(2)
        self.assertTrue(np.all((np.asarray(x) >= 13.0) & (np.asarray(x) < 15.0)))
        np.testing.assert_array_equal(tile.get_payload(2), x)
        np.testing.assert_array_equal(tile.get_payload(2, 1), 2 * np.asarray(y))
        # one particle per cell, its payload as the second element of the pair
        P = runko.pic.threeD.ParticleState
        tile = self.make_tile()
        tile.inject_to_each_cell(2, lambda x, y, z: ([P(pos=(x, y, z), vel=(0, 0, 0))], [[x], [z]]))
        x, _, z = tile.get_positions(2)
        self.assertEqual(len(x), 6 * 7 * 6)
        np.testing.assert_array_equal(tile.get_payload(2), x)
        np.testing.assert_array_equal(tile.get_payload(2, 1), z)

    def test_a_payload_is_required_iff_the_species_has_one(self):
        tile = self.make_tile()
        P = runko.pic.threeD.ParticleState
        states = [P(pos=(13.0, 8.0, 19.0), vel=(0, 0, 0))] * 3
        for species, payload in ((2, None), (2, [1.0, 2.0, 3.0]), (2, [[1.0, 2.0]] * 2), (0, [1.0, 2.0, 3.0])):
            with self.assertRaises(RuntimeError):
                tile.inject(species, states, payload=payload)
        tile.inject(2, states, payload=[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        tile.inject(0, states)
        np.testing.assert_array_equal(tile.get_payload(2, 1), [4.0, 5.0, 6.0])
        with self.assertRaises(IndexError):
            tile.get_payload(0)


if __name__ == "__main__":
    _ = runko.actions.RuntimeInstance()
    unittest.main()
