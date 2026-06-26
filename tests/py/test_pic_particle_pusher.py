# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import runko


def decimal_part(x: float):
    return x - int(x)


def make_test_tile(particle_pusher, field_interpolator):
    config = runko.Configuration(None)
    config.n_tiles = [4, 4, 4]
    config.n_cells_per_tile = [10, 11, 13]
    config.cfl = 1
    config.field_propagator = "fdtd2"
    config.q0 = -1
    config.m0 = 1
    config.particle_pusher = particle_pusher
    config.field_interpolator = field_interpolator
    config.current_depositer = "zigzag_1st"

    return config, runko.pic.threeD.Tile((3, 3, 3), config)


def in_middle_part(x, y, z, config):
    a = 3 < x
    b = config.n_cells_per_tile[0] - 3 > x
    c = 3 < y
    d = config.n_cells_per_tile[1] - 3 > y
    e = 3 < z
    f = config.n_cells_per_tile[2] - 3 > z
    return a and b and c and d and e and f


class _particle_pusher_tests:

    particle_pusher = None     # override in subclass
    field_interpolator = None  # override in subclass

    def test_Ex_only(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        E = lambda x, y, z: (0.1, 0, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        def still_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

        tile.inject_to_each_cell(0, still_particle)

        pos0 = tile.get_positions(0)
        vel0 = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos0, *vel0):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0)
            self.assertEqual(vz, 0)

        tile.push_particles()

        pos = tile.get_positions(0)
        vel = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos, *vel):
            self.assertNotEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertTrue(vx < 0)
            self.assertEqual(vy, 0)
            self.assertEqual(vz, 0)


    def test_Ey_only(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        E = lambda x, y, z: (0, 0.1, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        def still_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

        tile.inject_to_each_cell(0, still_particle)

        pos0 = tile.get_positions(0)
        vel0 = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos0, *vel0):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0)
            self.assertEqual(vz, 0)

        tile.push_particles()

        pos = tile.get_positions(0)
        vel = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos, *vel):
            self.assertEqual(decimal_part(x), 0)
            self.assertNotEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertTrue(vy < 0)
            self.assertEqual(vz, 0)


    def test_Ez_only(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        E = lambda x, y, z: (0, 0, 0.1)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        def still_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

        tile.inject_to_each_cell(0, still_particle)

        pos0 = tile.get_positions(0)
        vel0 = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos0, *vel0):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0)
            self.assertEqual(vz, 0)

        tile.push_particles()

        pos = tile.get_positions(0)
        vel = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos, *vel):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertNotEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0)
            self.assertTrue(vz < 0)


    def test_Bx_only(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        B = lambda x, y, z: (0.1, 0, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, B, zero)

        def moving_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(0, 0.1, 0))]

        tile.inject_to_each_cell(0, moving_particle)

        pos0 = tile.get_positions(0)
        vel0 = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos0, *vel0):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0.1)
            self.assertEqual(vz, 0)

        tile.push_particles()

        pos = tile.get_positions(0)
        vel = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos, *vel):
            self.assertEqual(decimal_part(x), 0)
            self.assertNotEqual(decimal_part(y), 0)
            self.assertNotEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertAlmostEqual(vy, 0.1, places=3)
            self.assertTrue(vz > 0)


    def test_By_only(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        B = lambda x, y, z: (0, 0.1, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, B, zero)

        def moving_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(0.1, 0, 0))]

        tile.inject_to_each_cell(0, moving_particle)

        pos0 = tile.get_positions(0)
        vel0 = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos0, *vel0):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0.1)
            self.assertEqual(vy, 0)
            self.assertEqual(vz, 0)

        tile.push_particles()

        pos = tile.get_positions(0)
        vel = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos, *vel):
            self.assertNotEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertNotEqual(decimal_part(z), 0)

            self.assertAlmostEqual(vx, 0.1, places=3)
            self.assertEqual(vy, 0)
            self.assertTrue(vz < 0)


    def test_Bz_only(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        B = lambda x, y, z: (0, 0, 0.1)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, B, zero)

        def moving_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(0, 0.1, 0))]

        tile.inject_to_each_cell(0, moving_particle)

        pos0 = tile.get_positions(0)
        vel0 = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos0, *vel0):
            self.assertEqual(decimal_part(x), 0)
            self.assertEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0.1)
            self.assertEqual(vz, 0)

        tile.push_particles()

        pos = tile.get_positions(0)
        vel = tile.get_velocities(0)

        for x, y, z, vx, vy, vz in zip(*pos, *vel):
            self.assertNotEqual(decimal_part(x), 0)
            self.assertNotEqual(decimal_part(y), 0)
            self.assertEqual(decimal_part(z), 0)

            self.assertTrue(vx < 0)
            self.assertAlmostEqual(vy, 0.1, places=3)
            self.assertEqual(vz, 0)


    def test_same_E_at_different_x_has_the_same_effect(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        E = lambda x, y, z: (0, 0.1, 0.1)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        i, j, k = tile.index
        xA = (i + 0.1) * config.n_cells_per_tile[0]
        xB = (i + 0.9) * config.n_cells_per_tile[0]
        y = (j + 0.5) * config.n_cells_per_tile[1]
        z = (k + 0.5) * config.n_cells_per_tile[2]

        particles = [runko.pic.threeD.ParticleState(pos=(xA, y, z), vel=(0, 0, 0)),
                     runko.pic.threeD.ParticleState(pos=(xB, y, z), vel=(0, 0, 0))]

        tile.inject(0, particles)

        def assertSameYZ():
            _, (yA, yB), (zA, zB) = tile.get_positions(0)
            self.assertAlmostEqual(yA, yB, places=5)
            self.assertAlmostEqual(zA, zB, places=5)

        assertSameYZ()
        tile.push_particles()
        assertSameYZ()
        tile.push_particles()
        assertSameYZ()


    def test_same_E_at_different_y_has_the_same_effect(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        E = lambda x, y, z: (0.1, 0, 0.1)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        i, j, k = tile.index
        x = (i + 0.5) * config.n_cells_per_tile[0]
        yA = (j + 0.1) * config.n_cells_per_tile[1]
        yB = (j + 0.9) * config.n_cells_per_tile[1]
        z = (k + 0.5) * config.n_cells_per_tile[2]

        particles = [runko.pic.threeD.ParticleState(pos=(x, yA, z), vel=(0, 0, 0)),
                     runko.pic.threeD.ParticleState(pos=(x, yB, z), vel=(0, 0, 0))]

        tile.inject(0, particles)

        def assertSameXZ():
            (xA, xB), _, (zA, zB) = tile.get_positions(0)
            self.assertAlmostEqual(xA, xB, places=5)
            self.assertAlmostEqual(zA, zB, places=5)

        assertSameXZ()
        tile.push_particles()
        assertSameXZ()
        tile.push_particles()
        assertSameXZ()


    def test_same_E_at_different_z_has_the_same_effect(self):
        config, tile = make_test_tile(self.particle_pusher, self.field_interpolator)

        E = lambda x, y, z: (0.1, 0.1, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        i, j, k = tile.index
        x = (i + 0.5) * config.n_cells_per_tile[0]
        y = (j + 0.5) * config.n_cells_per_tile[1]
        zA = (k + 0.1) * config.n_cells_per_tile[2]
        zB = (k + 0.9) * config.n_cells_per_tile[2]

        particles = [runko.pic.threeD.ParticleState(pos=(x, y, zA), vel=(0, 0, 0)),
                     runko.pic.threeD.ParticleState(pos=(x, y, zB), vel=(0, 0, 0))]

        tile.inject(0, particles)

        def assertSameXY():
            (xA, xB), (yA, yB), _ = tile.get_positions(0)
            self.assertAlmostEqual(xA, xB, places=5)
            self.assertAlmostEqual(yA, yB, places=5)

        assertSameXY()
        tile.push_particles()
        assertSameXY()
        tile.push_particles()
        assertSameXY()


# Generate concrete test classes for all (pusher, interpolator) combinations
_pushers = ["boris", "higuera_cary", "faraday"]
_interpolators = ["linear_1st", "linear_1st_unrolled"]

for _pusher in _pushers:
    for _interp in _interpolators:
        _name = f"pic_particle_pusher_{_pusher}_{_interp}"
        globals()[_name] = type(_name, (_particle_pusher_tests, unittest.TestCase), {
            "particle_pusher": _pusher,
            "field_interpolator": _interp,
        })


if __name__ == "__main__":
    unittest.main()
