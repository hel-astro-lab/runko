import unittest
import runko


def decimal_part(x: float):
    return x - int(x)


def make_test_tile():
    config = runko.Configuration(None)
    config.Nx = 4
    config.Ny = 4
    config.Nz = 4
    config.NxMesh = 10
    config.NyMesh = 11
    config.NzMesh = 13
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
    config.q0 = -1
    config.m0 = 1
    config.delgam = 1.0e-5
    config.temperature_ratio = 1.0
    config.sigma = 40
    config.c_omp = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"

    return config, runko.pic.Tile((0, 0, 0), config)


def in_middle_part(x, y, z, config):
    a = config.xmin + 3 < x
    b = config.xmin + config.NxMesh - 3 > x
    c = config.ymin + 3 < y
    d = config.ymin + config.NyMesh - 3 > y
    e = config.zmin + 3 < z
    f = config.zmin + config.NzMesh - 3 > z
    return a and b and c and d and e and f


class pic2_particle_pusher_boris(unittest.TestCase):

    def test_Ex_only(self):
        config, tile = make_test_tile()

        E = lambda x, y, z: (0.1, 0, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        def still_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

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
        config, tile = make_test_tile()

        E = lambda x, y, z: (0, 0.1, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        def still_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

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
        config, tile = make_test_tile()

        E = lambda x, y, z: (0, 0, 0.1)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(E, zero, zero)

        def still_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

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
        config, tile = make_test_tile()

        B = lambda x, y, z: (0.1, 0, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, B, zero)

        def moving_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.ParticleState(pos=(x, y, z), vel=(0, 0.1, 0))]

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
        config, tile = make_test_tile()

        B = lambda x, y, z: (0, 0.1, 0)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, B, zero)

        def moving_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.ParticleState(pos=(x, y, z), vel=(0.1, 0, 0))]

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
        config, tile = make_test_tile()

        B = lambda x, y, z: (0, 0, 0.1)
        zero = lambda x, y, z: (0, 0, 0)

        tile.set_EBJ(zero, B, zero)

        def moving_particle(x, y, z):
            # Add particles only to the middle parts,
            # because fields in non-halo regions might be uninitialized.
            if not in_middle_part(x, y, z, config):
                return []
            return [runko.ParticleState(pos=(x, y, z), vel=(0, 0.1, 0))]

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


if __name__ == "__main__":
    unittest.main()
