import unittest
import runko

def make_valid_emf2_config():
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

    return config


class pic2_tile(unittest.TestCase):

    def test_pic_tiles_are_initially_empty(self):

        config = make_valid_emf2_config()
        config.qe = 1
        config.me = 1
        config.delgam = 1.0e-5
        config.temperature_ratio = 1.0
        config.sigma = 40
        config.c_omp = 1
        config.particle_pusher = "boris"
        config.field_interpolator = "linear_1st"

        tile_grid_idx = (0, 1, 2)
        tile = runko.pic.Tile(tile_grid_idx, config)

        (pos0_x, pos0_y, pos0_z) = tile.get_positions(runko.particle.electron)
        (vel0_x, vel0_y, vel0_z) = tile.get_velocities(runko.particle.electron)

        self.assertEqual(0, len(pos0_x))
        self.assertEqual(0, len(pos0_y))
        self.assertEqual(0, len(pos0_z))
        self.assertEqual(0, len(vel0_x))
        self.assertEqual(0, len(vel0_y))
        self.assertEqual(0, len(vel0_z))


    def test_inject_to_each_cell_roundtrip(self):

        config = make_valid_emf2_config()
        config.qe = 1
        config.me = 1
        config.delgam = 1.0e-5
        config.temperature_ratio = 1.0
        config.sigma = 40
        config.c_omp = 1
        config.particle_pusher = "boris"
        config.field_interpolator = "linear_1st"

        tile_grid_idx = (0, 0, 0)
        tile = runko.pic.Tile(tile_grid_idx, config)

        ppc = 2

        # This should be called for each (non-halo) cell in the tile.
        def particle_generator(x, y, z):
            P = runko.ParticleState
            new_p = []
            for i in range(ppc):
                new_p.append( P(pos=(x, y, z), vel=(0, 0, i)))
            return new_p

        tile.inject_to_each_cell(runko.particle.electron, particle_generator)

        pos_x, pos_y, pos_z = tile.get_positions(runko.particle.electron)
        vel_x, vel_y, vel_z = tile.get_velocities(runko.particle.electron)

        expected_num_of_particles = ppc * config.NxMesh * config.NyMesh * config.NzMesh

        self.assertEqual(expected_num_of_particles, len(pos_x))
        self.assertEqual(expected_num_of_particles, len(pos_y))
        self.assertEqual(expected_num_of_particles, len(pos_z))
        self.assertEqual(expected_num_of_particles, len(vel_x))
        self.assertEqual(expected_num_of_particles, len(vel_y))
        self.assertEqual(expected_num_of_particles, len(vel_z))

        # Go through all particles, check their velocities
        # and store their positions if vz == 0.
        # Then go through all cells and make sure that every cell is present.

        seen_pos = set()

        for i in range(expected_num_of_particles):
            x, y, z = pos_x[i], pos_y[i], pos_z[i]
            vx, vy, vz = vel_x[i], vel_y[i], vel_z[i]

            if vz == 0:
                self.assertTrue((x, y, z) not in seen_pos)
                seen_pos.add((x, y, z))
            else:
                self.assertEqual(vz, 1)

            self.assertEqual(vx, 0)
            self.assertEqual(vy, 0)

        self.assertEqual(config.NxMesh * config.NyMesh * config.NzMesh, len(seen_pos))

        import itertools
        cell_index_space = itertools.product(range(config.NxMesh),
                                             range(config.NyMesh),
                                             range(config.NzMesh))
        for i, j, k in cell_index_space:
            self.assertTrue((i, j, k) in seen_pos)


    def test_inject_to_each_cell_multiple_times(self):

        config = make_valid_emf2_config()
        config.qe = 1
        config.me = 1
        config.delgam = 1.0e-5
        config.temperature_ratio = 1.0
        config.sigma = 40
        config.c_omp = 1
        config.particle_pusher = "boris"
        config.field_interpolator = "linear_1st"

        tile_grid_idx = (0, 0, 0)
        tile = runko.pic.Tile(tile_grid_idx, config)

        # This should be called for each (non-halo) cell in the tile.
        def make_gen(w):
            def particle_generator(x, y, z):
                return [runko.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]
            return particle_generator

        def assertLengths(expected_num_of_particles):
            pos_x, pos_y, pos_z = tile.get_positions(runko.particle.electron)
            vel_x, vel_y, vel_z = tile.get_velocities(runko.particle.electron)

            self.assertEqual(expected_num_of_particles, len(pos_x))
            self.assertEqual(expected_num_of_particles, len(pos_y))
            self.assertEqual(expected_num_of_particles, len(pos_z))
            self.assertEqual(expected_num_of_particles, len(vel_x))
            self.assertEqual(expected_num_of_particles, len(vel_y))
            self.assertEqual(expected_num_of_particles, len(vel_z))


        assertLengths(0)

        tile.inject_to_each_cell(runko.particle.electron, make_gen(1))
        tile.inject_to_each_cell(runko.particle.electron, make_gen(2))
        tile.inject_to_each_cell(runko.particle.electron, make_gen(3))

        N = config.NxMesh * config.NyMesh * config.NzMesh
        assertLengths(3 * N)


    def test_inject_to_each_cell_multiple_particle_types(self):

        config = make_valid_emf2_config()
        config.qe = 1
        config.me = 1
        config.qi = 1
        config.mi = 1
        config.delgam = 1.0e-5
        config.temperature_ratio = 1.0
        config.sigma = 40
        config.c_omp = 1
        config.particle_pusher = "boris"
        config.field_interpolator = "linear_1st"

        tile_grid_idx = (0, 0, 0)
        tile = runko.pic.Tile(tile_grid_idx, config)

        def assertLengths(type, expected_num_of_particles):
            pos_x, pos_y, pos_z = tile.get_positions(type)
            vel_x, vel_y, vel_z = tile.get_velocities(type)

            self.assertEqual(expected_num_of_particles, len(pos_x))
            self.assertEqual(expected_num_of_particles, len(pos_y))
            self.assertEqual(expected_num_of_particles, len(pos_z))
            self.assertEqual(expected_num_of_particles, len(vel_x))
            self.assertEqual(expected_num_of_particles, len(vel_y))
            self.assertEqual(expected_num_of_particles, len(vel_z))


        assertLengths(runko.particle.electron, 0)
        assertLengths(runko.particle.ion, 0)

        def particle_generator_electron(x, y, z):
            return [runko.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

        def particle_generator_ion(x, y, z):
            return 2 * [runko.ParticleState(pos=(x, y, z), vel=(0, 0, 0))]

        tile.inject_to_each_cell(runko.particle.electron, particle_generator_electron)
        tile.inject_to_each_cell(runko.particle.ion, particle_generator_ion)

        N = config.NxMesh * config.NyMesh * config.NzMesh
        assertLengths(runko.particle.electron, N)
        assertLengths(runko.particle.ion, 2 * N)


if __name__ == "__main__":
    unittest.main()
