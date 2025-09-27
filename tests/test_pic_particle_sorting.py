import unittest
import runko

def make_test_tile():
    config = runko.Configuration(None)
    config.Nx = 4
    config.Ny = 4
    config.Nz = 4
    config.NxMesh = 6
    config.NyMesh = 7
    config.NzMesh = 6
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
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

        assertNoParticles(0)
        assertNoParticles(1)
        tile.sort_particles()
        assertNoParticles(0)
        assertNoParticles(1)


    def test_sorted_particles_stay_same(self):
        conf, tile = make_test_tile()

        def pgen(x, y, z):
            return runko.pic.threeD.ParticleStateBatch(pos=(x, y, z), vel=(x, y, z))

        tile.batch_inject_to_cells(0, pgen)
        tile.batch_inject_to_cells(1, pgen)

        pos0set0 = set((x, y, z) for x, y, z in zip(*tile.get_positions(0)))
        vel0set0 = set((x, y, z) for x, y, z in zip(*tile.get_velocities(0)))
        pos0set1 = set((x, y, z) for x, y, z in zip(*tile.get_positions(1)))
        vel0set1 = set((x, y, z) for x, y, z in zip(*tile.get_velocities(1)))

        def assertNoChange():
            posset0 = set((x, y, z) for x, y, z in zip(*tile.get_positions(0)))
            velset0 = set((x, y, z) for x, y, z in zip(*tile.get_velocities(0)))
            posset1 = set((x, y, z) for x, y, z in zip(*tile.get_positions(1)))
            velset1 = set((x, y, z) for x, y, z in zip(*tile.get_velocities(1)))

            self.assertEqual(posset0, pos0set0)
            self.assertEqual(velset0, vel0set0)
            self.assertEqual(posset1, pos0set1)
            self.assertEqual(velset1, vel0set1)

        assertNoChange()
        tile.sort_particles()
        assertNoChange()
        tile.sort_particles()
        tile.sort_particles()
        assertNoChange()


if __name__ == "__main__":
    unittest.main()
