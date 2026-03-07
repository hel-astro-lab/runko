import unittest
import numpy as np
import runko


zero = lambda x, y, z: (0, 0, 0)


def make_reflector_tile():
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
    config.delgam = 1e-5
    config.temperature_ratio = 1.0
    config.sigma = 40
    config.c_omp = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"

    return config, runko.pic.threeD.Tile((0, 0, 0), config)


class pic_reflector_wall(unittest.TestCase):

    def test_stationary_wall_reflects_ux(self):
        _, tile = make_reflector_tile()

        wall = runko.pic.threeD.reflector_wall(walloc=5.0)
        tile.register_reflector_wall(wall)

        tile.set_EBJ(zero, zero, zero)

        P = runko.pic.threeD.ParticleState
        tile.inject(0, [P(pos=(5.5, 5.5, 6.5), vel=(-1.0, 0.0, 0.0))])

        tile.push_particles()
        tile.reflect_particles()

        pos_x, pos_y, pos_z = tile.get_positions(0)
        vel_x, vel_y, vel_z = tile.get_velocities(0)

        self.assertEqual(1, len(pos_x))

        # ux must have flipped to positive
        self.assertGreater(vel_x[0], 0)

        # magnitude preserved for stationary wall
        self.assertAlmostEqual(abs(vel_x[0]), 1.0, places=4)

        # particle must be back in front of the wall
        self.assertGreater(pos_x[0], 5.0)

        # transverse velocity unchanged
        self.assertAlmostEqual(vel_y[0], 0.0, places=5)
        self.assertAlmostEqual(vel_z[0], 0.0, places=5)


    def test_field_bc_zeros_Ey_Ez_behind_wall(self):
        _, tile = make_reflector_tile()

        wall = runko.pic.threeD.reflector_wall(walloc=5.0)
        tile.register_reflector_wall(wall)

        nonzero = lambda x, y, z: (1.0, 2.0, 3.0)
        tile.set_EBJ(nonzero, zero, zero)

        # verify E is nonzero before applying BC
        (Ex0, Ey0, Ez0), _, _ = tile.get_EBJ()
        self.assertTrue(np.any(Ey0 != 0))
        self.assertTrue(np.any(Ez0 != 0))

        tile.reflector_wall_field_bc()

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        # Wall at global x=5.0, tile mins=(0,0,0), halo_size=3.
        # Lattice index of wall: iw = int(5.0 - 0) + 3 = 8.
        # get_EBJ returns non-halo region only, so numpy index = lattice index - halo_size.
        # Cells at and behind wall: numpy x-indices 0..5 (lattice 3..8, i.e. <= iw).
        # Cells in front of wall: numpy x-indices 6..9 (lattice 9..12).
        iw_np = 6  # first numpy index in front of wall

        # Ey and Ez zeroed behind wall
        self.assertTrue(np.all(Ey[:iw_np, :, :] == 0))
        self.assertTrue(np.all(Ez[:iw_np, :, :] == 0))

        # Ey and Ez untouched in front of wall
        self.assertTrue(np.all(Ey[iw_np:, :, :] != 0))
        self.assertTrue(np.all(Ez[iw_np:, :, :] != 0))

        # Ex unchanged everywhere
        self.assertTrue(np.allclose(Ex, Ex0))


    def test_no_current_behind_wall_after_reflect_and_deposit(self):
        _, tile = make_reflector_tile()

        wall = runko.pic.threeD.reflector_wall(walloc=5.0)
        tile.register_reflector_wall(wall)

        tile.set_EBJ(zero, zero, zero)

        # inject one particle per cell near the wall (global x in [5, 6))
        P = runko.pic.threeD.ParticleState

        def particle_generator(x, y, z):
            if int(x) == 5:
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(-0.5, 0.0, 0.0))]
            return []

        tile.inject_to_each_cell(0, particle_generator)

        tile.push_particles()
        tile.reflect_particles()
        tile.deposit_current()

        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        tol = 1e-5

        # cells well behind the wall (numpy x-indices 0..3) should have ~0 current
        for ix in range(4):
            self.assertTrue(np.all(np.abs(Jx[ix, :, :]) < tol),
                            f"Jx nonzero behind wall at x-index {ix}")
            self.assertTrue(np.all(np.abs(Jy[ix, :, :]) < tol),
                            f"Jy nonzero behind wall at x-index {ix}")
            self.assertTrue(np.all(np.abs(Jz[ix, :, :]) < tol),
                            f"Jz nonzero behind wall at x-index {ix}")


    def test_overflow_particle_gets_parked(self):
        _, tile = make_reflector_tile()

        wall = runko.pic.threeD.reflector_wall(walloc=5.0)
        tile.register_reflector_wall(wall)

        tile.set_EBJ(zero, zero, zero)

        # inject particle already far behind the wall
        P = runko.pic.threeD.ParticleState
        tile.inject(0, [P(pos=(2.0, 5.5, 6.5), vel=(-0.1, 0.0, 0.0))])

        tile.push_particles()
        tile.reflect_particles()

        pos_x, pos_y, pos_z = tile.get_positions(0)
        vel_x, vel_y, vel_z = tile.get_velocities(0)

        self.assertEqual(1, len(pos_x))

        # particle should be parked near domain start
        self.assertAlmostEqual(pos_x[0], 1.0, places=4)

        # velocity zeroed
        self.assertAlmostEqual(vel_x[0], 0.0, places=5)


if __name__ == "__main__":
    unittest.main()
