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


    def test_particle_in_front_of_wall_unchanged(self):
        _, tile = make_reflector_tile()

        wall = runko.pic.threeD.reflector_wall(walloc=5.0)
        tile.register_reflector_wall(wall)

        tile.set_EBJ(zero, zero, zero)

        # particle well in front of the wall, moving away
        P = runko.pic.threeD.ParticleState
        tile.inject(0, [P(pos=(7.0, 5.5, 6.5), vel=(0.5, 0.2, -0.1))])

        tile.push_particles()

        # record post-push state
        px_pre, py_pre, pz_pre = tile.get_positions(0)
        vx_pre, vy_pre, vz_pre = tile.get_velocities(0)

        tile.reflect_particles()

        px, py, pz = tile.get_positions(0)
        vx, vy, vz = tile.get_velocities(0)

        # position and velocity must be untouched by reflect
        self.assertAlmostEqual(px[0], px_pre[0], places=5)
        self.assertAlmostEqual(py[0], py_pre[0], places=5)
        self.assertAlmostEqual(pz[0], pz_pre[0], places=5)
        self.assertAlmostEqual(vx[0], vx_pre[0], places=5)
        self.assertAlmostEqual(vy[0], vy_pre[0], places=5)
        self.assertAlmostEqual(vz[0], vz_pre[0], places=5)


    def test_moving_wall_reflects_ux(self):
        _, tile = make_reflector_tile()

        betawall = 0.5
        gammawall = 1.0 / np.sqrt(1.0 - betawall**2)
        wall = runko.pic.threeD.reflector_wall(
            walloc=5.0, betawall=betawall, gammawall=gammawall)
        tile.register_reflector_wall(wall)

        tile.set_EBJ(zero, zero, zero)

        # particle near the wall, moving toward it
        P = runko.pic.threeD.ParticleState
        ux0 = -0.8
        tile.inject(0, [P(pos=(5.5, 5.5, 6.5), vel=(ux0, 0.0, 0.0))])

        tile.push_particles()
        tile.reflect_particles()

        vx, vy, vz = tile.get_velocities(0)
        px, _, _ = tile.get_positions(0)

        # ux_new = gammawall^2 * gam * (2*beta - ux/gam * (1 + beta^2))
        gam = np.sqrt(1.0 + ux0**2)
        ux_expected = gammawall**2 * gam * (
            2.0 * betawall - ux0 / gam * (1.0 + betawall**2))

        self.assertAlmostEqual(vx[0], ux_expected, places=3)

        # particle must be in front of the wall
        self.assertGreater(px[0], 5.0)

        # transverse velocity unchanged
        self.assertAlmostEqual(vy[0], 0.0, places=5)
        self.assertAlmostEqual(vz[0], 0.0, places=5)


    def test_mixed_particles_reflect_park_skip(self):
        _, tile = make_reflector_tile()

        wall = runko.pic.threeD.reflector_wall(walloc=5.0)
        tile.register_reflector_wall(wall)

        tile.set_EBJ(zero, zero, zero)

        P = runko.pic.threeD.ParticleState
        tile.inject(0, [
            # particle 0: well in front of wall (skip)
            P(pos=(7.0, 5.5, 6.5), vel=(0.5, 0.0, 0.0)),
            # particle 1: close behind wall, genuine crossing (reflect)
            P(pos=(5.5, 5.5, 6.5), vel=(-1.0, 0.0, 0.0)),
            # particle 2: far behind wall, artifact (park)
            P(pos=(2.0, 5.5, 6.5), vel=(-0.1, 0.0, 0.0)),
        ])

        tile.push_particles()
        tile.reflect_particles()

        px, _, _ = tile.get_positions(0)
        vx, _, _ = tile.get_velocities(0)

        self.assertEqual(3, len(px))

        # particle 0 (skip): velocity unchanged, in front of wall
        self.assertAlmostEqual(vx[0], 0.5, places=4)
        self.assertGreater(px[0], 5.0)

        # particle 1 (reflect): ux flipped, back in front of wall
        self.assertGreater(vx[1], 0)
        self.assertAlmostEqual(abs(vx[1]), 1.0, places=4)
        self.assertGreater(px[1], 5.0)

        # particle 2 (park): parked at x=1.0, ux=0
        self.assertAlmostEqual(px[2], 1.0, places=4)
        self.assertAlmostEqual(vx[2], 0.0, places=5)


if __name__ == "__main__":
    unittest.main()
