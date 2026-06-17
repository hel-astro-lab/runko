# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

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


if __name__ == "__main__":
    unittest.main()
