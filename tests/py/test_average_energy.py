# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import shutil
import tempfile
import unittest

import numpy as np

import runko
from tests.mpiio_test_helpers import make_config, make_pic_config, setup_grid_with_fields


class TestAverageEnergy(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def load(self, name):
        return np.loadtxt(os.path.join(self.outdir, name), ndmin=2)

    def test_field_energy_densities(self):
        """2x1x1 grid with constant E=(1,2,3) and B=(4,5,6).
        Averages are E^2/2 = 7 and B^2/2 = 38.5 and are appended once per lap."""

        config = make_config(Nx=2, NxMesh=10, NyMesh=11, NzMesh=13, outdir=self.outdir)
        tile_grid = runko.TileGrid(config)
        setup_grid_with_fields(tile_grid,
                               config,
                               lambda x, y, z: (1, 2, 3),
                               lambda x, y, z: (4, 5, 6),
                               lambda x, y, z: (0, 0, 0))

        sim = tile_grid.configure_simulation(config)

        def f(x):
            x.io_average_E_energy_density()
            x.io_average_B_energy_density()

        sim.for_one_lap(f)
        sim.for_one_lap(f)

        E = self.load("average_E_energy_density.txt")
        B = self.load("average_B_energy_density.txt")

        self.assertEqual(E.shape, (2, 2))
        self.assertEqual(B.shape, (2, 2))
        self.assertEqual(E[1, 0], E[0, 0] + 1)
        self.assertEqual(B[1, 0], B[0, 0] + 1)

        self.assertTrue(np.allclose(E[:, 1], 7.0, rtol=1e-5))
        self.assertTrue(np.allclose(B[:, 1], 38.5, rtol=1e-5))


    def test_kinetic_energy(self):
        """2x1x1 grid with two species: species 0 at rest and
        species 1 with u=(0.3, 0.4, 0), so that <E_k> = sqrt(1 + 0.25) - 1.
        Species 2 is configured but has no particles, so its average is zero."""

        config = make_pic_config(Nx=2, outdir=self.outdir, nspecies=3)
        tile_grid = runko.TileGrid(config)

        P = runko.pic.threeD.ParticleState

        def at_rest(x, y, z):
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(0, 0, 0))]

        def moving(x, y, z):
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(0.3, 0.4, 0))]

        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            tile.inject_to_each_cell(0, at_rest)
            tile.inject_to_each_cell(1, moving)
            tile_grid.add_tile(tile, idx)

        sim = tile_grid.configure_simulation(config)
        sim.for_one_lap(lambda x: x.io_average_kinetic_energy())

        d = self.load("average_kinetic_energy.txt")

        self.assertEqual(d.shape, (1, 4))
        self.assertEqual(d[0, 1], 0.0)
        self.assertTrue(np.allclose(d[0, 2], np.sqrt(1.25) - 1, rtol=1e-5))
        self.assertEqual(d[0, 3], 0.0)


if __name__ == "__main__":
    _ = runko.actions.RuntimeInstance()
    unittest.main()
