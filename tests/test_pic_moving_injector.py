import unittest
import numpy as np
import runko


def make_config():
    config = runko.Configuration(None)
    config.Nx = 4
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 10
    config.NyMesh = 4
    config.NzMesh = 4
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 0.45
    config.field_propagator = "FDTD2"
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"
    return config


def identity_pgen(x, y, z):
    """Places one particle at each cell corner (no random offset)."""
    return runko.pic.threeD.ParticleStateBatch(pos=(x, y, z), vel=(x, y, z))


class TestBatchInjectInStripe(unittest.TestCase):

    def test_stripe_no_overlap(self):
        """Tile at (0,0,0) spans x=[0,10). Stripe [20,25) should inject nothing."""
        config = make_config()
        tile = runko.pic.threeD.Tile((0, 0, 0), config)
        tile.batch_inject_in_x_stripe(0, identity_pgen, 20.0, 25.0)

        posx, _, _ = tile.get_positions(0)
        self.assertEqual(len(posx), 0)

    def test_stripe_no_overlap_left(self):
        """Tile at (2,0,0) spans x=[20,30). Stripe [5,10) should inject nothing."""
        config = make_config()
        tile = runko.pic.threeD.Tile((2, 0, 0), config)
        tile.batch_inject_in_x_stripe(0, identity_pgen, 5.0, 10.0)

        posx, _, _ = tile.get_positions(0)
        self.assertEqual(len(posx), 0)

    def test_stripe_full_overlap(self):
        """Stripe covering entire tile should inject into all cells."""
        config = make_config()
        tile = runko.pic.threeD.Tile((0, 0, 0), config)
        tile.batch_inject_in_x_stripe(0, identity_pgen, -1.0, 11.0)

        posx, _, _ = tile.get_positions(0)
        expected = config.NxMesh * config.NyMesh * config.NzMesh
        self.assertEqual(len(posx), expected)

    def test_stripe_partial_overlap(self):
        """Stripe [3, 7) on tile (0,0,0) should inject into cells 3,4,5,6."""
        config = make_config()
        tile = runko.pic.threeD.Tile((0, 0, 0), config)
        tile.batch_inject_in_x_stripe(0, identity_pgen, 3.0, 7.0)

        posx, _, _ = tile.get_positions(0)
        n_x_cells = 4  # cells 3, 4, 5, 6
        expected = n_x_cells * config.NyMesh * config.NzMesh
        self.assertEqual(len(posx), expected)

        # all injected x-positions should be within [3, 7)
        for px in posx:
            self.assertGreaterEqual(px, 3.0)
            self.assertLess(px, 7.0)

    def test_pgen_receives_only_stripe_coords(self):
        """Verify pgen only gets x-coordinates within the stripe."""
        config = make_config()
        tile = runko.pic.threeD.Tile((0, 0, 0), config)

        received_x = []

        def capturing_pgen(x, y, z):
            received_x.append(np.array(x))
            return runko.pic.threeD.ParticleStateBatch(
                pos=(x, y, z), vel=(x, y, z))

        tile.batch_inject_in_x_stripe(0, capturing_pgen, 5.0, 8.0)

        self.assertEqual(len(received_x), 1)
        xs = received_x[0]
        for xval in xs:
            self.assertGreaterEqual(xval, 5.0)
            self.assertLess(xval, 8.0)

        # should be cells 5, 6, 7 => 3 x-cells
        unique_x = sorted(set(int(round(v)) for v in xs))
        self.assertEqual(unique_x, [5, 6, 7])

    def test_multiple_species(self):
        """Inject species 0 and 1 in same stripe independently."""
        config = make_config()
        tile = runko.pic.threeD.Tile((0, 0, 0), config)

        tile.batch_inject_in_x_stripe(0, identity_pgen, 2.0, 5.0)
        tile.batch_inject_in_x_stripe(1, identity_pgen, 2.0, 5.0)

        posx0, _, _ = tile.get_positions(0)
        posx1, _, _ = tile.get_positions(1)
        expected = 3 * config.NyMesh * config.NzMesh
        self.assertEqual(len(posx0), expected)
        self.assertEqual(len(posx1), expected)

    def test_stripe_exact_tile_boundary(self):
        """Stripe matching tile bounds exactly should inject all cells."""
        config = make_config()
        tile = runko.pic.threeD.Tile((1, 0, 0), config)
        # tile (1,0,0) spans x=[10,20)
        tile.batch_inject_in_x_stripe(0, identity_pgen, 10.0, 20.0)

        posx, _, _ = tile.get_positions(0)
        expected = config.NxMesh * config.NyMesh * config.NzMesh
        self.assertEqual(len(posx), expected)

    def test_batch_inject_to_cells_unchanged(self):
        """Refactored batch_inject_to_cells still produces same results."""
        config = make_config()

        tile_a = runko.pic.threeD.Tile((0, 0, 0), config)
        tile_b = runko.pic.threeD.Tile((0, 0, 0), config)

        tile_a.batch_inject_to_cells(0, identity_pgen)
        tile_b.batch_inject_in_x_stripe(0, identity_pgen, 0.0, 10.0)

        posx_a, posy_a, posz_a = tile_a.get_positions(0)
        posx_b, posy_b, posz_b = tile_b.get_positions(0)

        np.testing.assert_array_equal(posx_a, posx_b)
        np.testing.assert_array_equal(posy_a, posy_b)
        np.testing.assert_array_equal(posz_a, posz_b)

    def test_gap_free_injection_sequence(self):
        """Two injection rounds with flow correction should leave no x-gap."""
        config = make_config()
        # use a wide tile to fit everything
        config.NxMesh = 40
        tile = runko.pic.threeD.Tile((0, 0, 0), config)

        beta_inj = 0.5
        beta_flow = 0.3
        cfl = config.cfl
        n_inj = 2

        # first injection: stripe [5, 5 + beta_inj * n_inj * cfl]
        injloc = 5.0
        stride = n_inj * cfl
        x_left_1  = injloc
        x_right_1 = injloc + beta_inj * stride

        tile.batch_inject_in_x_stripe(0, identity_pgen, x_left_1, x_right_1)
        injloc = x_right_1

        # second injection: account for flow drift
        x_left_2  = injloc - beta_flow * stride
        x_right_2 = injloc + beta_inj * stride

        tile.batch_inject_in_x_stripe(0, identity_pgen, x_left_2, x_right_2)

        # the second injection's left edge should overlap with
        # the first injection's right edge (accounting for drift)
        # x_left_2 = x_right_1 - beta_flow * stride < x_right_1
        self.assertLess(x_left_2, x_right_1,
                        "Second injection should overlap with first to fill the flow gap")

        # check that combined particles cover [x_left_1, x_right_2) without gaps
        posx, _, _ = tile.get_positions(0)
        cell_xs = sorted(set(int(round(px)) for px in posx))

        expected_min = int(round(x_left_1))
        expected_max = int(round(x_right_2)) - 1
        expected_cells = list(range(expected_min, expected_max + 1))
        self.assertEqual(cell_xs, expected_cells,
                         "No x-cell gaps in combined injection")


class MockSim:
    """Minimal simulation mock for MovingInjector tests."""
    def __init__(self, lap=1):
        self.lap = lap
    def local_tiles(self):
        return []


class TestMovingInjectorClass(unittest.TestCase):

    def test_default_state(self):
        inj = runko.MovingInjector(
            injloc=10.0, beta_inj=0.5, beta_flow=0.3, cfl=0.45)
        self.assertEqual(inj.injloc, 10.0)
        self.assertEqual(inj.beta_inj, 0.5)
        self.assertEqual(inj.beta_flow, 0.3)
        self.assertEqual(inj.cfl, 0.45)
        self.assertEqual(inj.n_inj, 1)
        self.assertTrue(inj.moving)

    def test_moving_flag_stops_injection(self):
        inj = runko.MovingInjector(
            injloc=10.0, beta_inj=0.5, beta_flow=0.3, cfl=0.45)
        inj.moving = False
        original_loc = inj.injloc

        inj.inject(MockSim(), [(0, identity_pgen)], 1)
        self.assertEqual(inj.injloc, original_loc)

    def test_skip_lap_zero(self):
        inj = runko.MovingInjector(
            injloc=10.0, beta_inj=0.5, beta_flow=0.3, cfl=0.45, n_inj=1)
        original_loc = inj.injloc

        inj.inject(MockSim(lap=0), [(0, identity_pgen)], 1)
        self.assertEqual(inj.injloc, original_loc)

    def test_stops_at_domain_boundary(self):
        inj = runko.MovingInjector(
            injloc=35.0, beta_inj=1.0, beta_flow=0.0, cfl=1.0,
            n_inj=1, Lx=40.0, Lx_margin=10.0)

        inj.inject(MockSim(), [(0, identity_pgen)], 1)
        self.assertFalse(inj.moving)
        # injloc should NOT have been updated
        self.assertEqual(inj.injloc, 35.0)


if __name__ == "__main__":
    unittest.main()
