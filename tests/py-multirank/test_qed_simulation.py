# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""The QED tile on a 2x2x2 grid over several ranks: reactions at the top of the lap, photons
and their products pushed across tile boundaries and exchanged, the layout kept as tail
arrivals. Checks conservation across the exchange, unique ids, and every particle in its tile."""

import mpi_unittest
import tempfile
import numpy as np

import runko
from tests.qed_test_helpers import WL, make_conf, rows, iso_dirs, mj_momenta, spread

NC = 8


def broadcast_outdir():
    runtime = runko.default_runtime()
    outdir = tempfile.mkdtemp(prefix="runko-qed-mpi-") if runtime.mpi_rank() == 0 else None
    return runtime.mpi_bcast(outdir, root=0)


def make_grid(seed=1):
    # the closed-box species and sampler settings; no splitting: photon count and weight are conserved exactly
    config = make_conf(seed, 0.45 / NC, nc=NC, nsort=5, slack=0.2, channels=["compton-", "compton+"])
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.n_laps = 20
    config.n_tiles = [2, 2, 2]
    config.io_outdir = broadcast_outdir()
    return config, runko.TileGrid(config)


def totals(simulation):
    """summed over every tile of every rank: weight and count per species"""
    loc = np.zeros(6)
    for tile in simulation.local_tiles():
        for i, s in enumerate(("ph", "e-", "e+")):
            tt = tile.totals(s); loc[i] += tt["W"]; loc[3 + i] += tt["N"]
    runtime = runko.default_runtime()
    return np.sum(runtime.mpi_bcast(runtime.mpi_gather(loc)), axis=0)


def test_qed_across_tiles():
    conf, tile_grid = make_grid()
    rank = runko.default_runtime().mpi_rank()
    rng = np.random.default_rng(100 + rank)
    for idx in tile_grid.local_tile_indices():
        tile = runko.qed.threeD.Tile(idx, conf)
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)
        nph, npm = 2000, 400
        x = rng.gamma(3.0, 0.05, nph)
        tile.set_particles("ph", rows(x[:, None] * iso_dirs(nph, rng)), list(spread(nph, rng)))
        for s in ("e-", "e+"):
            tile.set_particles(s, rows(mj_momenta(npm, 2.0, rng)), [WL] * npm)
        tile_grid.add_tile(tile, idx)
    simulation = tile_grid.configure_simulation(conf)
    t0 = totals(simulation)

    def lap(x):
        x.prtcl_interact()
        x.prtcl_push()
        x.prtcl_pack_outgoing()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

    simulation.for_each_lap(lap)
    t1 = totals(simulation)

    # weights and counts are conserved across all the exchanges (Compton moves no weight
    # between species and, without splitting, makes no new photon LPs)
    asserts = []
    for i in range(6):
        asserts.append(mpi_unittest.assertEqualDeferred(True, bool(abs(t1[i] - t0[i]) <= 1e-9 * max(t0[i], 1.0))))
    # events happened, on every rank
    n_ev = sum(sum(tile.counts().values()) for tile in simulation.local_tiles())
    asserts.append(mpi_unittest.assertEqualDeferred(True, n_ev > 0))
    # every particle inside its own tile, with the layout intact
    for tile in simulation.local_tiles():
        mins, maxs = np.array(tile.mins), np.array(tile.maxs)
        for s in (0, 1, 2):
            pos = np.stack([np.asarray(a) for a in tile.get_positions(s)], 1)
            inside = np.all((pos >= mins) & (pos < maxs)) if len(pos) else True
            asserts.append(mpi_unittest.assertEqualDeferred(True, bool(inside)))
        asserts.append(mpi_unittest.assertEqualDeferred(True, tile.laid()))
    mpi_unittest.assertDeferredResults(asserts)
    # ids unique per species over the whole grid
    for s in (0, 1, 2):
        ids = np.concatenate([np.asarray(tile.get_ids(s)) for tile in simulation.local_tiles()] or [np.zeros(0)])
        all_ids = runko.default_runtime().mpi_gather(ids, root=0)
        if rank == 0:
            a = np.concatenate(all_ids)
            mpi_unittest.assertEqual(len(a), len(np.unique(a)))
        else:
            mpi_unittest.assertEqual(0, 0)


if __name__ == "__main__":
    test_qed_across_tiles()
