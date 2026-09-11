# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import itertools
import runko_cpp_bindings.actions as actions
from runko_cpp_bindings.emf.threeD import antenna_mode
from .configuration import Configuration


class ProxyTile:
    def __init__(self, tile_id, sim_context, runtime=None):
        self.tile_id = tile_id
        self.sim_context = sim_context
        self._runtime = runtime
        self._ensure_yee = (actions.ensure_constructed_yee_lattices, actions.current_context)


    def get_EBJ(self):
        self.sim_context.eval(self._ensure_yee)
        return self.sim_context.get_EBJ(self.tile_id)


    def get_EBJ_with_halo(self):
        self.sim_context.eval(self._ensure_yee)
        return self.sim_context.get_EBJ_with_halo(self.tile_id)


    def set_EBJ(self, E, B, J):
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_yee),
                               (actions.quote, (actions.set_EBJ, actions.current_context, E, B, J))))


    def batch_set_EBJ(self, Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz):
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_yee),
                               (actions.quote, (actions.batch_set_EBJ, actions.current_context, Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz))))


    def add_current(self):
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_yee),
                               (actions.quote, (actions.add_current, actions.current_context))))


    def global_coordinate_map(self):
        return self.sim_context.global_coordinate_map(self.tile_id)


    @property
    def index(self) -> tuple[int, int, int]:
        return self.sim_context.index(self.tile_id)


    def register_antenna(self, mode: antenna_mode):
        return self.sim_context.eval((actions.register_antenna, actions.current_context, mode))


    def deposit_antenna_current(self):
        return self.sim_context.eval((actions.deposit_antenna_current, actions.current_context))




def make_independent_tile(idx, conf) -> ProxyTile:
    """
    Constructs a simulation context with a single tile.
    Rerturns a ProxyTile tile corresponding to the tile.
    """

    idx_arr = np.array(idx)
    check0 = np.zeros(len(idx)) <= idx_arr
    check1 = np.array(idx) < conf.n_tiles
    if not np.all(np.logical_and(check0, check1)):
        msg = "Given tile index is not in range defined by n_tiles configuration parameter."
        raise Exception(msg)

    runtime = actions.RuntimeInstance()
    sim = actions.SimulationContext(conf)
    sim.add_tile(idx)

    ids = sim.get_local_tile_ids()
    if len(ids) != 1:
        # Sanity check.
        raise Exception("runko internal logic error")

    return ProxyTile(ids[0], sim, runtime)
