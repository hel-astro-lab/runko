# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import itertools
import runko_cpp_bindings.actions as actions
import runko_cpp_bindings.qed as qed
from .runko_logging import runko_logger, on_main_rank
from runko_cpp_bindings.emf.threeD import antenna_mode
from runko_cpp_bindings.pic.threeD import ParticleState, ParticleStateBatch
from .configuration import Configuration


class ProxyTile:
    def __init__(self, tile_id, sim_context, runtime=None):
        self.tile_id = tile_id
        self.sim_context = sim_context
        self._runtime = runtime
        self._logger = runko_logger("ProxyTile")
        self._ensure_yee = (actions.ensure_constructed_yee_lattices, actions.current_context)
        self._ensure_prtcls = (actions.ensure_constructed_particle_containers, actions.current_context)


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


    @property
    def mins(self) -> tuple[int, int, int]:
        return self.sim_context.mins(self.tile_id)


    @property
    def maxs(self) -> tuple[int, int, int]:
        return self.sim_context.maxs(self.tile_id)


    def register_antenna(self, mode: antenna_mode):
        return self.sim_context.eval((actions.register_antenna, actions.current_context, mode))


    def deposit_antenna_current(self):
        return self.sim_context.eval((actions.deposit_antenna_current, actions.current_context))


    def push_e(self):
        return self.sim_context.eval((actions.push_e, actions.current_context))


    def push_half_b(self):
        return self.sim_context.eval((actions.push_half_b, actions.current_context))


    def filter_current(self):
        return self.sim_context.eval((actions.filter_current, actions.current_context))


    def register_edge_bc(self, edge_bc):
        return self.sim_context.eval((actions.register_edge_bc, actions.current_context, edge_bc))


    def apply_edge_bc(self, edge_bc, mode):
        return self.sim_context.eval((actions.apply_edge_bc, actions.current_context, edge_bc, mode))


    def apply_edge_bcs(self, mode):
        return self.sim_context.eval((actions.apply_edge_bcs, actions.current_context, mode))


    def get_positions(self, ptype: int):
        self.sim_context.eval(self._ensure_prtcls)
        return self.sim_context.get_positions(self.tile_id, ptype)


    def get_velocities(self, ptype: int):
        self.sim_context.eval(self._ensure_prtcls)
        return self.sim_context.get_velocities(self.tile_id, ptype)


    def get_ids(self, ptype: int):
        self.sim_context.eval(self._ensure_prtcls)
        return self.sim_context.get_ids(self.tile_id, ptype)


    def get_payload(self, ptype: int, c: int = 0):
        """payload component c of the alive particles of a species"""
        self.sim_context.eval(self._ensure_prtcls)
        return self.sim_context.get_payload(self.tile_id, ptype, c)


    def inject_to_each_cell(self, ptype: int, pgen):
        """pgen(x, y, z) returns the states of a cell, or (states, payload) for a species
        with a payload"""
        prg = (actions.inject_to_each_cell, actions.current_context, ptype, pgen)
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_prtcls),
                               (actions.quote, prg)))


    def inject(self, ptype: int, particles: list[ParticleState], payload=None):
        """payload: shape (n_payload, N), or (N,) when n_payload is 1; required iff the
        species carries one"""
        prg = (actions.inject, actions.current_context, ptype, particles)
        if payload is not None:
            import numpy as np
            prg += (np.asarray(payload, dtype=float).tolist(),)
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_prtcls),
                               (actions.quote, prg)))


    def batch_inject_to_cells(self, ptype: int, batch_pgen):
        prg = (actions.batch_inject_to_cells, actions.current_context, ptype, batch_pgen)
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_prtcls),
                               (actions.quote, prg)))


    def batch_inject_in_x_stripe(self, ptype: int, batch_pgen, x_left: float, x_right: float):
        prg = (actions.batch_inject_in_x_stripe, actions.current_context, ptype, batch_pgen, x_left, x_right)
        self.sim_context.eval((actions.sequence,
                               (actions.quote, self._ensure_prtcls),
                               (actions.quote, prg)))


    def push_particles(self):
        self.sim_context.eval((actions.push_particles, actions.current_context))


    def deposit_current(self):
        self.sim_context.eval((actions.deposit_current, actions.current_context))


    def sort_particles(self):
        self.sim_context.eval((actions.sort_particles, actions.current_context))


    def register_reflector_wall(self, wall):
        self.sim_context.eval((actions.register_reflector_wall, actions.current_context, wall))


    def reflect_particles(self):
        self.sim_context.eval((actions.reflect_particles, actions.current_context))


    def advance_reflector_walls(self):
        self.sim_context.eval((actions.advance_reflector_walls, actions.current_context))


    def pack_outgoing_particles(self):
        if on_main_rank():
            self._logger.warn("pack_outgoing_particles is depricated and implemented as no-op.")


    # ---- the QED Monte Carlo of this tile (a tile of runko.qed.threeD.Tile); species are
    # named "ph", "e-", "e+"
    def species(self, s: str) -> int:
        """the container index of a QED species"""
        return qed.species(self.sim_context, self.tile_id, s)

    def box(self):
        return qed.box(self.sim_context, self.tile_id)

    def step(self):
        return qed.step(self.sim_context, self.tile_id)

    def laid(self):
        return qed.laid(self.sim_context, self.tile_id)

    def layout(self):
        qed.layout(self.sim_context, self.tile_id)

    def interact(self):
        qed.interact(self.sim_context, self.tile_id)

    def get_weights(self, ptype: int):
        """the weights of the live particles of container ptype: payload component 0 of a
        species with a payload (the photons), the species weight w{i} of one without"""
        return qed.get_weights(self.sim_context, self.tile_id, ptype)

    def set_particles(self, s: str, u, w):
        qed.set_particles(self.sim_context, self.tile_id, s, u, w)

    def reset_momenta(self, s: str, u):
        qed.reset_momenta(self.sim_context, self.tile_id, s, u)

    def get_slot_keys(self, s: str):
        return qed.get_slot_keys(self.sim_context, self.tile_id, s)

    def get_tables(self, s: str):
        return qed.get_tables(self.sim_context, self.tile_id, s)

    def totals(self, s: str):
        return qed.totals(self.sim_context, self.tile_id, s)

    def tau_T(self):
        return qed.tau_T(self.sim_context, self.tile_id)

    def tau_ssa(self, x: float):
        return qed.tau_ssa(self.sim_context, self.tile_id, x)

    def pair_sum(self, channel: str, weight: bool = False):
        return qed.pair_sum(self.sim_context, self.tile_id, channel, weight)

    def counts(self):
        return qed.counts(self.sim_context, self.tile_id)

    def moved(self):
        return qed.moved(self.sim_context, self.tile_id)

    def report(self):
        return qed.report(self.sim_context, self.tile_id)

    def inject_pairs_plaw(self, ell: float, p: float, gmin: float, gmax: float):
        qed.inject_pairs_plaw(self.sim_context, self.tile_id, ell, p, gmin, gmax)

    def inject_photons(self, ell: float, theta: float, w: float, planck: bool = True):
        qed.inject_photons(self.sim_context, self.tile_id, ell, theta, w, planck)

    def heat_diffusive(self, ell: float):
        qed.heat_diffusive(self.sim_context, self.tile_id, ell)

    def escape_photons(self, tau_T: float = -1.0, dt_esc: float = -1.0):
        qed.escape_photons(self.sim_context, self.tile_id, tau_T, dt_esc)

    def budget(self):
        return qed.budget(self.sim_context, self.tile_id)

    def reset_budget(self):
        qed.reset_budget(self.sim_context, self.tile_id)

    def escaped_spectrum(self):
        return qed.escaped_spectrum(self.sim_context, self.tile_id)

    def escaped_mu_spectrum(self):
        return qed.escaped_mu_spectrum(self.sim_context, self.tile_id)




def make_independent_tile(idx, conf) -> ProxyTile:
    """
    Constructs a simulation context with a single tile.
    Rerturns a ProxyTile tile corresponding to the tile.
    """

    import numpy as np
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


def make_independent_qed_tile(idx, conf) -> ProxyTile:
    """
    The tile of make_independent_tile whose species react through the QED Monte Carlo:
    its lattice, particle containers and QED state built from conf.
    """

    tile = make_independent_tile(idx, conf)
    tile.sim_context.eval((actions.ensure_constructed_qed, actions.current_context))
    return tile
