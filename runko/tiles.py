# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import itertools

class InitTile:
    """
    Represents data that is used to initialize some tile.
    """

    def __init__(self, idx, conf):
        self._extents_wout_halo = np.array(conf.n_cells_per_tile, dtype=int)
        self._idx = idx

        def non_halo_index_space_():
            i = np.arange(self._extents_wout_halo[0], dtype=int)
            j = np.arange(self._extents_wout_halo[1], dtype=int)
            k = np.arange(self._extents_wout_halo[2], dtype=int)
            return itertools.product(i, j, k)

        self.non_halo_index_space = non_halo_index_space_

        self._Ex = np.zeros(self._extents_wout_halo)
        self._Ey = np.zeros(self._extents_wout_halo)
        self._Ez = np.zeros(self._extents_wout_halo)
        self._Bx = np.zeros(self._extents_wout_halo)
        self._By = np.zeros(self._extents_wout_halo)
        self._Bz = np.zeros(self._extents_wout_halo)
        self._Jx = np.zeros(self._extents_wout_halo)
        self._Jy = np.zeros(self._extents_wout_halo)
        self._Jz = np.zeros(self._extents_wout_halo)


    def set_EBJ(self, E, B, J):
        for idx in self.non_halo_index_space():
            self._Ex[*idx], self._Ey[*idx], self._Ez[*idx] = E(*idx)
            self._Bx[*idx], self._By[*idx], self._Bz[*idx] = B(*idx)
            self._Jx[*idx], self._Jy[*idx], self._Jz[*idx] = J(*idx)


    def get_EBJ_with_halo(self):
        return ((self._Ex, self._Ey, self._Ez),
                (self._Bx, self._By, self._Bz),
                (self._Jx, self._Jy, self._Jz))


class ProxyTile:
    def __init__(self, tile_id, sim_context):
        self.tile_id = tile_id
        self.sim_context = sim_context

    def get_EBJ_with_halo(self):
        return self.sim_context.get_EBJ_with_halo(self.tile_id)
