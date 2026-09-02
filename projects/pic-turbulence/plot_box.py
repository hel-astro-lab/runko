# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Three-dimensional rendering of the simulation box peripherals for a field variable."""

import argparse
import re
from itertools import combinations, product

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from runko.postprocessing import read_config
from runko.mpiio_reader import read_field
from turb_utils import (guide_field, snapshot_files, set_plot_style,
                        save_figure)


default_values = { 'cmap':  "viridis", 'vmin':  -1.0, 'vmax':  +1.0, 'title': '', 'derived': False, }

default_turbulence_values = {
    'n0':    {'title': r"$n_-/n_0$",  'cmap': "viridis", 'vmin': 0.0, 'vmax': 4.0},
    'n1':    {'title': r"$n_+/n_0$",  'cmap': "viridis", 'vmin': 0.0, 'vmax': 4.0},
    'ntot':  {'title': r"$n/n_0$",    'cmap': "viridis", 'vmin': 0.0, 'vmax': 4.0, 'derived': True},

    'ex':    {'title': r"$E_x/B_0$",  'cmap': "RdBu", 'vmin': -0.5, 'vmax': 0.5},
    'ey':    {'title': r"$E_y/B_0$",  'cmap': "RdBu", 'vmin': -0.5, 'vmax': 0.5},
    'ez':    {'title': r"$E_z/B_0$",  'cmap': "RdBu", 'vmin': -0.5, 'vmax': 0.5},

    'bx':    {'title': r"$B_x/B_0$",  'cmap': "RdBu", 'vmin': -1.0, 'vmax': 1.0},
    'by':    {'title': r"$B_y/B_0$",  'cmap': "RdBu", 'vmin': -1.0, 'vmax': 1.0},
    'bz':    {'title': r"$B_z/B_0$",  'cmap': "RdBu", 'vmin':  0.0, 'vmax': 2.0},
    'bperp': {'title': r"$B_\perp/B_0$", 'cmap': "magma", 'vmin': 0.0, 'vmax': 1.0, 'derived': True},

    'jx':    {'title': r"$J_x/J_0$",  'cmap': "RdBu", 'vmin': -1.0, 'vmax': 1.0},
    'jy':    {'title': r"$J_y/J_0$",  'cmap': "RdBu", 'vmin': -1.0, 'vmax': 1.0},
    'jz':    {'title': r"$J_z/J_0$",  'cmap': "RdBu", 'vmin': -1.0, 'vmax': 1.0},
    'jpar':  {'title': r"$J_\parallel/J_0$", 'cmap': "RdBu", 'vmin': -1.0, 'vmax': 1.0, 'derived': True},
    'jperp': {'title': r"$J_\perp/J_0$", 'cmap': "magma", 'vmin': 0.0, 'vmax': 1.0, 'derived': True},
    'j':     {'title': r"$|J|/J_0$",  'cmap': "magma", 'vmin': 0.0, 'vmax': 1.0, 'derived': True},

    'b2':    {'title': r"$U_B/U_{B_0}$", 'cmap': "magma", 'vmin': 0.0, 'vmax': 2.0, 'derived': True},
    'e2':    {'title': r"$U_E/U_{B_0}$", 'cmap': "magma", 'vmin': 0.0, 'vmax': 0.2, 'derived': True},
}


def preset(var):
    """Plotting settings for var: the defaults, overwritten by the per-variable entries."""
    args = dict(default_values)
    args.update(default_turbulence_values[var])
    return args


def get_normalization(var, conf):
    """Code units per plot unit. Densities in n_0, fields in B_0, currents in q_e n_0 c^2."""
    n0 = 2.0 * conf.ppc              # total number density per cell
    qe = abs(conf.q0)
    B0 = guide_field(conf)

    if var in ('n0', 'n1'):
        return conf.ppc
    if var == 'ntot':
        return n0
    if var in ('ex', 'ey', 'ez', 'bx', 'by', 'bz', 'bperp'):
        return B0
    if var in ('b2', 'e2'):
        return B0**2
    return qe * n0 * conf.cfl**2     # jx jy jz jpar jperp j


def read_variable(path, var):
    """The (nz, ny, nx) cube of var, assembled from the raw fields when derived."""
    g = lambda name: read_field(path, name)

    if not preset(var)['derived']:
        return g(var)

    if var == 'ntot':
        return g('n0') + g('n1')
    if var == 'bperp':                                  # perpendicular to the guide field
        return np.sqrt(g('bx')**2 + g('by')**2)
    if var == 'b2':
        return g('bx')**2 + g('by')**2 + g('bz')**2
    if var == 'e2':
        return g('ex')**2 + g('ey')**2 + g('ez')**2
    if var == 'j':
        return np.sqrt(g('jx')**2 + g('jy')**2 + g('jz')**2)
    if var == 'jpar':                                   # guide field lies along z
        return g('jz')
    if var == 'jperp':
        return np.sqrt(g('jx')**2 + g('jy')**2)


def draw_box(ax, data, cmap, norm):
    """Draw the three camera-facing faces of the cube on a unit box.

    data is indexed (z, y, x), so the visible faces are the far z plane and the near
    y and x planes; each is mapped onto the unit cube with its own meshgrid.
    """
    nz, ny, nx = data.shape
    lin = lambda n: np.linspace(0.0, 1.0, n)

    # top, z = 1
    X, Y = np.meshgrid(lin(nx), lin(ny))
    ax.plot_surface(X, Y, np.ones_like(X), rstride=1, cstride=1, shade=False,
                    antialiased=True, facecolors=cmap(norm(data[-1, :, :])))

    # front, y = 0
    X, Z = np.meshgrid(lin(nx), lin(nz))
    ax.plot_surface(X, np.zeros_like(X), Z, rstride=1, cstride=1, shade=False,
                    antialiased=True, facecolors=cmap(norm(data[:, 0, :])))

    # left, x = 0
    Y, Z = np.meshgrid(lin(ny), lin(nz))
    ax.plot_surface(np.zeros_like(Y), Y, Z, rstride=1, cstride=1, shade=False,
                    antialiased=True, facecolors=cmap(norm(data[:, :, 0])))

    # outline every cube edge except the three hidden behind the box
    hidden = np.array([1.0, 1.0, 0.0])
    corners = [np.array(c) for c in product([0.0, 1.0], repeat=3)]
    for p0, p1 in combinations(corners, 2):
        if np.sum(np.abs(p0 - p1)) != 1.0:
            continue
        if np.array_equal(p0, hidden) or np.array_equal(p1, hidden):
            continue
        ax.plot(*zip(p0, p1), color='k', linestyle='dashed', lw=0.4, zorder=20)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Turbulence box peripherals in 3D")
    parser.add_argument("--dir", required=True, help="simulation output directory")
    parser.add_argument("--var", required=True, choices=sorted(default_turbulence_values),
                        help="variable to plot")
    parser.add_argument("--lap", type=int, default=None,
                        help="lap to plot (default: every snapshot in the directory)")
    args = parser.parse_args()

    conf = read_config(args.dir)
    opts = preset(args.var)
    norm = matplotlib.colors.Normalize(vmin=opts['vmin'], vmax=opts['vmax'])
    cmap = matplotlib.colormaps[opts['cmap']]
    scale = get_normalization(args.var, conf)

    paths = snapshot_files(args.dir, "flds_*.bin")
    if args.lap is not None:
        paths = [p for p in paths
                 if int(re.search(r"_(\d+)\.bin$", p).group(1)) == args.lap]
        if not paths:
            raise SystemExit(f"no flds snapshot for lap {args.lap} in {args.dir}")

    set_plot_style()

    for path in paths:
        lap = int(re.search(r"_(\d+)\.bin$", path).group(1))
        data = read_variable(path, args.var) / scale

        fig = plt.figure(1, figsize=(3.45, 3.45), dpi=200)
        ax = fig.add_subplot(111, projection="3d")

        draw_box(ax, data, cmap, norm)

        ax.view_init(55.0, -120.0)
        ax.set_box_aspect((1.0, 1.0, 1.0))
        ax.set_axis_off()
        fig.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)

        cax = fig.add_axes([0.3, 0.06, 0.4, 0.01])
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
                          cax=cax, orientation="horizontal", ticklocation="bottom")
        cb.ax.tick_params(labelsize=6)
        for spine in cb.ax.spines.values():
            spine.set_linewidth(0.5)
        fig.text(0.22, 0.06, opts['title'], fontsize=8, ha='center', va='center')

        save_figure(fig, args.dir, f"plot_box_{args.var}_{lap:04d}", ext="png")
        fig.clear()
