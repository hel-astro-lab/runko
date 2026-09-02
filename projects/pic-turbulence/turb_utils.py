# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared helpers for the pic-turbulence analysis scripts."""

import glob
import re
from math import atan2, degrees
import numpy as np
import matplotlib.pyplot as plt
from runko.oscillating_langevin_antenna import sample_oscillating_langevin_antenna

# snapshot arrays come out of the readers as (nz, ny, nx)
Z, Y, X = 0, 1, 2


def laps_per_eddy(conf):
    """Laps in one light crossing time of the driving-scale eddy L / modes_perp."""
    L = conf.n_tiles[0] * conf.n_cells_per_tile[0]
    return (L / conf.modes_perp) / conf.cfl


def numerical_mass(conf, species):
    """Particle mass in code units, m_s |q_s|, as the simulation sets it up."""
    return getattr(conf, f"m{species}") * abs(getattr(conf, f"q{species}"))


def guide_field(conf):
    """Uniform B_z the run was initialized with, from the cold species-0 magnetization."""
    return np.sqrt(conf.ppc * conf.cfl**2 * conf.sigma * numerical_mass(conf, 0))


def d_forward(f, axis):
    """f(n+1) - f(n) across a periodic box."""
    return np.roll(f, -1, axis=axis) - f


def d_backward(f, axis):
    """f(n) - f(n-1) across a periodic box."""
    return f - np.roll(f, 1, axis=axis)


def curl_forward(vx, vy, vz):
    """Yee curl on forward differences, the one that advances B from E."""
    return (d_forward(vz, Y) - d_forward(vy, Z),
            d_forward(vx, Z) - d_forward(vz, X),
            d_forward(vy, X) - d_forward(vx, Y))


def curl_backward(vx, vy, vz):
    """Yee curl on backward differences, the one that advances E from B."""
    return (d_backward(vz, Y) - d_backward(vy, Z),
            d_backward(vx, Z) - d_backward(vz, X),
            d_backward(vy, X) - d_backward(vx, Y))


def div_forward(vx, vy, vz):
    """Divergence conjugate to curl_forward; annihilates whatever advances B."""
    return d_forward(vx, X) + d_forward(vy, Y) + d_forward(vz, Z)


def div_backward(vx, vy, vz):
    """Divergence conjugate to curl_backward; this is the Gauss operator for E."""
    return d_backward(vx, X) + d_backward(vy, Y) + d_backward(vz, Z)


def langevin_antenna(conf):
    """Rebuild the driving antenna, replaying exactly what pic.py hands to runko.

    The random stream is a function of the seed alone, so drawing the phases and the
    Langevin series in pic.py's order reproduces the realization that was run.
    Returns the vector potential amplitude, one wave vector per mode, and the complex
    per-lap coefficients indexed [mode, lap].
    """
    rng = np.random.default_rng(seed=42)

    L = conf.n_tiles[0] * conf.n_cells_per_tile[0]
    kpar = 2.0 * np.pi * conf.modes_par / L
    kperp = 2.0 * np.pi * conf.modes_perp / L
    modes = [(conf.modes_perp, 0, conf.modes_par),
             (conf.modes_perp, 0, -conf.modes_par),
             (0, conf.modes_perp, conf.modes_par),
             (0, conf.modes_perp, -conf.modes_par)]

    A = np.array([0.0, 0.0, kpar * guide_field(conf) * conf.Bperp_to_B0
                  / (kperp**2 * np.sqrt(len(modes)))])

    # Alfven speed from the enthalpy-corrected magnetization, as pic.py derives it
    d0 = conf.theta0
    d1 = conf.theta1_to_theta0 * d0
    enth = lambda d: (1 + 3.89 * d + 5.56 * d**2) / (1 + 1.39 * d)
    sigma = guide_field(conf)**2 / (conf.ppc * conf.cfl**2 * (
        enth(d0) * numerical_mass(conf, 0) + enth(d1) * numerical_mass(conf, 1)))
    valf = np.sqrt(sigma / (1.0 + sigma))
    linear_freq = kpar * conf.cfl * valf

    coeffs = []
    for _ in modes:
        phase = np.exp(2j * np.pi * rng.random())
        coeffs.append(phase * sample_oscillating_langevin_antenna(
            size=conf.n_laps,
            characteristic_freq=conf.antenna_freq * linear_freq,
            decorrelation_rate=-conf.antenna_decorr * linear_freq,
            gen=rng))

    return A, 2.0 * np.pi * np.array(modes, dtype=float) / L, np.array(coeffs)


def antenna_current(antenna, lap, shape, cfl):
    """Antenna current on the Yee grid at one lap, as runko deposits it.

    The vector potential is laid down on the same staggered points as E, then
    J = cfl curl_backward(curl_forward(A)) exactly as Tile::deposit_antenna_current
    does it, so this is the discrete current the run actually saw, not its continuum
    limit.
    """
    A, kvecs, coeffs = antenna
    nz, ny, nx = shape
    z = np.arange(nz)[:, None, None]
    y = np.arange(ny)[None, :, None]
    x = np.arange(nx)[None, None, :]

    pot = [np.zeros(shape) for _ in range(3)]
    for k, w in zip(kvecs, coeffs[:, lap]):
        # each component sits half a cell forward along its own axis
        for c, (dx, dy, dz) in enumerate(((0.5, 0.0, 0.0),
                                          (0.0, 0.5, 0.0),
                                          (0.0, 0.0, 0.5))):
            if A[c] == 0.0:
                continue
            phi = k[0] * (x + dx) + k[1] * (y + dy) + k[2] * (z + dz)
            pot[c] += A[c] * (w.real * np.cos(phi) - w.imag * np.sin(phi))

    return tuple(cfl * j for j in curl_backward(*curl_forward(*pot)))


def snapshot_files(outdir, pattern):
    """Paths of the snapshots matching pattern in outdir, ordered by lap."""
    paths = glob.glob(f"{outdir}/{pattern}")
    if not paths:
        raise SystemExit(f"no {pattern} snapshots found in {outdir}")
    return sorted(paths, key=lambda p: int(re.search(r"_(\d+)\.bin$", p).group(1)))


def set_plot_style():
    """Publication styling shared by every figure; see plotting conventions."""
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    plt.rc('xtick', top=True, direction='out', labelsize=7)
    plt.rc('ytick', right=True, direction='out', labelsize=7)
    plt.rc('axes', labelsize=8)
    plt.rc('legend', handlelength=4.0)


def label_lines(lines, xvals, fontsize=6, color=None, rotate=True, **kwargs):
    """Write each line's label onto the line itself, in place of a legend.

    The label sits at the given x on an opaque patch that punches a hole in the
    line. With rotate it follows the slope of the curve there as it appears on
    screen; without it every label stays horizontal.
    """
    for line, x in zip(lines, xvals):
        ax = line.axes
        xd, yd = np.asarray(line.get_xdata()), np.asarray(line.get_ydata())

        i = np.searchsorted(xd, x) - 1
        i = min(max(i, 0), len(xd) - 2)
        xa, xb, ya, yb = xd[i], xd[i + 1], yd[i], yd[i + 1]
        y = ya + (yb - ya) * (x - xa) / (xb - xa)

        # slope in display coordinates, so the text follows the drawn curve
        rotation = 0.0
        if rotate:
            (sxa, sya), (sxb, syb) = ax.transData.transform([(xa, ya), (xb, yb)])
            rotation = (degrees(atan2(syb - sya, sxb - sxa)) + 90.0) % 180.0 - 90.0

        ax.text(x, y, line.get_label(),
                color=color or line.get_color(),
                rotation=rotation, rotation_mode='anchor',
                fontsize=fontsize, ha='center', va='center', zorder=2.5,
                bbox=dict(boxstyle='square', pad=0.1,
                          fc=ax.get_facecolor(), ec='none'),
                **kwargs)


def save_figure(fig, outdir, name):
    """Write the figure into the simulation directory it was made from."""
    outfile = f"{outdir}/{name}.pdf" if outdir else f"{name}.pdf"
    fig.savefig(outfile)
    print(f"saved {outfile}")
