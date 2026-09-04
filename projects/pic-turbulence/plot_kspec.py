# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Perpendicular power spectra of the field cubes, binned in k_perp across the xy plane."""

import argparse
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from runko.postprocessing import read_config
from runko.mpiio_reader import read_field
from turb_utils import (laps_per_eddy, guide_field, current_unit, snapshot_files,
                        set_plot_style, save_figure)


components = {
    'b':     ('bx', 'by', 'bz'),
    'bperp': ('bx', 'by'),
    'bpar':  ('bz',),
    'e':     ('ex', 'ey', 'ez'),
    'j':     ('jx', 'jy', 'jz'),
}

titles = {
    'b':     r"$E_B(k_\perp)$",
    'bperp': r"$E_{B_\perp}(k_\perp)$",
    'bpar':  r"$E_{B_\parallel}(k_\perp)$",
    'e':     r"$E_E(k_\perp)$",
    'j':     r"$E_J(k_\perp)$",
}


def kperp_spectrum(fields):
    """Shell-summed power of a vector field, binned in k_perp = |(kx, ky)|.

    Returns the perpendicular wavenumber in inverse cells and E(k_perp) such that
    sum(E) * dk equals half the box mean of the summed squared components. The shell
    sum already carries the ~2 pi k mode count of the annulus, so no Jacobian is
    applied on top of it; dividing by N^2 is what makes Parseval hold.
    """
    nz, ny, nx = fields[0].shape

    power = np.zeros((nz, ny, nx))
    for f in fields:
        power += np.abs(np.fft.fftn(f))**2
    power /= (nz * ny * nx)**2

    kx = np.fft.fftfreq(nx) * nx                  # cycles across the box
    ky = np.fft.fftfreq(ny) * ny
    KY, KX = np.meshgrid(ky, kx, indexing='ij')
    rperp = np.broadcast_to(np.rint(np.hypot(KX, KY)).astype(int), (nz, ny, nx))

    dk = 2.0 * np.pi / nx                         # cubic box of nx cells
    shell = np.bincount(rperp.ravel(), weights=power.ravel()) / (2.0 * dk)
    k = np.arange(len(shell)) * dk

    # beyond the Nyquist radius the annuli are only partly inside the sampled cube
    keep = min(nx, ny) // 2 + 1
    return k[:keep], shell[:keep], dk


def local_slope(k, spec, window=5):
    """d ln E / d ln k, smoothed over a few bins to damp the shell-to-shell scatter."""
    slope = np.full(len(k), np.nan)
    ok = spec > 0
    if ok.sum() < 2:                       # an empty snapshot has no slope to show
        return slope

    s = np.gradient(np.log(spec[ok]), np.log(k[ok]))
    if len(s) >= window:                   # convolve 'same' only when it preserves length
        kernel = np.ones(window) / window
        s = (np.convolve(s, kernel, mode='same')
             / np.convolve(np.ones(len(s)), kernel, mode='same'))

    slope[ok] = s
    return slope


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Turbulence perpendicular power spectra")
    parser.add_argument("--dir", required=True, help="simulation output directory")
    parser.add_argument("--var", default="b", choices=sorted(components),
                        help="field whose spectrum to plot")
    parser.add_argument("--lap", type=int, default=None,
                        help="lap to plot (default: every snapshot, overlaid)")
    args = parser.parse_args()

    conf = read_config(args.dir)
    teddy = laps_per_eddy(conf)
    unit = current_unit(conf) if args.var == 'j' else guide_field(conf)

    paths = snapshot_files(args.dir, "flds_*.bin")
    if args.lap is not None:
        paths = [p for p in paths
                 if int(re.search(r"_(\d+)\.bin$", p).group(1)) == args.lap]
        if not paths:
            raise SystemExit(f"no flds snapshot for lap {args.lap} in {args.dir}")

    times, spectra = [], []
    for path in paths:
        fields = [read_field(path, c) for c in components[args.var]]
        k, spec, dk = kperp_spectrum(fields)
        spectra.append(spec / unit**2)
        times.append(int(re.search(r"_(\d+)\.bin$", path).group(1)) / teddy)

    times = np.array(times)
    kd = k * conf.n_cells_per_skindepth          # wavenumber in units of 1/d_e
    print(f"{len(times)} snapshots, k_perp d_e from {kd[1]:.4f} to {kd[-1]:.4f}")

    set_plot_style()

    # spectrum panel at the golden ratio, slope panel half its height; margins in
    # inches so the panel shapes hold when the margins are retuned
    golden = 0.5 * (1.0 + np.sqrt(5.0))
    axleft, axright = 0.20, 0.96
    bottom_in, top_in = 0.42, 0.52
    width_in = 3.25
    spec_in = (axright - axleft) * width_in / golden
    height_in = 1.5 * spec_in + bottom_in + top_in
    axbottom = bottom_in / height_in
    axtop = 1.0 - top_in / height_in

    fig = plt.figure(1, figsize=(width_in, height_in))
    gs = plt.GridSpec(3, 1)
    gs.update(hspace=0.0)
    axE = plt.subplot(gs[0:2, 0])
    axS = plt.subplot(gs[2, 0], sharex=axE)
    for ax in (axE, axS):
        ax.minorticks_on()

    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=max(times[-1], 1e-12))
    cmap = matplotlib.colormaps['turbo']

    for t, spec in zip(times, spectra):
        col = cmap(norm(t)) if len(times) > 1 else 'k'
        axE.plot(kd[1:], np.where(spec[1:] > 0, spec[1:], np.nan),
                 color=col, linewidth=0.8)
        axS.plot(kd[1:], local_slope(kd[1:], spec[1:]), color=col, linewidth=0.8)

    axE.set_xscale('log')
    axE.set_yscale('log')
    axE.set_xlim(1.0e-2, 3.0)
    # normalizing by the field unit pins the peak near unity, so the limits can be
    # fixed; the floor cuts the filter-suppressed tail below the last decade of physics
    axE.set_ylim(1.0e-6, 1.0e1)
    # skip the decade sitting on the shared spine so it misses the slope panel's 0
    axE.set_yticks([1.0e-4, 1.0e-2, 1.0e0])
    axE.set_ylabel(titles[args.var])
    axE.tick_params(labelbottom=False)

    axS.set_ylim(-4.0, 0.0)
    axS.set_ylabel(r"$\dfrac{\mathrm{d}\ln E}{\mathrm{d}\ln k_\perp}$")
    axS.set_xlabel(r"$k_\perp d_e$")

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    if len(times) > 1:
        cax = fig.add_axes([axleft, axtop + 0.02, axright - axleft,
                            (axtop - axbottom) * 0.03])
        cb = matplotlib.colorbar.ColorbarBase(
            cax, cmap=cmap, norm=norm, orientation='horizontal', ticklocation='top')
        cb.set_label(r"$t / t_0$")
        save_figure(fig, args.dir, f"plot_kspec_{args.var}")
    else:
        save_figure(fig, args.dir, f"plot_kspec_{args.var}_{args.lap:04d}")

    plt.show()
