"""k-space Fourier spectrum of upstream EM fluctuations at one snapshot.

Locate the shock front (same algorithm as `plot_em_emission.py`), select the
slab `[+win_min, +win_max]` skin depths upstream, and compute the 1-D power
spectrum

    P(k_x) = integral P_3D(k_x, k_y, k_z) dk_y dk_z,

i.e. the marginal of the 3-D power spectrum over the transverse wavenumbers,
for the magnetic-field fluctuations carrying the X-mode (delta B_z) and the
O-mode (delta B_y).  By Parseval along y and z this equals the (y, z)-mean of
the squared 1-D rFFT in x, which is what is implemented here.

The compensated spectrum k_x P(k_x) / B_0^2 is plotted log-log so the peak
marks the energy-carrying scale.
"""

import argparse
import json
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import runko
from runko.auto_outdir import resolve_outdir
from runko.mpiio_reader import read_field, read_header

from plot_energy_timeline import compute_upstream
from plot_shock_2d import find_shock_index


CACHE_FILENAME = "em_kspectrum_cache.json"


def load_cache(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def save_cache(path, cache):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)
    os.replace(tmp, path)


def compute_entry(path, conf, up, dx, win_min, win_max, threshold):
    """Read snapshot at `path`, return cache entry dict (or None if no shock)."""
    B0_sq = up["binit"] ** 2

    hdr = read_header(path)
    lap = int(hdr["lap"])
    time_omp = lap * conf.cfl / conf.n_cells_per_skindepth

    n0 = read_field(path, "n0")
    n1 = read_field(path, "n1")
    norm_rho = 2.0 * conf.ppc * conf.io_grid_stride**3
    rho_1d = ((n0 + n1) / norm_rho).mean(axis=(0, 1))
    nx = rho_1d.shape[0]
    del n0, n1

    ix_shock = find_shock_index(rho_1d, threshold)
    if ix_shock == 0:
        return None

    ix_lo = max(0, min(ix_shock + int(math.ceil(win_min / dx)), nx))
    ix_hi = max(0, min(ix_shock + int(math.ceil(win_max / dx)), nx))
    nslab = ix_hi - ix_lo
    if nslab < 2:
        return None

    spectra = {}
    for c in ("by", "bz"):
        delta = read_field(path, c).astype(np.float64)[:, :, ix_lo:ix_hi] - up[c]
        spectra[c] = compute_marginal_spectrum(delta, axis=2) / B0_sq
        del delta

    k = 2.0 * np.pi * np.fft.rfftfreq(nslab, d=dx)

    return {
        "time_omp": time_omp,
        "ix_shock": int(ix_shock),
        "ix_lo": int(ix_lo),
        "ix_hi": int(ix_hi),
        "nslab": int(nslab),
        "k": k.tolist(),
        "zeta_x": spectra["bz"].tolist(),
        "zeta_o": spectra["by"].tolist(),
    }


def compute_marginal_spectrum(delta, axis=2):
    """Marginal 1-D power spectrum along `axis`, normalised so sum_k P_k = <delta^2>.

    Equivalent to integrating |F_3D|^2 over the two transverse wavenumbers
    (Parseval).  Uses np.fft.rfft along `axis`; the one-sided correction
    doubles interior bins but leaves the DC and (when n is even) Nyquist bins
    alone.
    """
    n = delta.shape[axis]
    F = np.fft.rfft(delta, axis=axis) / n
    other = tuple(i for i in range(delta.ndim) if i != axis)
    P = (np.abs(F) ** 2).mean(axis=other)
    P[1:-1] *= 2.0
    if n % 2 != 0:
        P[-1] *= 2.0
    return P


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conf", required=True, help="Path to .ini config file")
    parser.add_argument("--lap", type=int, required=True,
                        help="lap of the flds_<lap>.bin snapshot to analyse")
    parser.add_argument("--win_min", type=float, default=5.0,
                        help="upstream slab inner edge, skin depths (default 5)")
    parser.add_argument("--win_max", type=float, default=25.0,
                        help="upstream slab outer edge, skin depths (default 25)")
    parser.add_argument("--shock_threshold", type=float, default=1.5,
                        help="density threshold n/n_0 for shock localization")
    args = parser.parse_args()

    if args.win_max <= args.win_min:
        sys.exit(f"win_max ({args.win_max}) must be > win_min ({args.win_min})")

    conf = runko.Configuration(args.conf)
    outdir = resolve_outdir(conf)
    up = compute_upstream(conf)
    dx = conf.io_grid_stride / conf.n_cells_per_skindepth  # output cell width in skin depths

    print(f"outdir = {outdir}")
    print(f"binit (= B_0) = {up['binit']:.6g}")
    print(f"upstream Bx={up['bx']:.4g} By={up['by']:.4g} Bz={up['bz']:.4g}")
    print(f"slab = [{args.win_min:.2f}, {args.win_max:.2f}] c/wp upstream of shock "
          f"(dx={dx:.4g} c/wp per coarse cell)")

    cache_path = os.path.join(outdir, CACHE_FILENAME)
    cache = load_cache(cache_path)
    key = str(args.lap)

    if key in cache:
        entry = cache[key]
        print(f"cache hit at {cache_path}: reusing entry for lap {args.lap}")
    else:
        path = os.path.join(outdir, f"flds_{args.lap}.bin")
        if not os.path.exists(path):
            sys.exit(f"snapshot not found: {path}")
        print(f"cache miss: computing from {path}")

        hdr = read_header(path)
        if int(hdr["lap"]) != args.lap:
            print(f"  warning: header lap={hdr['lap']} != --lap={args.lap}")

        entry = compute_entry(
            path, conf, up, dx,
            args.win_min, args.win_max, args.shock_threshold,
        )
        if entry is None:
            sys.exit(f"no usable shock/slab at lap {args.lap} "
                     f"(threshold={args.shock_threshold})")

        cache[key] = entry
        save_cache(cache_path, cache)
        print(f"updated {cache_path}: now {len(cache)} laps")

    k = np.asarray(entry["k"])
    zeta_x = np.asarray(entry["zeta_x"])
    zeta_o = np.asarray(entry["zeta_o"])
    nslab = entry["nslab"]
    L_slab = nslab * dx
    print(f"lap={args.lap}  t={entry['time_omp']:.4g} omega_p^-1  "
          f"ix_shock={entry['ix_shock']}  "
          f"slab=[{entry['ix_lo']}:{entry['ix_hi']}] "
          f"({nslab} cells, ~{L_slab:.1f} c/wp)  "
          f"sum P_O={zeta_o.sum():.4e}  sum P_X={zeta_x.sum():.4e}")

    fig = plt.figure(1, figsize=(3.25, 2.2))
    plt.rc('font', family='serif')
    plt.rc('text', usetex=False)
    plt.rc('xtick', top=True, direction='out', labelsize=7)
    plt.rc('ytick', right=True, direction='out', labelsize=7)
    plt.rc('axes', labelsize=8)
    plt.rc('legend', handlelength=4.0)

    ax = fig.add_subplot(111)
    ax.minorticks_on()
    fig.subplots_adjust(left=0.18, bottom=0.16, right=0.96, top=0.96)

    ax.plot(k[1:], k[1:] * zeta_o[1:],
            color="C0", linewidth=1.0, label=r"O-mode $\delta B_y$")
    ax.plot(k[1:], k[1:] * zeta_x[1:],
            color="C1", linewidth=1.0, label=r"X-mode $\delta B_z$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 30.0)
    ax.set_ylim(1e-5, 1e-1)
    ax.set_xlabel(r"Wavenumber $k\,c/\omega_p$")
    ax.set_ylabel(r"$k_x P(k_x) / B_0^2$")
    ax.legend(frameon=False, fontsize=7)

    slap = str(args.lap).rjust(7, '0')
    base = os.path.join(outdir, f"plot_em_kspectrum_{slap}")
    plt.savefig(base + ".pdf")
    plt.savefig(base + ".png", dpi=300)
    print(f"saved {base}.pdf and .png")


if __name__ == "__main__":
    main()
