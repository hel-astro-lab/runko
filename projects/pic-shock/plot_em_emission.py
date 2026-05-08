"""EM emission timeline measured in an upstream slab in front of the shock.

For each `flds_<lap>.bin` snapshot, locate the shock front (same algorithm as
`plot_shock_2d.py --track_shock`), select a slab `[+win_min, +win_max]` skin
depths upstream of it, and compute:

  - per-component fluctuation energies <(F - F_up)^2> for F in {Ex,Ey,Ez,Bx,By,Bz}
  - X-mode magnetic energy fraction zeta_X = <(Bz - Bz_up)^2> / B_0^2
  - O-mode magnetic energy fraction zeta_O = <(By - By_up)^2> / B_0^2

Results are cached as a CSV in the job dir; reruns skip already-computed laps
and only process new snapshots.
"""

import argparse
import glob
import math
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import runko
from runko.auto_outdir import resolve_outdir
from runko.mpiio_reader import read_field, read_header

from plot_energy_timeline import compute_upstream
from plot_shock_2d import find_shock_index


# CSV columns (header line uses '#' so np.loadtxt ignores it)
CSV_COLUMNS = [
    "lap", "time_omp", "ix_shock", "ix_lo", "ix_hi", "nslab",
    "dEx2", "dEy2", "dEz2", "dBx2", "dBy2", "dBz2",
    "zeta_X", "zeta_O",
]

EM_COMPONENTS = ["ex", "ey", "ez", "bx", "by", "bz"]


def lap_from_path(path):
    return int(os.path.basename(path).split("_")[1].split(".")[0])


def load_csv(path):
    """Return (rows, cached_laps) where rows is a (N, len(CSV_COLUMNS)) array."""
    if not os.path.exists(path):
        return np.empty((0, len(CSV_COLUMNS))), set()
    rows = np.loadtxt(path, ndmin=2)
    if rows.size == 0:
        return np.empty((0, len(CSV_COLUMNS))), set()
    if rows.shape[1] != len(CSV_COLUMNS):
        raise ValueError(
            f"{path}: expected {len(CSV_COLUMNS)} columns, got {rows.shape[1]} — "
            f"delete the file to rebuild from scratch"
        )
    return rows, {int(lap) for lap in rows[:, 0]}


def save_csv(path, rows):
    """Sort by lap and write atomically."""
    rows = rows[np.argsort(rows[:, 0])]
    tmp = path + ".tmp"
    header = " ".join(CSV_COLUMNS)
    fmt = ["%d", "%.6g", "%d", "%d", "%d", "%d"] + ["%.6e"] * 8
    np.savetxt(tmp, rows, header=header, fmt=fmt)
    os.replace(tmp, path)


def analyze_snapshot(path, conf, up, dx, win_min, win_max, threshold):
    """Compute one row of CSV data from a flds_<lap>.bin snapshot.

    Returns the row as a list, or None if the slab is unusable
    (no shock yet, or shock too close to the upstream boundary).
    """
    hdr = read_header(path)
    lap = int(hdr["lap"])
    time_omp = lap * conf.cfl / conf.c_omp

    # Density profile — used only for shock localization.
    n0 = read_field(path, "n0")
    n1 = read_field(path, "n1")
    norm_rho = 2.0 * conf.ppc * conf.stride**3
    rho_1d = ((n0 + n1) / norm_rho).mean(axis=(0, 1))
    nx = rho_1d.shape[0]
    del n0, n1

    ix_shock = find_shock_index(rho_1d, threshold)
    if ix_shock == 0:
        print(f"  lap {lap}: no shock above n/n_0 > {threshold}, skipping")
        return None

    ix_lo = ix_shock + int(math.ceil(win_min / dx))
    ix_hi = ix_shock + int(math.ceil(win_max / dx))
    ix_lo = max(0, min(ix_lo, nx))
    ix_hi = max(0, min(ix_hi, nx))
    nslab = ix_hi - ix_lo
    if nslab < 1:
        print(f"  lap {lap}: slab past upstream boundary "
              f"(ix_shock={ix_shock}, nx={nx}), skipping")
        return None

    fluct = {}
    for c in EM_COMPONENTS:
        arr = read_field(path, c).astype(np.float64)
        fluct[c] = float(np.mean((arr[:, :, ix_lo:ix_hi] - up[c])**2))
        del arr

    B0_sq = up["binit"]**2
    zeta_X = fluct["bz"] / B0_sq  # <(Bz - Bz_up)^2> / B_0^2
    zeta_O = fluct["by"] / B0_sq  # <(By - By_up)^2> / B_0^2

    print(f"  lap {lap}: ix_shock={ix_shock} slab=[{ix_lo}:{ix_hi}] "
          f"({nslab} cells, ~{nslab*dx:.1f} c/wp)  "
          f"zeta_X={zeta_X:.3e} zeta_O={zeta_O:.3e}")

    return [lap, time_omp, ix_shock, ix_lo, ix_hi, nslab,
            fluct["ex"], fluct["ey"], fluct["ez"],
            fluct["bx"], fluct["by"], fluct["bz"],
            zeta_X, zeta_O]


def make_plot(rows, outdir, U0):
    """Single panel: per-component fluctuations and X/O mode energies."""
    if rows.shape[0] == 0:
        print("no data to plot")
        return

    cols = {name: rows[:, i] for i, name in enumerate(CSV_COLUMNS)}

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

    t = cols["time_omp"]

    field_color = {"e": "C0", "b": "C1"}
    comp_lstyle = {"x": "solid", "y": "dashed", "z": "dashdot"}
    for c in EM_COMPONENTS:
        fb, comp = c[0], c[1]
        key = f"d{fb.upper()}{comp}2"
        ax.plot(t, cols[key] / U0,
                color=field_color[fb], linestyle=comp_lstyle[comp], linewidth=0.8,
                label=rf"$\langle\delta {fb.upper()}_{comp}^2\rangle$")

    ax.plot(t, cols["zeta_X"], color="C2", linewidth=1.2,
            label=r"$\zeta_X$")
    ax.plot(t, cols["zeta_O"], color="C3", linewidth=1.2,
            label=r"$\zeta_O$")

    ax.set_yscale("log")
    ax.set_xlabel(r"Time $t~(\omega_p^{-1})$")
    ax.set_ylabel(r"$\langle\delta F^2\rangle / B_0^2$")
    ax.set_xlim(t.min(), t.max())
    ax.legend(frameon=False, fontsize=5, ncol=2)

    base = os.path.join(outdir, "plot_em_emission")
    plt.savefig(base + ".pdf")
    plt.savefig(base + ".png", dpi=300)
    print(f"saved {base}.pdf and .png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conf", required=True, help="Path to .ini config file")
    parser.add_argument("--win_min", type=float, default=5.0,
                        help="upstream slab inner edge, skin depths (default 5)")
    parser.add_argument("--win_max", type=float, default=25.0,
                        help="upstream slab outer edge, skin depths (default 25)")
    parser.add_argument("--shock_threshold", type=float, default=1.5,
                        help="density threshold n/n_0 for shock localization")
    parser.add_argument("--csv", default="em_emission_timeline.csv",
                        help="cache filename inside outdir")
    args = parser.parse_args()

    if args.win_max <= args.win_min:
        sys.exit(f"win_max ({args.win_max}) must be > win_min ({args.win_min})")

    conf = runko.Configuration(args.conf)
    outdir = resolve_outdir(conf)
    up = compute_upstream(conf)
    U0 = up["binit"]**2
    dx = conf.stride / conf.c_omp  # output cell width in skin depths

    print(f"outdir = {outdir}")
    print(f"binit (= B_0) = {up['binit']:.6g}")
    print(f"upstream Bx={up['bx']:.4g} By={up['by']:.4g} Bz={up['bz']:.4g}")
    print(f"slab = [{args.win_min:.2f}, {args.win_max:.2f}] c/wp upstream of shock "
          f"(dx={dx:.4g} c/wp per coarse cell)")

    csv_path = os.path.join(outdir, args.csv)
    rows, cached_laps = load_csv(csv_path)
    print(f"cache: {csv_path} ({len(cached_laps)} laps)")

    fld_files = sorted(
        glob.glob(os.path.join(outdir, "flds_*.bin")),
        key=lap_from_path,
    )
    new_files = [f for f in fld_files if lap_from_path(f) not in cached_laps]
    print(f"snapshots: {len(fld_files)} total, {len(new_files)} new, "
          f"{len(fld_files) - len(new_files)} cached")

    new_rows = []
    for fpath in new_files:
        row = analyze_snapshot(
            fpath, conf, up, dx,
            args.win_min, args.win_max, args.shock_threshold,
        )
        if row is not None:
            new_rows.append(row)

    if new_rows:
        rows = np.vstack([rows, np.array(new_rows)])
        save_csv(csv_path, rows)
        print(f"updated {csv_path}: now {rows.shape[0]} rows")
    else:
        print("no new rows to write")

    make_plot(rows, outdir, U0)


if __name__ == "__main__":
    main()
