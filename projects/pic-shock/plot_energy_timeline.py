"""Plot time evolution of average energy density in each EM field component."""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import glob
import os
import math

import runko
from runko.auto_outdir import resolve_outdir
from runko.mpiio_reader import read_header, read_field


# -- physics helpers -----------------------------------------------------------

def compute_upstream(conf):
    """Compute upstream field values and binit (same formulas as pic.py)."""
    oppc = 2 * conf.ppc
    omp  = conf.cfl / conf.c_omp
    qe   = -(omp**2 * conf.upstream_gamma) / (0.5 * oppc * (1.0 + conf.m0 / conf.m1))
    m0   = conf.m0 * abs(qe)
    m1   = conf.m1 * abs(qe)
    binit = math.sqrt(conf.upstream_gamma * oppc * 0.5 * conf.cfl**2 * m0 * (1.0 + m0/m1) * conf.sigma)
    beta  = math.sqrt(1.0 - 1.0 / conf.upstream_gamma**2)
    return {
        "binit": binit,
        "gamma": conf.upstream_gamma,
        "ex": 0.0,
        "ey": -beta * binit * conf.bz_proj,
        "ez": +beta * binit * conf.by_proj,
        "bx": binit * conf.bx_proj,
        "by": binit * conf.by_proj,
        "bz": binit * conf.bz_proj,
    }


# -- main ----------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Energy density timeline for PIC shock")
    parser.add_argument("--conf", required=True, help="Path to .ini config file")
    args = parser.parse_args()

    # -- config ----------------------------------------------------------------
    conf = runko.Configuration(args.conf)
    outdir = resolve_outdir(conf)

    cfl   = conf.cfl
    c_omp = conf.c_omp
    up    = compute_upstream(conf)
    U0    = up["binit"]**2  # reference energy density

    print(f"outdir = {outdir}")
    print(f"binit  = {up['binit']:.6g},  U0 = {U0:.6g}")
    print(f"upstream: Ex={up['ex']:.4g}  Ey={up['ey']:.4g}  Ez={up['ez']:.4g}  "
          f"Bx={up['bx']:.4g}  By={up['by']:.4g}  Bz={up['bz']:.4g}")

    # -- load field snapshots --------------------------------------------------
    fld_files = sorted(glob.glob(os.path.join(outdir, "flds_*.bin")),
                       key=lambda f: int(os.path.basename(f).split("_")[1].split(".")[0]))
    print(f"found {len(fld_files)} field snapshots")

    components = ["ex", "ey", "ez", "bx", "by", "bz"]
    snap_laps = []
    snap_energy = {c: [] for c in components}

    for fpath in fld_files:
        hdr = read_header(fpath)
        snap_laps.append(hdr["lap"])
        for c in components:
            arr = read_field(fpath, c).astype(np.float64)
            snap_energy[c].append(np.mean((arr - up[c])**2))

    snap_laps  = np.array(snap_laps)
    snap_time  = snap_laps * cfl / c_omp  # plasma time ω_p⁻¹
    for c in components:
        snap_energy[c] = np.array(snap_energy[c]) / U0

    # -- load CSV time series --------------------------------------------------
    ke_data = np.loadtxt(os.path.join(outdir, "average_kinetic_energy.txt"))
    B_data  = np.loadtxt(os.path.join(outdir, "average_B_energy_density.txt"))
    E_data  = np.loadtxt(os.path.join(outdir, "average_E_energy_density.txt"))

    # CSV field energies are: sum(F²) / (2 * N_cells)
    # i.e., volume-averaged electromagnetic energy density in code units (F²/2).
    csv_factor = 1.0 / 2.0  # CSV_value = csv_factor * <F²>

    E2_up = up["ex"]**2 + up["ey"]**2 + up["ez"]**2
    B2_up = up["bx"]**2 + up["by"]**2 + up["bz"]**2

    csv_time_E = E_data[:, 0] * cfl / c_omp
    csv_E      = np.abs(E_data[:, 1] - csv_factor * E2_up) / (csv_factor * U0)

    csv_time_B = B_data[:, 0] * cfl / c_omp
    csv_B      = np.abs(B_data[:, 1] - csv_factor * B2_up) / (csv_factor * U0)

    # CSV KE is <γ-1> per species (no mass factor); upstream value is γ_up - 1
    KE_up = (up["gamma"] - 1.0) * 2  # both species identical, sum of 2
    csv_time_ke = ke_data[:, 0] * cfl / c_omp
    csv_ke      = np.abs((ke_data[:, 1] + ke_data[:, 2]) - KE_up) / U0

    # -- figure setup ----------------------------------------------------------
    fig = plt.figure(1, figsize=(3.25, 2.2))

    plt.rc('font', family='serif')
    plt.rc('text', usetex=False)
    plt.rc('xtick', top=True, direction='out', labelsize=7)
    plt.rc('ytick', right=True, direction='out', labelsize=7)
    plt.rc('axes', labelsize=8)
    plt.rc('legend', handlelength=4.0)

    ax = fig.add_subplot(111)
    ax.minorticks_on()
    fig.subplots_adjust(left=0.18, bottom=0.16, right=0.96, top=0.92)

    # -- plot per-component δF energy from snapshots ---------------------------
    lstyles = {"x": "solid", "y": "dashed", "z": "dotted"}

    colors = {"e": "C0", "b": "C1"}

    for c in components:
        fb   = c[0]          # 'e' or 'b'
        comp = c[1]          # 'x', 'y', or 'z'
        ax.plot(snap_time, snap_energy[c],
                color=colors[fb], linestyle=lstyles[comp], linewidth=0.8,
                label=rf"$\delta {fb.upper()}_{comp}$")

    # -- plot CSV average energies ---------------------------------------------
    ax.plot(csv_time_E, csv_E, color="C0",   linewidth=1.3, alpha=0.7, label=r"$\delta\langle E^2 \rangle$")
    ax.plot(csv_time_B, csv_B, color="C1",   linewidth=1.3, alpha=0.7, label=r"$\delta\langle B^2 \rangle$")
    ax.plot(csv_time_ke, csv_ke, color="C2", linewidth=1.3, alpha=0.7, label=r"$|\delta K|$")

    # -- axes ------------------------------------------------------------------
    ax.set_yscale("log")
    ax.set_xlabel(r"Time $t~(\omega_p^{-1})$")
    ax.set_ylabel(r"$\langle \delta F^2 \rangle / B_0^2$")
    ax.set_ylim((1e-11, 1e2))
    ax.set_xlim(snap_time[0], snap_time[-1])
    ax.legend(frameon=False, fontsize=5, ncol=3)

    # -- save ------------------------------------------------------------------
    base = os.path.join(outdir, "plot_energy_timeline")
    plt.savefig(base + ".pdf")
    plt.savefig(base + ".png", dpi=300)
    print(f"saved {base}.pdf and .png")
