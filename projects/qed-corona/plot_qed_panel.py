# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The diagnostic panel of a qed-pic run, from its saved records: the photon and pair spectra and
weight distributions of one record, the escaping spectrum and its angular map, and the time
series of the depth, the energies, the counts, the channel rates and the energy balance of
the antenna against the escaping radiation.

    python plot_qed_panel.py [--outdir corona-pic] [--lap N | --all]
"""

import argparse
import glob
import os
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

SPECIES = ("ph", "e-", "e+")
COLORS = {"ph": "C2", "e-": "C0", "e+": "C1", "esc": "C3"}


def centers(edges):
    return np.sqrt(edges[:-1] * edges[1:])


def pos(a):
    """positive values only, for a log axis; the rest is left out of the line"""
    a = np.asarray(a, dtype=float)
    return np.where(a > 0, a, np.nan)


def record_laps(outdir):
    return sorted(int(os.path.basename(f)[4:-4]) for f in glob.glob(os.path.join(outdir, "qed_*.npz"))
                  if os.path.basename(f)[4:-4].isdigit())


def load_series(outdir):
    """the time series qed_timeseries.npz, written when the run finished; None before that"""
    path = os.path.join(outdir, "qed_timeseries.npz")
    return dict(np.load(path, allow_pickle=True)) if os.path.exists(path) else None


def energy_series(outdir, conf, w_unit):
    """field and kinetic energies per record lap, in the depth-weight units of the QED records:
    the code-unit energy of a numerical particle is m0_num, of weight w_lep"""
    n_t = conf.n_tiles[0]; n_c = conf.n_cells_per_tile[0]
    vol = (n_t * n_c) ** 3
    m0_num = conf.m0 * abs(conf.q0)
    scale = w_unit / m0_num / n_t ** 2                  # code energy -> depth-weight units per tile edge
    out = {}
    for key, fn in (("U_B", "average_B_energy_density.txt"), ("U_E", "average_E_energy_density.txt")):
        d = np.loadtxt(os.path.join(outdir, fn), ndmin=2)
        out[key] = (d[:, 0].astype(int), d[:, 1] * vol * scale)
    d = np.loadtxt(os.path.join(outdir, "average_kinetic_energy.txt"), ndmin=2)
    n_lep = conf.ppc * vol
    out["U_kin"] = (d[:, 0].astype(int), (d[:, 1] + d[:, 2]) * n_lep * scale)   # leptons only: photons carry their own record
    return out


def panel(outdir, lap, conf, series, energies, dt_box, save):
    R = np.load(os.path.join(outdir, f"qed_{lap}.npz"))
    x = centers(R["x_bins"]); p = centers(R["p_bins"]); w = centers(R["w_bins"]); mu = 0.5 * (R["mu_bins"][1:] + R["mu_bins"][:-1])
    dlx = np.diff(np.log(R["x_bins"])); dlp = np.diff(np.log(R["p_bins"])); dlw = np.diff(np.log(R["w_bins"]))
    names = [str(n) for n in R["names"]]
    n_t = conf.n_tiles[0]

    fig, axs = plt.subplots(5, 4, figsize=(16, 17))
    for ax in axs.flat:
        ax.minorticks_on(); ax.set_yscale("log")

    ax = axs[0, 0]; ax.set_title("photons"); ax.set_xlabel(r"$\log x$"); ax.set_ylabel(r"$x\,\ell_x \propto x\,dn_x/d\log x$")
    ax.plot(np.log10(x), R["ph_xw"] / dlx, color=COLORS["ph"]); ax.set_ylim(1e-6, 1e3)
    ax = axs[0, 1]; ax.set_title("pairs"); ax.set_xlabel(r"$\log p$"); ax.set_ylabel(r"$p\,d\tau/dp$")
    ax.plot(np.log10(p), R["lep_pw"] / dlp, color=COLORS["e-"]); ax.set_ylim(1e-4, 1e2)
    ax = axs[0, 2]; ax.set_title("escaping photons"); ax.set_xlabel(r"$\log x$"); ax.set_ylabel(r"$x\,\ell_x$ per crossing")
    ax.plot(np.log10(x), R["esc"] * x / dlx, color=COLORS["esc"]); ax.set_ylim(1e-6, 1e3)
    ax = axs[1, 0]; ax.set_xlabel(r"$\log x$"); ax.set_ylabel(r"$dN_x/d\log x$"); ax.plot(np.log10(x), R["ph_xn"] / dlx, color=COLORS["ph"])
    ax = axs[1, 1]; ax.set_xlabel(r"$\log p$"); ax.set_ylabel(r"$dN_p/d\log p$"); ax.plot(np.log10(p), R["lep_pn"] / dlp, color=COLORS["e-"])
    ax = axs[2, 0]; ax.set_xlabel(r"$\log w$"); ax.set_ylabel(r"$dN_x/d\log w$"); ax.plot(np.log10(w), R["ph_wn"] / dlw, color=COLORS["ph"])
    ax.axvline(np.log10(conf.w0), color=COLORS["e-"], lw=0.8, label=r"$w_\pm$"); ax.legend(loc="upper left", fontsize=8)
    ax = axs[2, 1]; ax.set_xlabel(r"$\log x$"); ax.set_ylabel(r"$\langle w \rangle(x)$"); ax.set_title("mean photon weight")
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(np.log10(x), R["ph_xw"] / x / np.maximum(R["ph_xn"], 1) * n_t ** 2, color=COLORS["ph"])
    ax.axhline(conf.w0, color=COLORS["e-"], lw=0.8)

    ax = axs[3, 0]; ax.set_yscale("linear"); ax.set_xlabel(r"$\log w$"); ax.set_ylabel(r"$\log x$")
    H = R["ph_xw2"]
    if H.max() > 0:
        ax.pcolormesh(np.log10(R["w_bins"]), np.log10(R["x_bins"]), H, norm=LogNorm(vmin=H[H > 0].min(), vmax=H.max()), cmap="viridis")
    ax = axs[3, 1]; ax.set_yscale("linear"); ax.set_xlabel(r"$\mu = k_z/x$"); ax.set_ylabel(r"$\log x$"); ax.set_title(r"escaping $d(x\ell_x)/d\mu$")
    H = R["esc_mu"] * x[:, None]
    if H.max() > 0:
        ax.pcolormesh(R["mu_bins"], np.log10(R["x_bins"]), H, norm=LogNorm(vmin=H[H > 0].min(), vmax=H.max()), cmap="magma")

    if series is not None:
        S = series; t = S["t"]; laps = S["lap"]
        ax = axs[0, 3]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$\tau_T$")
        ax.plot(t, S["tau"], color="C0", label="box"); ax.plot(t, S["tau_min"], color="C1", lw=0.7, label="tile min/max"); ax.plot(t, S["tau_max"], color="C1", lw=0.7)
        ax.legend(fontsize=8)
        ax = axs[1, 2]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$U_x/U_\pm$")
        ax.plot(t, S["E"][:, 0] / np.maximum(S["E"][:, 1] + S["E"][:, 2], 1e-300), color="C0")
        ax = axs[1, 3]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"energy $\ell$ (per tile edge)")
        for i, s in enumerate(SPECIES):
            ax.plot(t, S["E"][:, i], color=COLORS[s], label=s)
        ax.legend(fontsize=8)
        ax = axs[2, 2]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$W/W_0$")
        W0 = S["W"][0, 1] + S["W"][0, 2]
        for i, s in enumerate(SPECIES):
            ax.plot(t, S["W"][:, i] / W0, color=COLORS[s], label=s)
        ax.plot(t, pos(np.cumsum(S["esc_W"]) / W0), color=COLORS["esc"], ls="dashed", label="escaped"); ax.legend(fontsize=8)
        ax = axs[2, 3]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$N_{LP}$")
        for i, s in enumerate(SPECIES):
            ax.plot(t, S["N"][:, i], color=COLORS[s], label=s)
        ax.legend(fontsize=8)
        ax = axs[3, 2]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel("events per crossing")
        for k, nm in enumerate(names):
            ax.plot(t, pos(S["counts"][:, k]), label=nm, lw=0.9)
        ax.legend(fontsize=6, ncol=2)
        ax = axs[4, 2]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel("weight moved per crossing")
        for k, nm in enumerate(names):
            ax.plot(t, pos(S["moved"][:, k]), label=nm, lw=0.9)
        ax = axs[4, 1]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$\langle\gamma\rangle$, $\langle x\rangle$")
        for i, s in enumerate(SPECIES):
            ax.plot(t, S["E"][:, i] / np.maximum(S["W"][:, i], 1e-300), color=COLORS[s], label=s)
        ax.legend(fontsize=8)

        # the energy balance: the antenna feeds the fields and the plasma; what does not
        # accumulate leaves as radiation. ell_inj = dU_tot/dt + ell_esc, sampled at the records
        ax = axs[4, 0]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"field energy (per tile edge)")
        for key, c in (("U_B", "C0"), ("U_E", "C1"), ("U_kin", "C2")):
            L, U = energies[key]; ax.plot(L * dt_box, U, color=c, label=key)
        ax.legend(fontsize=8)
        U_tot = np.zeros_like(t)
        for key in ("U_B", "U_E"):
            L, U = energies[key]; U_tot += np.interp(laps, L, U)
        U_tot += S["E"].sum(axis=1)
        ell_inj = np.gradient(U_tot, t) + S["ell_esc"] if len(t) > 1 else np.zeros_like(t)
        ax = axs[3, 3]; ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$\ell_\mathrm{inj}$, $\ell_\mathrm{esc}$")
        ax.plot(t, pos(ell_inj), color="C0", label=r"$\ell_\mathrm{inj}$ (antenna)"); ax.plot(t, pos(S["ell_esc"]), color=COLORS["esc"], label=r"$\ell_\mathrm{esc}$")
        ax.legend(fontsize=8)
        ax = axs[4, 3]; ax.set_yscale("linear"); ax.set_xlabel(r"$t$ ($H/c$)"); ax.set_ylabel(r"$\ell_\mathrm{inj}/\ell_\mathrm{esc}$")
        ax.plot(t, ell_inj / np.maximum(S["ell_esc"], 1e-300), color="C0"); ax.axhline(1.0, lw=0.5); ax.set_ylim(0, 4)
        for ax in (axs[0, 3], axs[1, 2], axs[1, 3], axs[2, 2], axs[2, 3], axs[3, 2], axs[3, 3], axs[4, 0], axs[4, 1], axs[4, 2], axs[4, 3]):
            ax.axvline(lap * dt_box, color="k", lw=0.5, ls="dotted")

    fig.suptitle(f"lap {lap}, t = {lap * dt_box:.3f} H/c")
    fig.tight_layout()
    fig.savefig(save, dpi=110); plt.close(fig)
    print(f"saved {save}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="diagnostic panel of a qed-pic run")
    parser.add_argument("--outdir", default="corona-pic")
    parser.add_argument("--lap", type=int, default=None, help="record lap (default: the last)")
    parser.add_argument("--all", action="store_true", help="one panel per record")
    args = parser.parse_args()

    with open(os.path.join(args.outdir, "config.pkl"), "rb") as f:
        conf = pickle.load(f)
    n_t = conf.n_tiles[0]
    dt_box = conf.cfl / conf.n_cells_per_tile[0] / n_t
    laps = record_laps(args.outdir)
    series = load_series(args.outdir)
    energies = energy_series(args.outdir, conf, conf.w0)
    todo = laps if args.all else [args.lap if args.lap is not None else laps[-1]]
    for lap in todo:
        panel(args.outdir, lap, conf, series, energies, dt_box, os.path.join(args.outdir, f"qed_panel_{lap}.png"))
