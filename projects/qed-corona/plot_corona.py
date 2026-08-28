# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Six panels of a corona run: pair balance, energy budget, LP counts, channel rates, photon spectra,
lepton distribution against a Maxwell-Juettner fit. Reads the corona.npz the driver writes."""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import kv
from scipy.optimize import brentq


def setup_rc():
    plt.rc("font", family="serif")
    plt.rc("text", usetex=False)
    plt.rc("xtick", top=True, direction="out", labelsize=7)
    plt.rc("ytick", right=True, direction="out", labelsize=7)
    plt.rc("axes", labelsize=8)
    plt.rc("legend", handlelength=3.0)


def load(fname):
    d = np.load(fname)
    rec = {k: d[k] for k in ("t", "tau", "Wg", "Eg", "Epm", "N", "ell_esc", "ell_in", "gbar", "dtau_cre", "dtau_ann", "dtau_inj")}
    rec["rates"] = {k: d["rates"][i] for i, k in enumerate(d["rate_names"])}
    sc = list(d["scalars"])
    return dict(rec=rec, cx=d["cx"], ellx=d["ellx"], lep_p=d["lep_p"], lep_w=d["lep_w"], ph_x=d["ph_x"], ph_w=d["ph_w"],
                tau=sc[0], ell_esc=sc[1], ell_inj=sc[2], x_peak=sc[3], x_ssa=sc[4], slope=sc[5], gbar=sc[6], nlate=int(sc[7]),
                nsnap=int(sc[8]), w_lep=sc[9], x_min=sc[10], x_unit=sc[11], tau_ssa_min=sc[12], cre=sc[13], ann=sc[14], inj=sc[15],
                a_res=sc[16], x_ref=sc[17])


def mean_gamma(T):
    return 3 * T + kv(1, 1 / T) / kv(2, 1 / T)


def plot_corona(res, fname):
    setup_rc()
    rec = res["rec"]
    t = rec["t"]
    fig = plt.figure(figsize=(7.0, 4.6))
    gs = plt.GridSpec(2, 3)
    gs.update(wspace=0.36, hspace=0.34)
    axs = [[plt.subplot(gs[i, j]) for j in range(3)] for i in range(2)]
    fig.subplots_adjust(left=0.07, bottom=0.10, right=0.98, top=0.97)
    # pair budget: the depth, and the three terms of d tau_T/dt -- injection, gamma-gamma creation,
    # annihilation -- with their net against the measured derivative, all per crossing
    ax = axs[0][0]
    ax.plot(t, rec["tau"], "k", lw=1, label=r"$\tau_\mathrm{T}$")
    nk = min(20, len(t)); k = np.ones(nk) / nk
    sm = lambda y: np.convolve(y, k, "same")
    inj, cre, ann = rec["dtau_inj"], rec["dtau_cre"], rec["dtau_ann"]
    ax.plot(t, sm(inj), "C0", lw=0.8, label="injection")
    ax.plot(t, sm(cre), "C2", lw=0.8, label=r"$\gamma\gamma$ creation")
    ax.plot(t, -sm(ann), "C3", lw=0.8, label=r"$-$annihilation")
    ax.plot(t, sm(inj + cre - ann), "C4", lw=0.8, label="net")
    if len(t) > 2:
        ax.plot(t, sm(np.gradient(rec["tau"], t)), "k", lw=0.6, ls="--", label=r"measured $\mathrm{d}\tau_\mathrm{T}/\mathrm{d}t$")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlabel(r"Time $ct/H$"); ax.set_ylabel(r"$\tau_\mathrm{T}$ and $\mathrm{d}\tau_\mathrm{T}/\mathrm{d}t \times H/c$")
    top = max(1.05 * rec["tau"].max(), 1e-3); ax.set_ylim(-0.4 * top, top)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2)
    # energy budget
    ax = axs[0][1]
    ax.plot(t, rec["ell_esc"] / res["ell_inj"], "C0", lw=1, label=r"$\ell_\mathrm{esc}/\ell_\mathrm{inj}$")
    ax.plot(t, rec["Eg"] / np.maximum(rec["tau"], 1e-300) / 10, "C1", lw=1, label=r"$\eta_x/10$")
    ax.plot(t, rec["gbar"], "C4", lw=1, label=r"$\langle\gamma\rangle$")
    ax.axhline(1.0, color="k", ls=":", lw=0.6)
    ax.set_xlabel(r"Time $ct/H$"); ax.set_ylabel("Energy budget"); ax.set_ylim(0, 3)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    # LP counts
    ax = axs[0][2]
    for i, (lab, col) in enumerate((("photons", "C0"), ("electrons", "C1"), ("positrons", "C2"))):
        ax.plot(t, rec["N"][:, i], col, lw=1, label=lab)
    ax.set_yscale("log"); ax.set_xlabel(r"Time $ct/H$"); ax.set_ylabel("Number of LPs")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    # channel rates per crossing
    ax = axs[1][0]
    r = rec["rates"]
    for kk, col in (("compton-", "C0"), ("compton+", "C0"), ("photann", "C2"), ("pairann", "C3"), ("ssa-", "C1"), ("ssa+", "C1"), ("syn-", "C4"), ("syn+", "C4")):
        if kk in r:
            ax.plot(t, np.maximum(r[kk], 1e-3), col, lw=0.8, ls="-" if kk.endswith("-") or kk in ("photann", "pairann") else "--",
                    label=kk if not kk.endswith("+") else None)
    ax.set_yscale("log"); ax.set_xlabel(r"Time $ct/H$"); ax.set_ylabel("Events per crossing")
    ax.legend(frameon=False, fontsize=6, loc="lower right", ncol=2)
    # photon spectra: escaping (late average) and inside the box (end)
    ax = axs[1][1]
    ok = res["ellx"] > 0
    top = res["ellx"].max() if ok.any() else 1.0
    ax.plot(res["cx"][ok], res["ellx"][ok], "k", lw=1, label=r"escaping $x\,\ell_x$")
    bins = np.logspace(-8, 2, 101)
    h, e = np.histogram(res["ph_x"], bins=bins, weights=res["ph_x"] * res["ph_w"])
    c = np.sqrt(e[1:] * e[:-1]); hin = h / np.diff(np.log(e)); okk = hin > 0
    scale = top / max(hin.max(), 1e-300)
    ax.plot(c[okk], hin[okk] * scale, "C0", lw=0.8, ls="--", label=r"inside, $x\,\mathrm{d}E_\gamma/\mathrm{d}\ln x$ (scaled)")
    ax.axvline(res["x_min"], color="k", lw=0.4, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-8, 1e2); ax.set_ylim(max(top * 1e-6, 1e-10), top * 3)
    ax.set_xlabel(r"Photon energy $x$"); ax.set_ylabel(r"$x\,\ell_x$")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    # lepton distribution against a Maxwell-Juettner fit at the late mean energy
    ax = axs[1][2]
    bins = np.logspace(-2, 2, 60)
    h, e = np.histogram(res["lep_p"], bins=bins, weights=res["lep_w"])
    c = np.sqrt(e[1:] * e[:-1]); d = c * h / np.diff(np.log(e)) / res["nsnap"]
    ax.plot(c, d, "k", lw=1, label="simulation")
    gb = res["gbar"]
    if 1.002 < gb < 1e3:
        T = brentq(lambda T: mean_gamma(T) - gb, 1e-2, 3e2)
        ax.plot(c, res["tau"] * c ** 3 * np.exp(-np.sqrt(1 + c * c) / T) / (T * kv(2, 1 / T)), "k", lw=0.8, ls=":",
                label=rf"Maxwell–Jüttner, $\theta = {T:.2f}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-2, 1e2); ax.set_ylim(max(d.max() * 1e-5, 1e-8), d.max() * 3)
    ax.set_xlabel(r"Lepton momentum $p$"); ax.set_ylabel(r"$p^2\,\mathrm{d}\tau_\mathrm{T}/\mathrm{d}p$")
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    for row in axs:
        for ax in row:
            ax.minorticks_on()
    plt.savefig(fname)
    plt.close(fig)
    print("saved", fname)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="corona")
    args = ap.parse_args()
    res = load(os.path.join(args.outdir, "corona.npz"))
    print(f"tau_T = {res['tau']:.3f}, ell_esc/ell_inj = {res['ell_esc']/res['ell_inj']:.3f}, x_peak = {res['x_peak']:.3g}, "
          f"photon index {2-res['slope']:.2f}, x_SSA = {res['x_unit']:.2e}, <gamma> = {res['gbar']:.3f}, "
          f"pair budget (inj + gg)/ann = {(res['inj']+res['cre'])/max(res['ann'], 1e-300):.3f}")
    plot_corona(res, os.path.join(args.outdir, "corona.pdf"))
