# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""The photon LPs inside the box at the end of a corona run: weight against energy, and the LP
count, weight and LPs per unit weight per energy bin against the resolution law."""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_corona import load, setup_rc


def plot_weights(res, fname):
    setup_rc()
    x, w = res["ph_x"], res["ph_w"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6))
    fig.subplots_adjust(left=0.08, bottom=0.18, right=0.98, top=0.95, wspace=0.3)
    ax = axs[0]
    lx, lw = np.log10(np.maximum(x, 1e-300)), np.log10(np.maximum(w, 1e-300))
    h = ax.hist2d(lx, lw, bins=[np.linspace(-6, 1, 71), np.linspace(lw.min() - 0.1, lw.max() + 0.1, 76)],
                  weights=w, cmap="viridis", norm=matplotlib.colors.LogNorm())
    plt.colorbar(h[3], ax=ax, label="weight per bin")
    ax.set_xlabel(r"$\log_{10} x$"); ax.set_ylabel(r"$\log_{10} w$")
    ax = axs[1]
    edges = np.logspace(-6, 1, 15); c = np.sqrt(edges[1:] * edges[:-1])
    NN = np.histogram(x, bins=edges)[0]; WW = np.histogram(x, bins=edges, weights=w)[0]; ok = WW > 0
    ax.plot(c[ok], NN[ok] / WW[ok], "k", lw=1, label=r"$N_\mathrm{LP}/W$")
    if ok.any():
        ref = (NN[ok] / WW[ok])[len(c[ok]) // 2] / (c[ok][len(c[ok]) // 2] / res["x_ref"]) ** res["a_res"]
        ax.plot(c, ref * (c / res["x_ref"]) ** res["a_res"], "k:", lw=0.8, label=rf"$\propto x^{{{res['a_res']:g}}}$ (resolution function)")
    ax.plot(c[ok], NN[ok], "C0", lw=1, label=r"$N_\mathrm{LP}$")
    ax.plot(c[ok], WW[ok], "C1", lw=1, label=r"$W$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"Photon energy $x$")
    ax.legend(frameon=False, fontsize=6)
    for ax in axs:
        ax.minorticks_on()
    plt.savefig(fname)
    plt.close(fig)
    print("saved", fname)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="corona")
    args = ap.parse_args()
    plot_weights(load(os.path.join(args.outdir, "corona.npz")), os.path.join(args.outdir, "corona_weights.pdf"))
