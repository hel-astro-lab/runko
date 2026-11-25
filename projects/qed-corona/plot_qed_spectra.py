# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The weighted momentum spectra of the pairs and the photons of a qed-pic run over time, from
the pspectra_<lap>.bin snapshots summed over the tiles: u dN/du per species, colored by lap.

    python plot_qed_spectra.py [--outdir corona-pic]
"""

import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from runko.mpiio_spectra_reader import (read_spectra_header, read_spectra_snapshot,
                                        u_bin_edges, u_bin_centers)


SPECIES = {0: "e-", 1: "e+", 2: "ph"}


def spectra_laps(outdir):
    return sorted(int(os.path.basename(f)[9:-4])
                  for f in glob.glob(os.path.join(outdir, "pspectra_*.bin")))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="momentum spectra of a qed-pic run over time")
    parser.add_argument("--outdir", default="corona-pic")
    parser.add_argument("--every", type=int, default=1, help="plot every n-th snapshot")
    args = parser.parse_args()

    laps = spectra_laps(args.outdir)[::args.every]

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    cmap = plt.get_cmap("viridis")
    for k, lap in enumerate(laps):
        path = os.path.join(args.outdir, f"pspectra_{lap}.bin")
        hdr = read_spectra_header(path)
        S = read_spectra_snapshot(path)
        u = u_bin_centers(hdr)
        dlu = np.diff(np.log(u_bin_edges(hdr)))
        for i, name in SPECIES.items():
            h = S[f"s{i}_u"].sum(axis=(0, 1, 2))
            axs[i].plot(u, h / dlu, color=cmap(k / max(len(laps) - 1, 1)), lw=0.9,
                        label=f"lap {lap}" if k in (0, len(laps) - 1) else None)

    for i, name in SPECIES.items():
        axs[i].set_xscale("log")
        axs[i].set_yscale("log")
        axs[i].set_title(name)
        axs[i].set_xlabel(r"$u$ ($x$ for photons)")
        axs[i].set_ylabel(r"$u\,dW/du$")
        axs[i].legend(fontsize=8)

    fig.tight_layout()
    save = os.path.join(args.outdir, "spectra_time.png")
    fig.savefig(save, dpi=110)
    plt.close(fig)
    print(f"saved {save}")
