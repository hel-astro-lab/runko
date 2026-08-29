# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Mid-plane slices of a qed-pic field snapshot: the fields, the current and the number densities
of the electrons, the positrons and the photons; bz is shown about its mean.

    python plot_qed_fields.py [--outdir corona-pic] [--lap N] [--axis z]
"""

import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from runko.mpiio_reader import read_header, read_field_snapshot

PANELS = (("dbz", "RdBu_r", True), ("bx", "RdBu_r", True), ("by", "RdBu_r", True),
          ("ex", "RdBu_r", True), ("ey", "RdBu_r", True), ("jz", "RdBu_r", True),
          ("n0", "viridis", False), ("n1", "viridis", False), ("n2", "inferno", False))


def field_laps(outdir):
    return sorted(int(os.path.basename(f)[5:-4]) for f in glob.glob(os.path.join(outdir, "flds_*.bin")))


def plot(outdir, lap, axis, save):
    path = os.path.join(outdir, f"flds_{lap}.bin")
    hdr = read_header(path); F = read_field_snapshot(path)
    b0 = np.mean(F["bz"])
    F["dbz"] = F["bz"] - b0                     # the perturbation about the mean field
    sl = {"z": lambda a: a[a.shape[0] // 2, :, :], "y": lambda a: a[:, a.shape[1] // 2, :], "x": lambda a: a[:, :, a.shape[2] // 2]}[axis]
    fig, axs = plt.subplots(3, 3, figsize=(13, 12))
    for ax, (name, cmap, signed) in zip(axs.flat, PANELS):
        a = sl(F[name])
        if signed:
            v = np.max(np.abs(a)) or 1.0
            im = ax.imshow(a, origin="lower", cmap=cmap, vmin=-v, vmax=v)
            ax.set_title(f"{name}: max |{name}|/B0 = {v / b0:.3g}" if b0 else name)
        else:
            im = ax.imshow(a, origin="lower", cmap=cmap)
            ax.set_title(name)
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(f"lap {lap}, slice through the middle along {axis}")
    fig.tight_layout(); fig.savefig(save, dpi=110); plt.close(fig)
    print(f"saved {save}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="field slices of a qed-pic run")
    parser.add_argument("--outdir", default="corona-pic")
    parser.add_argument("--lap", type=int, default=None, help="snapshot lap (default: the last)")
    parser.add_argument("--axis", default="z", choices=("x", "y", "z"))
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    laps = field_laps(args.outdir)
    todo = laps if args.all else [args.lap if args.lap is not None else laps[-1]]
    for lap in todo:
        plot(args.outdir, lap, args.axis, os.path.join(args.outdir, f"flds_{lap}_{args.axis}.png"))
