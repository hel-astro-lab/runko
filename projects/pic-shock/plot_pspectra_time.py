"""Plot particle spectra (e- and e+) as a function of time."""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import os

import runko
from runko.auto_outdir import resolve_outdir
from runko.mpiio_spectra_reader import (
    read_spectra_header,
    read_spectra_snapshot,
    u_bin_centers,
)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Particle spectra time evolution")
    parser.add_argument("--conf", required=True, help="Path to .ini config file")
    args = parser.parse_args()

    # -- config ----------------------------------------------------------------
    conf = runko.Configuration(args.conf)
    outdir = resolve_outdir(conf)
    cfl = conf.cfl
    c_omp = conf.c_omp

    print(f"outdir = {outdir}")

    # -- load spectra snapshots ------------------------------------------------
    times = []
    spectra_s0 = []
    spectra_s1 = []
    u = None

    for lap in range(0, conf.Nt + 1, conf.output_interval):
        fpath = os.path.join(outdir, f"pspectra_{lap}.bin")
        if not os.path.exists(fpath):
            continue

        hdr = read_spectra_header(fpath)
        if u is None:
            u = u_bin_centers(hdr)

        fields = read_spectra_snapshot(fpath)
        # sum over spatial axes (z, y, x) to get total f(u)
        spectra_s0.append(fields["s0_u"].sum(axis=(0, 1, 2)))
        spectra_s1.append(fields["s1_u"].sum(axis=(0, 1, 2)))
        times.append(lap * cfl / c_omp)

    times = np.array(times)
    print(f"loaded {len(times)} spectra snapshots")

    # -- figure setup ----------------------------------------------------------
    fig = plt.figure(1, figsize=(3.25, 4.0))

    plt.rc('font', family='serif')
    plt.rc('text', usetex=False)
    plt.rc('xtick', top=True, direction='out', labelsize=7)
    plt.rc('ytick', right=True, direction='out', labelsize=7)
    plt.rc('axes', labelsize=8)
    plt.rc('legend', handlelength=4.0)

    nrow_fig = 2
    ncol_fig = 1
    gs = plt.GridSpec(nrow_fig, ncol_fig)
    gs.update(wspace=0.25, hspace=0.12)

    axs = np.empty((nrow_fig, ncol_fig), dtype=object)
    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i, j] = plt.subplot(gs[i, j])
            axs[i, j].minorticks_on()

    axleft = 0.18; axbottom = 0.10; axright = 0.96; axtop = 0.82
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    ax_em = axs[0, 0]  # e- (species 0)
    ax_ep = axs[1, 0]  # e+ (species 1)

    # -- colormap --------------------------------------------------------------
    tmax = conf.Nt * cfl / c_omp
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=tmax)
    cmap = matplotlib.colormaps['turbo']

    # -- plot spectra ----------------------------------------------------------
    for i, t in enumerate(times):
        col = cmap(norm(t))
        ax_em.plot(u, spectra_s0[i], color=col, linewidth=0.8)
        ax_ep.plot(u, spectra_s1[i], color=col, linewidth=0.8)

    # -- axes ------------------------------------------------------------------
    for ax in [ax_em, ax_ep]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(u[0], u[-1])
        ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(
            base=10.0, subs=np.arange(2, 10), numticks=100))

    # ylim: tight to data with floor
    all_s0 = np.array(spectra_s0)
    all_s1 = np.array(spectra_s1)
    ymin = 0.5
    ymax_s0 = all_s0.max() * 3.0
    ymax_s1 = all_s1.max() * 3.0
    ax_em.set_ylim(ymin, ymax_s0)
    ax_ep.set_ylim(ymin, ymax_s1)

    ax_em.set_xticklabels([])
    ax_ep.set_xlabel(r"$\beta \gamma$")
    ax_em.set_ylabel(r"$u\,\mathrm{d}N/\mathrm{d}u$")
    ax_ep.set_ylabel(r"$u\,\mathrm{d}N/\mathrm{d}u$")

    # panel labels
    ax_em.text(0.95, 0.92, r"$e^-$", transform=ax_em.transAxes,
               fontsize=8, ha='right', va='top')
    ax_ep.text(0.95, 0.92, r"$e^+$", transform=ax_ep.transAxes,
               fontsize=8, ha='right', va='top')

    # -- colorbar --------------------------------------------------------------
    axwidth = axright - axleft
    axheight = (axtop - axbottom) * 0.03
    axpad = 0.02
    cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])

    cb = matplotlib.colorbar.ColorbarBase(
        cax, cmap=cmap, norm=norm,
        orientation='horizontal', ticklocation='top')
    cb.set_label(r"Time $t~(\omega_p^{-1})$")

    # -- save ------------------------------------------------------------------
    base = os.path.join(outdir, "plot_pspectra_time")
    plt.savefig(base + ".pdf")
    #plt.savefig(base + ".png", dpi=300)
    print(f"saved {base}.pdf")
