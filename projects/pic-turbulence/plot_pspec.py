# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Particle momentum spectra of turbulence runs, histogrammed from the sampled particles."""

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from runko.postprocessing import read_config
from runko.mpiio_prtcl_reader import read_prtcl_header, read_prtcl_snapshot
from turb_utils import laps_per_eddy, snapshot_files, set_plot_style, save_figure


def maxwell_juttner(u, theta):
    """Isotropic Maxwell-Juttner dN/dln(u), unnormalized; theta = kT/mc^2."""
    return u**3 * np.exp(-np.sqrt(1.0 + u**2) / theta)


def log_broken_powerlaw(u, lgA, alpha1, alpha2, delta, lg_ubreak):
    """log10 of a smoothly broken power law f ~ u^-alpha1 below the break, u^-alpha2
    above it, delta setting the width of the transition. Built from logaddexp so the
    high-u limb cannot overflow for small delta."""
    ln10 = np.log(10.0)
    lnr = np.log(u) - lg_ubreak * ln10
    smooth = (np.logaddexp(0.0, lnr / delta) - np.log(2.0)) / ln10
    return lgA - alpha1 * lnr / ln10 + (alpha1 - alpha2) * delta * smooth


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Turbulence particle spectra")
    parser.add_argument("--dir", required=True, help="simulation output directory")
    parser.add_argument("--species", type=int, default=0, help="particle species to plot")
    parser.add_argument("--nbins", type=int, default=48, help="number of momentum bins")
    args = parser.parse_args()

    conf = read_config(args.dir)
    teddy = laps_per_eddy(conf)

    # four-velocity magnitude u = |beta gamma| of every sampled particle, lap by lap
    times, moms = [], []
    for path in snapshot_files(args.dir, f"prtcls_{args.species}_*.bin"):
        f = read_prtcl_snapshot(path)
        moms.append(np.sqrt(f["ux"]**2 + f["uy"]**2 + f["uz"]**2))
        times.append(read_prtcl_header(path)["lap"] / teddy)

    times = np.array(times)
    print(f"{len(times)} snapshots, {len(moms[0])} particles each, "
          f"t/t_0 up to {times[-1]:.2f}")

    # bins shared by all laps so the curves are directly comparable. Counts are divided
    # by the sample size and the log bin width and scaled to the initial particles per
    # cell, so every spectrum integrates to ppc no matter how many particles were
    # sampled, which pins the peak near unity and lets the axis limits stay fixed.
    umin = min(m.min() for m in moms)
    umax = max(m.max() for m in moms)
    edges = np.logspace(np.log10(umin), np.log10(umax), args.nbins + 1)
    u = np.sqrt(edges[:-1] * edges[1:])
    dlnu = np.log(edges[1] / edges[0])
    spectra = [conf.ppc * np.histogram(m, bins=edges)[0] / (len(m) * dlnu)
               for m in moms]

    set_plot_style()

    # single-column width, panel sized to the golden ratio; the margins are absolute
    # so the panel shape does not shift when they are retuned
    golden = 0.5 * (1.0 + np.sqrt(5.0))
    axleft, axright = 0.20, 0.96
    bottom_in, top_in = 0.42, 0.52
    width_in = 3.25
    height_in = (axright - axleft) * width_in / golden + bottom_in + top_in
    axbottom = bottom_in / height_in
    axtop = 1.0 - top_in / height_in

    fig = plt.figure(1, figsize=(width_in, height_in))
    ax = fig.add_subplot(1, 1, 1)
    ax.minorticks_on()

    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=times[-1])
    cmap = matplotlib.colormaps['turbo']

    for t, spec in zip(times, spectra):
        ax.plot(u, np.where(spec > 0, spec, np.nan),
                color=cmap(norm(t)), linewidth=0.8)

    # initial thermal distribution, peak-matched to the first snapshot
    theta = conf.theta0 * (conf.theta1_to_theta0 if args.species else 1.0)
    mj = maxwell_juttner(u, theta)
    mj *= spectra[0].max() / mj.max()
    ax.plot(u, mj, color='k', linewidth=1, linestyle='dashed')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.1, 100.0)
    ax.set_ylim(1.0e-2, 1.0e1)
    ax.set_xlabel(r"$\beta\gamma$")
    ax.set_ylabel(r"$\mathrm{d}N/\mathrm{d}\ln \gamma\beta$")

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    cax = fig.add_axes([axleft, axtop + 0.03, axright - axleft, (axtop - axbottom) * 0.04])
    cb = matplotlib.colorbar.ColorbarBase(
        cax, cmap=cmap, norm=norm, orientation='horizontal', ticklocation='top')
    cb.set_label(r"$t / t_0$")

    save_figure(fig, args.dir, f"plot_pspec_s{args.species}")

    if True:  # smoothly broken power law fitted to the final spectrum
        spec = spectra[-1]
        counts = np.histogram(moms[-1], bins=edges)[0]
        ok = counts > 0
        sigma = 1.0 / (np.log(10.0) * np.sqrt(counts[ok]))  # Poisson error in log10
        guess = [np.log10(spec.max()), -3.0, 3.0, 0.3, np.log10(u[spec.argmax()])]
        popt, pcov = curve_fit(
            log_broken_powerlaw, u[ok], np.log10(spec[ok]), p0=guess, sigma=sigma,
            bounds=([-5.0, -10.0, -10.0, 0.02, -2.0], [5.0, 10.0, 10.0, 3.0, 3.0]),
            maxfev=40000)
        perr = np.sqrt(np.diag(pcov))
        ubreak = 10.0**popt[4]
        resid = (np.log10(spec[ok]) - log_broken_powerlaw(u[ok], *popt)) / sigma

        print(f"broken power law fit at t/t_0 = {times[-1]:.2f}")
        print(f"  alpha_1     = {popt[1]:+.3f} +- {perr[1]:.3f}")
        print(f"  alpha_2     = {popt[2]:+.3f} +- {perr[2]:.3f}")
        print(f"  Delta       = {popt[3]:.3f} +- {perr[3]:.3f}")
        print(f"  gamma_break = {np.sqrt(1.0 + ubreak**2):.3f}"
              f"   (u_break = {ubreak:.3f} +- {perr[4] * ubreak * np.log(10.0):.3f})")
        print(f"  chi2/dof    = {np.sum(resid**2) / (ok.sum() - len(popt)):.2f}")

    plt.show()
