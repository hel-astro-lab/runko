# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Energy history of turbulence runs, in units of the guide-field energy density."""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from runko.postprocessing import read_config
from runko.mpiio_reader import read_header, read_field
from runko.mpiio_prtcl_reader import read_prtcl_snapshot
from turb_utils import (laps_per_eddy, numerical_mass, guide_field, snapshot_files,
                        div_forward, div_backward, langevin_antenna, antenna_current,
                        set_plot_style, label_lines, save_figure)


def field_energies(path, B0, U0):
    """Perpendicular, parallel and electric energy densities of one field snapshot.

    The uniform guide field is taken out of B_z, so the parallel entry is the energy
    of the compressive fluctuation alone and starts from zero.
    """
    bperp = np.mean(read_field(path, "bx")**2 + read_field(path, "by")**2)
    bpar = np.mean((read_field(path, "bz") - B0)**2)
    efld = np.mean(read_field(path, "ex")**2 + read_field(path, "ey")**2
                   + read_field(path, "ez")**2)
    return bperp / (2.0 * U0), bpar / (2.0 * U0), efld / (2.0 * U0)


def kinetic_energy(path, conf, species, U0):
    """Kinetic energy density of one species, rebuilt from the sampled particles.

    The sample is drawn with per-tile quotas proportional to the live particle count
    and the particles carry no weights, so the sample mean of gamma - 1 is the mean
    over the whole population. The box is periodic and every cell was filled with the
    same ppc, so the number density is ppc exactly.
    """
    f = read_prtcl_snapshot(path)
    gam = np.sqrt(1.0 + f["ux"]**2 + f["uy"]**2 + f["uz"]**2)
    n = conf.ppc * numerical_mass(conf, species) * conf.cfl**2
    return n * np.mean(gam - 1.0) / U0


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Turbulence energy history")
    parser.add_argument("--dir", required=True, help="simulation output directory")
    parser.add_argument("--linear", action="store_true",
                        help="linear y axis; the default log one spans the whole budget")
    args = parser.parse_args()

    conf = read_config(args.dir)
    teddy = laps_per_eddy(conf)
    B0 = guide_field(conf)
    U0 = 0.5 * B0**2  # guide-field energy density; the unit of every curve below

    times, bperp, bpar, efld = [], [], [], []
    for path in snapshot_files(args.dir, "flds_*.bin"):
        times.append(read_header(path)["lap"] / teddy)
        for column, value in zip((bperp, bpar, efld), field_energies(path, B0, U0)):
            column.append(value)

    kin = np.array([[kinetic_energy(path, conf, s, U0)
                     for path in snapshot_files(args.dir, f"prtcls_{s}_*.bin")]
                    for s in (0, 1)])

    if kin.shape[1] != len(times):
        raise SystemExit("field and particle snapshots were written at different laps")

    times = np.array(times)
    bperp, bpar, efld = np.array(bperp), np.array(bpar), np.array(efld)

    # only the energy the particles have gained is part of the free-energy budget;
    # the rest is the initial thermal content, 2 <gamma - 1> / sigma in these units
    dkin = kin - kin[:, :1]
    total = bperp + bpar + efld + dkin.sum(axis=0)

    print(f"{'driven' if conf.use_antenna else 'decaying'} run, "
          f"{len(times)} snapshots, t/t_0 up to {times[-1]:.2f}")
    print(f"  B_0             = {B0:.5f}   (B_0^2/2 = {U0:.5e})")
    print(f"  K_0(0), K_1(0)  = {kin[0, 0]:.3f}, {kin[1, 0]:.3f}")
    print(f"  budget sum      = {total[0]:.3f} -> {total[-1]:.3f}")
    if not conf.use_antenna:  # nothing injects, so the free energy is conserved at 1
        print(f"  closure         = {total.mean():.3f} +- {total.std():.3f}")

    set_plot_style()

    # single-column width, panel sized to the golden ratio; the margins are absolute
    # so the panel shape does not shift when they are retuned
    golden = 0.5 * (1.0 + np.sqrt(5.0))
    axleft, axright = 0.16, 0.96
    bottom_in, top_in = 0.42, 0.12
    width_in = 3.25
    height_in = (axright - axleft) * width_in / golden + bottom_in + top_in
    axbottom = bottom_in / height_in
    axtop = 1.0 - top_in / height_in

    fig = plt.figure(1, figsize=(width_in, height_in))
    ax = fig.add_subplot(1, 1, 1)
    ax.minorticks_on()

    # the decaying run starts with exactly one unit of free energy, so a closing
    # budget lands on this line; the driven run climbs away from it as the antenna works
    ax.axhline(1.0, color='k', linewidth=0.5, alpha=0.3, zorder=0)

    curves = [(bperp, r"$\delta B_\perp^2/2$", 'C0', 'solid'),
              (bpar, r"$\delta B_\parallel^2/2$", 'C1', 'solid'),
              (efld, r"$E^2/2$", 'C2', 'solid'),
              (kin[0], r"$K_{e^-}$", 'C3', 'solid'),
              (kin[1], r"$K_{e^+}$", 'C4', 'solid'),
              (total, r"$\Sigma$", 'k', 'dashed')]

    ytop = max(kin.max(), total.max())
    ybot = 1.0e-3 * ytop  # log floor; bperp and bpar are exactly zero at t = 0

    # values under the floor are dropped rather than clipped, so the curves do not
    # grow a spurious vertical limb down to the axis
    lines = [ax.plot(times, y if args.linear else np.where(y > ybot, y, np.nan),
                     color=col, linewidth=1, linestyle=ls, label=lab)[0]
             for y, lab, col, ls in curves]

    ax.set_xlim(0.0, times[-1])
    if args.linear:
        ax.set_ylim(0.0, 1.08 * ytop)
    else:
        ax.set_yscale('log')
        ax.set_ylim(ybot, 3.0 * ytop)
    ax.set_xlabel(r"$t / t_0$")
    ax.set_ylabel(r"$U / (B_0^2/2)$")

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    # inline labels instead of a legend; spread them out so they do not collide
    label_lines(lines, times[-1] * np.array([0.55, 0.30, 0.80, 0.42, 0.68, 0.85]),
                rotate=False)

    save_figure(fig, args.dir, "plot_energies")

    if True:  # divergence of B and E; the Yee scheme should hold div B at round-off

        dB, dE = [], []
        for path in snapshot_files(args.dir, "flds_*.bin"):
            b = [read_field(path, n) for n in ("bx", "by", "bz")]
            e = [read_field(path, n) for n in ("ex", "ey", "ez")]
            # each divergence is the one conjugate to the curl that advances its field,
            # so div B is annihilated exactly and div E is the charge density
            dB.append(np.sqrt(np.mean(div_forward(*b)**2)) / B0)
            dE.append(np.sqrt(np.mean(div_backward(*e)**2)) / B0)
        dB, dE = np.array(dB), np.array(dE)

        print(f"  rms div B / B_0 = {dB.max():.2e} (max over laps)")
        print(f"  rms div E / B_0 = {dE.max():.2e} (max over laps; this is the charge)")

        figd = plt.figure(2, figsize=(width_in, height_in))
        axd = figd.add_subplot(1, 1, 1)
        axd.minorticks_on()

        # div B is identically zero while the field is still uniform
        dlines = [axd.plot(times, np.where(dB > 0.0, dB, np.nan), color='C0',
                           linewidth=1, label=r"$\nabla\cdot B$")[0],
                  axd.plot(times, dE, color='C1', linewidth=1,
                           label=r"$\nabla\cdot E$")[0]]

        axd.set_yscale('log')
        axd.set_xlim(0.0, times[-1])
        axd.set_ylim(0.2 * dB[dB > 0.0].min(), 5.0 * dE.max())
        axd.set_xlabel(r"$t / t_0$")
        axd.set_ylabel(r"$\mathrm{rms}(\nabla\cdot X) \,/\, B_0$")
        figd.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)
        label_lines(dlines, times[-1] * np.array([0.5, 0.75]), rotate=False)

        save_figure(figd, args.dir, "plot_energies_div")

    if True:  # antenna injection power, from the antenna replayed out of pic.py

        if not conf.use_antenna:
            print("  no antenna in this run, skipping the injection power")
        else:
            antenna = langevin_antenna(conf)
            paths = snapshot_files(args.dir, "flds_*.bin")
            laps = [read_header(p)["lap"] for p in paths]

            power = []
            for path, lap in zip(paths, laps):
                e = [read_field(path, n) for n in ("ex", "ey", "ez")]
                jant = antenna_current(antenna, lap, e[0].shape, conf.cfl)
                # add_current does E <- E - J, and the snapshot holds E after that
                # subtraction, so the work the antenna did on the field over the lap
                # is -(E + J/2).J with the midpoint field
                power.append(-sum(np.mean((ec + 0.5 * jc) * jc)
                                  for ec, jc in zip(e, jant)) / U0)
            power = np.array(power)

            # trapezoid over the snapshot laps; the antenna decorrelates on ~ 1/gamma_0
            # laps, which the output interval only just resolves
            injected = np.concatenate(
                ([0.0], np.cumsum(0.5 * (power[1:] + power[:-1]) * np.diff(laps))))

            print(f"  antenna power   = {power.min():+.2e} to {power.max():+.2e} "
                  f"per lap in B_0^2/2")
            # E is only on disk every io_output_interval laps, so the trapezoid
            # under-resolves the antenna and recovers ~85 per cent of the budget;
            # the agreement is what confirms the replayed antenna is the right one
            print(f"  injected energy = {injected[-1]:.3f} vs budget {total[-1]:.3f}"
                  f"  ({100.0 * injected[-1] / total[-1]:.0f}%)")

            figa = plt.figure(3, figsize=(width_in, height_in))
            axa = figa.add_subplot(1, 1, 1)
            axa.minorticks_on()
            axa.axhline(0.0, color='k', linewidth=0.5, alpha=0.3, zorder=0)

            alines = [
                axa.plot(times, power * teddy, color='C0', linewidth=1,
                         label=r"$-\langle E\cdot J_\mathrm{ant}\rangle\, t_0$")[0],
                axa.plot(times, injected, color='C3', linewidth=1,
                         label=r"$\int P\,\mathrm{d}t$")[0],
                axa.plot(times, total, color='k', linewidth=1, linestyle='dashed',
                         label=r"$\Sigma$")[0]]

            axa.set_xlim(0.0, times[-1])
            axa.set_xlabel(r"$t / t_0$")
            axa.set_ylabel(r"$U / (B_0^2/2)$")
            figa.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)
            label_lines(alines, times[-1] * np.array([0.35, 0.75, 0.60]), rotate=False)

            save_figure(figa, args.dir, "plot_energies_antenna")

    plt.show()
