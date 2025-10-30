"""
usage: <this-script> <filepath> [shift-x [shift-y]]

Plots the output csv file from pic-beam-instabilities/beam.py.

filepath: path to the csv file.
shift-x:  how much growth rate lines are shifted in x-direction.
shift-x:  how much growth rate lines are shifted in y-direction.
"""


import numpy as np
import sys
import argparse

import matplotlib
import matplotlib.pyplot as plt


def setup_figure(columns=1,
                 nrow_fig=1,
                 ncol_fig=1,
                 # control these (in units of [0,1]) to position the figure
                 axleft=0.18,
                 axbottom=0.16,
                 axright=0.96,
                 axtop=0.92,
                 figsize_scale=1,
):
    """
    Based on: https://natj.github.io/teaching/figs/
    """

    if columns == 1:
        fig = plt.figure(1, figsize=(figsize_scale * 3.25, figsize_scale * 2.2))
    elif columns == 2:
        fig = plt.figure(1, figsize=(figsize_scale * 7.0, figsize_scale * 2.5))
    else:
        raise RuntimeError("This only supports single- or two-column figures")

    # add ticks to both sides
    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    # make labels slightly smaller
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=8)
    plt.rc('legend', handlelength=4.0)

    gs = plt.GridSpec(nrow_fig, ncol_fig)
    gs.update(wspace = 0.25)
    gs.update(hspace = 0.35)

    axs = np.empty((nrow_fig, ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    return fig, axs


if __name__ == "__main__":

    data = np.loadtxt(sys.argv[1], delimiter=",")

    t = data[:, 0]
    B2_parallel = data[:, 1]
    B2_perpendicular = data[:, 2]

    E2_parallel = data[:, 3]
    E2_perpendicular = data[:, 4]

    two_stream_growth = data[:, 5]
    oblique_growth = data[:, 6]
    filamentation_growth = data[:, 7]

    fig, axs = setup_figure(figsize_scale=1)

    axs[0, 0].plot(t, B2_parallel, label="$B_T$")
    axs[0, 0].plot(t, B2_perpendicular, label="$B_F$")
    axs[0, 0].plot(t, E2_parallel, label="$E_T$")
    axs[0, 0].plot(t, E2_perpendicular, label="$E_F$")

    axs[0, 0].plot(t, E2_parallel + E2_perpendicular, label="$E$")
    axs[0, 0].plot(t, B2_parallel + B2_perpendicular, label="$B$")

    ymin = 1e-5

    plot_shift_x = 0
    if len(sys.argv) >= 3:
        plot_shift_x = float(sys.argv[2])
    plot_shift_y = ymin
    if len(sys.argv) >= 4:
        plot_shift_y = float(sys.argv[3])


    axs[0, 0].plot(t + plot_shift_x, plot_shift_y * two_stream_growth, '--', label="two-stream")
    axs[0, 0].plot(t + plot_shift_x, plot_shift_y * oblique_growth, '--', label="oblique")
    axs[0, 0].plot(t + plot_shift_x, plot_shift_y * filamentation_growth, '--', label="filamentation")

    axs[0,0].set(yscale="log",
                 ylim=(ymin, 1))

    axs[0, 0].legend()

    plt.show()
