# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

markers = [
    'o',   # circle
    's',   # square
    '^',   # triangle up
    'v',   # triangle down
    '<',   # triangle left
    '>',   # triangle right
    'D',   # diamond
    'X',   # filled X
    '*',   # star
    '+',   # plus
    'x',   # x
    '1',   # tri down
    '2',   # tri up
    '3',   # tri left
    '4',   # tri right
]

def string_to_color_n_marker(s,
                             special_prefixes=set(["comm_local", "comm_external"]),
                             cmap_name="tab10"):
    """
    Split the string at the first '_' or after entry from special_prefixes
    to prefix and tail part. Give same color to same prefixes
    and same marker for same prefix-tail pair.
    """

    this = string_to_color_n_marker

    try:
        _ = this.first_time
    except:
        this.next_color = 0
        this.next_marker = dict()
        this.first_time = False
        this.seen_prefixes = dict()
        this.seen_tails = dict()

    prefix, tail = None, None
    for p in special_prefixes:
        if s.startswith(p):
            prefix = p
            tail = s[len(p):]
            break

    if not prefix:
        prefix, tail = s.split("_", 1)

    if prefix not in this.seen_prefixes:
        this.seen_prefixes[prefix] = this.next_color
        this.next_color += 1

        this.next_marker[prefix] = 0
        this.seen_tails[prefix] = dict()

    if tail not in this.seen_tails[prefix]:
        this.seen_tails[prefix][tail] = this.next_marker[prefix]
        this.next_marker[prefix] += 1

    color = this.seen_prefixes[prefix]
    marker = this.seen_tails[prefix][tail]

    cmap = plt.get_cmap(cmap_name)
    return cmap(color), markers[marker % len(markers)]


def sorted_by_key(key: list, *vals: list) -> None:
    zipped = list(zip(key, *vals))
    zipped.sort(key=lambda x: x[0])
    return zip(*zipped)



def main(argv: list[str]):
    parser = argparse.ArgumentParser(description="""
        Plot runko weak scaling based on the timing data from given output directories.
        The output directories must be populated with timing data
        by calling `pickle_timer_statistics` method of `runko.Simulation` class.
        Timer statistics contains measurements from N latest laps
        (based on configuration parameter `laps_in_timer_statistics`).
        Some components are not run on every lap,
        so average time contribution per one lap is plotted:
        average_duration * count / laps_in_timer_statistics.
        Unique color-marker-pair is given to each component.
        By default, slowest and fastest ranks are plotted.
        """)

    parser.add_argument("outdir",
                        nargs="+",
                        type=str,
                        help="Paths to runko output directories.")

    parser.add_argument("--push-time",
                        action="store_true",
                        help="""
                        Use push time [s] as units in y-axis.
                        Push time is calculated from duration dt as
                        dt * number_of_ranks / number_of_particles.
                        This require that the simulations have stored
                        the total number of particles in the configuration object
                        as `n_particles` provided to `runko.TileGrid`.""")

    parser.add_argument("--regex-filter",
                        action="store",
                        type=str,
                        help="Show only components which names match given regex.")

    parser.add_argument("--shade-bounds",
                        action="store_true",
                        help="""
                        Shade area between minimum and maximum time for each component.
                        """)

    parser.add_argument("--draw-all",
                        action="store_true",
                        help="Draw all data points.")

    parser.add_argument("--plot-totals",
                        action="store_true",
                        help="Plot sum of slowest components.")

    parser.add_argument("--ymin",
                        action="store",
                        type=float)

    parser.add_argument("--savefig",
                        action="store",
                        type=str,
                        help="Store resulting figure to given path.")

    parser.add_argument("--dpi",
                        action="store",
                        type=int)

    args = parser.parse_args(argv[1:])

    outdirs = args.outdir

    # Load heavy modules only after argument parsing.
    import runko

    fig, (ax_legend, ax) = plt.subplots(2, 1,
                                        layout="constrained",
                                        figsize=(8, 6),
                                        dpi=args.dpi)

    ax.grid(which="both",
            linewidth=1,
            alpha=0.2)

    # (ranks, name) -> (n_particles, list[TimingStatistic])
    data = dict()
    seen_ranks = set()

    for outdir in outdirs:
        print(f"Reading: {outdir}")
        conf = runko.read_config(outdir)
        ranks = conf.ranks
        if ranks in seen_ranks:
            print(f"error: multiple outdirs have same number of ranks ({ranks})")
            quit()
        seen_ranks.add(ranks)
        print(f"ranks: {ranks}")
        all_stats = runko.read_timer_statistics(outdir)

        for i in range(ranks):
            for name, stats in all_stats[i].items():
                if (ranks, name) not in data:
                    data[(ranks, name)] = (conf.n_particles, [stats])
                else:
                    data[(ranks, name)][1].append(stats)

    lines = dict()
    totals = dict()

    regex_filter = None
    if args.regex_filter:
        import re
        regex_filter = re.compile(args.regex_filter)

    for (ranks, name), (n_particles, stats) in data.items():
        if regex_filter:
            if not regex_filter.search(name):
                continue

        l = 1e99
        h = 0

        c, m = string_to_color_n_marker(name)

        for s in stats:
            if args.push_time:
                t = s.average * s.count / s.measured_laps * ranks / n_particles
            else:
                t = s.average * s.count / s.measured_laps

            l = min(l, t)
            h = max(h, t)

            if args.draw_all:
                ax.scatter([ranks],
                           [t],
                           label=name,
                           color=c,
                           marker=m)

        if not args.draw_all:
            ax.scatter([ranks, ranks],
                       [l, h],
                       label=name,
                       color=c,
                       marker=m)

        if ranks not in totals:
            totals[ranks] = h
        else:
            totals[ranks] += h

        if name not in lines:
            lines[name] = ([ranks], [l], [h])
        else:
            lines[name][0].append(ranks)
            lines[name][1].append(l)
            lines[name][2].append(h)


    if len(lines) == 0:
        errmsg ="""Could not find any timing data point to plot.
The outdirs might not contain any timing data"""

        if regex_filter:
            errmsg += " or given regex did not match anything."
        else:
            errmsg += "."

        print(errmsg)
        exit()


    if args.plot_totals:
        x = []
        y = []

        for r, t in totals.items():
            x.append(r)
            y.append(t)

        x, y = sorted_by_key(x, y)

        ax.plot(x, y, '--', label="total", color="black", linewidth=2)


    if args.shade_bounds:
        for name, (r, l, h) in lines.items():
            c, _ = string_to_color_n_marker(name)
            r, l, h = sorted_by_key(r, l, h)
            ax.fill_between(r, l, h, alpha=0.1, color=c)

    # Filter out duplicate entries from the legend.
    handles, labels = ax.get_legend_handles_labels()
    unique_labels_set = set()
    unique_labels = []
    unique_handles = []

    for h, l in zip(handles, labels):
        if l not in unique_labels_set:
            unique_labels_set.add(l)
            unique_labels.append(l)
            unique_handles.append(h)

    # Sort labels based on the name.
    zipped = list(zip(unique_labels, unique_handles))
    zipped.sort(key=lambda x: x[0])
    unique_labels, unique_handles = zip(*zipped)

    ax_legend.legend(unique_handles,
                     unique_labels,
                     mode="expand",
                     loc="lower left",
                     borderaxespad=0,
                     ncol=3)

    ax_legend.set(frame_on=False,
                  xticks=[],
                  yticks=[])

    if args.push_time:
        ylabel = "average push time [s]"
    else:
        ylabel = "average time [s]"

    ax.set(xscale="log",
           yscale="log",
           xticks=list(seen_ranks),
           xlabel="ranks",
           ylabel=ylabel,
           ylim=(args.ymin, None))

    # Required for xticks to work correctly with log scale.
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    if args.savefig:
        fig.savefig(args.savefig)
        print(f"plot stored as: {args.savefig}")
    else:
        plt.show()
