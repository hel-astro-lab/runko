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
    ret_val = []
    for x in zip(*zipped):
        ret_val.append(list(x))
    return ret_val


def calculate_max_splits(data: list):
    """
    Takes list of runko.TimerStatistics which each represent
    data from different rank. Assumes that each TimerStatistics.count
    is the same.
    """

    if len(data) == 0:
        return np.array([])

    laps = data[0].count
    mins = 1e99 * np.ones(laps, dtype=float)
    maxs = np.zeros(laps, dtype=float)

    for lap_times in data:
        assert(laps == lap_times.count)

        for lap, time in enumerate(lap_times.data):
            mins[lap] = min(mins[lap], time)
            maxs[lap] = max(maxs[lap], time)

    return np.array(maxs) - np.array(mins)


def main(argv: list[str]):
    parser = argparse.ArgumentParser(
        prog="runko weak-scaling",
        description="""
        Plot runko weak scaling based on the timing data from given output directories.
        The output directories must be populated with timing data
        by calling `pickle_timer_statistics` method of `runko.Simulation` class.
        Timer statistics contains measurements from N latest laps
        (based on configuration parameter `io_n_laps_in_timer_stats`).
        Some components are not run on every lap,
        so average time contribution per one lap is plotted:
        average_duration * count / io_n_laps_in_timer_stats.
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

    parser.add_argument("--max-splits",
                        action="store_true",
                        help="""
                        Plots the maximum splits for each component.
                        Maximum split is defined to be the difference between the slowest
                        and the fastest rank on the same lap.
                        All available laps are searched for the largest maximum split.
                        Incompatible with: --draw-all, --shade-bounds and --plot-totals.
                        """)

    parser.add_argument("--stackplot",
                        action="store_true",
                        help="""
                        Plots the slowest times accross ranks as a stackplot.
                        Incompatible with: --draw-all, --shade-bounds and --plot-totals.
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

    # Maps args.* to what argument --* it actually represents.
    name_map = lambda x: "--" + x.replace("_", "-")

    arg_conflicts = {"max_splits": ["draw_all", "shade_bounds", "plot_totals"],
                     "stackplot": ["draw_all", "shade_bounds", "plot_totals"]}

    # This assumes that arg_confligs only contains action=store_true args.
    arg_is_given = lambda x: hasattr(args, x) and getattr(args, x)
    for arg, conflicts in arg_conflicts.items():
        if arg_is_given(arg):
            for x in conflicts:
                if arg_is_given(x):
                    print(f"error: {name_map(x)} does not make sense with {name_map(arg)}")
                    exit(1)

    outdirs = args.outdir

    # Load heavy modules only after argument parsing.
    import runko

    disable_legend = args.stackplot

    if disable_legend:
        fig, (ax) = plt.subplots(1, 1,
                                 layout="constrained",
                                 figsize=(1.6 * 6, 6),
                                 dpi=args.dpi)
    else:
        fig, (ax_legend, ax) = plt.subplots(2, 1,
                                            layout="constrained",
                                            figsize=(8, 6),
                                            dpi=args.dpi)

    grid_alpha = 0.2
    if args.stackplot:
        grid_alpha = 0.4
    ax.grid(which="both",
            linewidth=1,
            alpha=grid_alpha)

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

        if args.max_splits:
            max_split = np.max(calculate_max_splits(stats))

            if args.push_time:
                max_split *= ranks / n_particles

            if not args.stackplot:
                ax.scatter([ranks],
                           [max_split],
                           label=name,
                           color=c,
                           marker=m)

            # max_splits uses h.
            h = max_split
        else:
            for s in stats:
                if args.push_time:
                    t = s.average * s.count / s.measured_laps * ranks / n_particles
                else:
                    t = s.average * s.count / s.measured_laps

                l = min(l, t)
                h = max(h, t)

                if args.draw_all and not args.stackplot:
                    ax.scatter([ranks],
                               [t],
                               label=name,
                               color=c,
                               marker=m)

            if not args.draw_all and not args.stackplot:
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

    if args.max_splits and not args.stackplot:
        for name, (r, _, h) in lines.items():
            c, _ = string_to_color_n_marker(name)
            r, h = sorted_by_key(r, h)
            ax.plot(r, h, alpha=0.2, color=c)

    if args.stackplot:
        shared_r = None
        names = []
        data = []
        colors = []

        for name, (r, l, h) in lines.items():
            colors.append(string_to_color_n_marker(name)[0])

            r, h = sorted_by_key(r, h)
            if shared_r and shared_r != r:
                raise RuntimeError("Internal logic error. Please send bug report :)")
            else:
                shared_r = r

            names.append(name)
            data.append(h)

        if not shared_r or len(shared_r) <= 1:
            raise RuntimeError("stackplot error: there is not enough data to plot")

        key = [-x[-1] for x in data]
        _, data, colors, names = sorted_by_key(key, data, colors, names)
        ax.stackplot(shared_r, data,
                     colors=colors,
                     edgecolor="black",
                     alpha=0.6)

        # Write labels manually:
        height_rhs = sum([x[-1] for x in data])
        next_anno_y = data[0][-1] * 0.5
        anno_x = shared_r[-1]

        for n, name in enumerate(names):

            # if data[n][-1] / height_rhs < 0.05:

            number_of_lines = 0
            if not n < number_of_lines:
                arrow = None
            else:
                total_arm = 20
                armB = 2 + int(total_arm * n / number_of_lines)
                armA = 2 + total_arm - armB
                arrow = dict(arrowstyle="-",
                             connectionstyle=f"arc,angleA=-180,angleB=0,armA={armA},armB={armB},rad=0",
                             relpos=(0, 0.5),
                             color=colors[n])

            ax.annotate(name,
                        xy=(anno_x, next_anno_y),
                        xycoords="data",
                        xytext=(1.02, n / len(names)),
                        textcoords=ax.transAxes,
                        color=colors[n],
                        arrowprops=arrow)

            if n + 1 < len(names):
                next_anno_y += 0.5 * (data[n][-1] + data[n + 1][-1])

    def setup_legend():
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

        unique_labels, unique_handles = sorted_by_key(unique_labels, unique_handles)
        ax_legend.legend(unique_handles,
                         unique_labels,
                         mode="expand",
                         loc="lower left",
                         borderaxespad=0,
                         ncol=3)

        ax_legend.set(frame_on=False, xticks=[], yticks=[])

    if not disable_legend:
        setup_legend()

    if args.max_splits:
        if args.push_time:
            ylabel = "max split (push time) [s]"
        else:
            ylabel = "max split [s]"
    else:
        if args.push_time:
            ylabel = "average push time [s]"
        else:
            ylabel = "average time [s]"

    if args.stackplot:
        ax.set(xscale="log",
               yscale="linear",
               xlabel="ranks",
               ylabel=ylabel)
    else:
        ax.set(xscale="log",
               yscale="log",
               xlabel="ranks",
               ylabel=ylabel,
               ylim=(args.ymin, None))


    ax.get_xaxis().set_major_locator(mpl.ticker.FixedLocator(list(seen_ranks)))
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_locator(mpl.ticker.NullLocator())

    if args.savefig:
        fig.savefig(args.savefig)
        print(f"plot stored as: {args.savefig}")
    else:
        plt.show()
