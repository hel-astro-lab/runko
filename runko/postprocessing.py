# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import pickle
import pathlib
from .configuration import Configuration
from .runko_timer import TimerStatistic


def read_config(outdir: str) -> Configuration:
    with open(f"{outdir}/config.pkl", "rb") as f:
         conf = pickle.load(f)

    return conf


def read_timer_statistics(outdir: str) -> dict:
    """
    Returns dict that maps rank integer to dict read from outdir
    that maps component name to runko.TimerStatistic.
    """

    stats = dict()
    path = pathlib.Path(f"{outdir}/timer-statistics")

    for f in path.iterdir():
        rank = int(f.stem)
        with open(f, "rb") as f:
            stats[rank] = pickle.load(f)

    return stats
