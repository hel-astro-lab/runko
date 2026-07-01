# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import time
import numpy as np
from dataclasses import dataclass
from .runko_logging import runko_logger


@dataclass
class TimeMeasurement:
    begin: float
    end: float | None = None


class Timer:
    def __init__(self):
        self.__times = dict()
        self.__logger = runko_logger("Timer")


    def __uniq_name(self, name: str):
        if name not in self.__times:
            return name

        n = 1
        while True:
            name_candidate = f"{name}-{n}"
            if name_candidate not in self.__times:
                return name_candidate
            n += 1


    def start(self, name: str):
        """
        Starts timer with given name.
        If the name is already present, use `name-n` as the name,
        where n is integer incremented for each used name.
        """

        uniq_name = self.__uniq_name(name)
        self.__logger.debug(uniq_name)
        self.__times[uniq_name] = TimeMeasurement(begin=time.time())


    def stop(self, name: str):
        """
        Stops the timer with given name.
        If it is already stopped, then look for names `name-n` in order
        for n = 1, 2, ... until first unstopped timer is found which is then stopped
        or if `name-n` is not started. In the latter case exception is raised.
        """

        stop_time = time.time()

        if name in self.__times and not self.__times[name].end:
            self.__times[name].end = stop_time
            return

        n = 1
        while True:
            name_candidate = f"{name}-{n}"
            if name_candidate in self.__times:
                if self.__times[name_candidate].end is None:
                    self.__times[name_candidate].end = stop_time
                    return
                else:
                    n += 1
                    continue
            else:
                raise RuntimeError("Trying to stop time which has not been started.")


    def get_elapsed_times(self):
        """
        Returns elapsed times in seconds.
        """

        return {name: t.end - t.begin for name, t in self.__times.items()}


    @property
    def time_measurements(self):
        """
        Dictionary that maps strings to runko.TimeMeasurement.
        """
        return self.__times


@dataclass
class TimerStatistics:
    def __init__(self, data, measured_laps: int):
        self.data = np.array(data)
        self.measured_laps = measured_laps

    @property
    def total(self):
        return np.sum(self.data)

    @property
    def average(self):
        return np.mean(self.data)

    @property
    def minimum(self):
        return np.min(self.data)

    @property
    def maximum(self):
        return np.max(self.data)

    @property
    def std_dev(self):
        return np.std(self.data)

    @property
    def count(self):
        return len(self.data)


def timer_statistics(timers: list[Timer]):
    all_elapsed_times = dict()

    for timer in timers:
        elapsed_times = timer.get_elapsed_times()

        for name, dt in elapsed_times.items():
            if name in all_elapsed_times:
                all_elapsed_times[name].append(dt)
            else:
                all_elapsed_times[name] = [dt]

    stats = dict()
    for name, elapsed_times in all_elapsed_times.items():
        stats[name] = TimerStatistics(data=elapsed_times, measured_laps=len(timers))

    return stats
