import time
import numpy as np
from dataclasses import dataclass


@dataclass
class TimeMeasurement:
    begin: float
    end: float | None = None


class Timer:
    def __init__(self):
        self.__times = dict()


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


@dataclass
class TimerStatistic:
    total: float
    average: float
    minimum: float
    maximum: float
    std_dev: float
    count: int


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
        et = np.array(elapsed_times)

        stats[name] = TimerStatistic(total=np.sum(et),
                                     average=np.mean(et),
                                     minimum=np.min(et),
                                     maximum=np.max(et),
                                     std_dev=np.std(et),
                                     count=len(elapsed_times))

    return stats
