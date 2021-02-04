# -*- coding: utf-8 -*-

from __future__ import print_function

import time
import numpy as np


class Timer:

    times = []
    ctimes = []

    do_print = True
    error_index = 99

    verbose = 0

    def __init__(self, names=[], components=[]):
        # self.names = names
        # self.n = len(names)

        # for i in range(self.n):
        #    self.times.append( [] )

        # self.components = components
        # self.nc = len(self.components)

        # for i in range(self.nc):
        #    self.ctimes.append( [] )

        self.names = {}
        for name in names:
            self.names[name] = []

        self.components = {}
        for name in components:
            self.components[name] = []

    # def _look_name(self, name):
    #    #for i, ref_name in enumerate(self.names):
    #    #    if name == ref_name:
    #    #        return i
    #    #return self.error_index
    #    return self.names[name]

    # def _look_comp(self, name):
    #    for i, ref_name in enumerate(self.components):
    #        if name == ref_name:
    #            return i
    #    return self.error_index

    def start(self, name):
        t0 = time.time()  # stop before lookup

        # add name if it does not exist
        # indx = self._look_name(name)
        # if indx == self.error_index:
        #    self.names.append(name)
        #    self.n += 1
        #    self.times.append( [] )
        #    indx = self._look_name(name)

        # self.times[indx] = [] #purge old list
        # self.times[indx].append( t0 )

        # if not(name in self.names):
        #    self.names[name] = []

        # purge old list
        self.names[name] = []
        self.names[name].append(t0)



    def purge(self, name):
        # indx = self._look_name(name)
        # self.times[indx] = [] #purge old list
        # self.times[indx].append( time.time() )

        self.names[name] = []

    def purge_comps(self):
        # for name in self.components:
        #    indx = self._look_comp(name)
        #    self.ctimes[indx] = []
        for name in self.components:
            self.components[name] = []

    def stop(self, name):
        t0 = time.time()  # stop before lookup
        # indx = self._look_name(name)
        # self.times[indx].append( t0 )

        self.names[name].append(t0)

    def lap(self, name):
        t0 = time.time()  # stop before lookup
        # indx = self._look_name(name)
        # self.times[indx].append( t0 )

        if not (name in self.names):
            self.names[name] = []

        self.names[name].append(t0)

    def _calc_mean(self, arr):
        x = np.zeros(len(arr) - 1)
        for ai in range(0, len(arr) - 1):
            x[ai] = arr[ai + 1] - arr[ai]
        return x

    def start_comp(self, name):
        t0 = time.time()  # stop before lookup

        # add comp name if it does not exist
        # indx = self._look_comp(name)
        # if indx == self.error_index:
        #    self.components.append(name)
        #    self.nc += 1
        #    self.ctimes.append( [] )
        #    indx = self._look_comp(name)
        #    print("adding {} into {}".format(name, indx))
        #    print(self.components)
        #    print(self.ctimes)

        # self.ctimes[indx].append( t0 )

        if not (name in self.components):
            self.components[name] = []

        self.components[name].append(t0)

        if self.verbose > 0: print(name) 

        return name

    def stop_comp(self, name):
        t0 = time.time()  # stop before lookup
        # indx = self._look_comp(name)
        # print("stopping {} indx {}".format(name, indx))
        # self.ctimes[indx][-1] =  t0 - self.ctimes[indx][-1]

        self.components[name][-1] = t0 - self.components[name][-1]

    def stats(self, name):

        # print("timer len:", len(self.times), len(self.ctimes) )
        # indx = self._look_name(name)
        # ts = np.array( self.times[indx] )

        ts = np.array(self.names[name])
        cnts = len(ts) - 1
        if cnts < 1:
            return

        t0 = ts[-1] - ts[0]

        if self.do_print:
            if t0 < 1.0e-6:
                print(
                    "--- {:.20}   time: {:8.5f} mus/{:3d}  ({})".format(
                        name.ljust(20), t0 / 1.0e6, cnts, t0
                    )
                )
            elif 1.0e-6 < t0 < 1.0e-3:
                print(
                    "--- {:.20}   time: {:8.5f} ms /{:3d}  ({})".format(
                        name.ljust(20), t0 / 1.0e3, cnts, t0
                    )
                )
            elif 1.0e-3 < t0:
                print(
                    "--- {:.20}   time: {:8.5f} s  /{:3d}  ({})".format(
                        name.ljust(20), t0, cnts, t0
                    )
                )

            # print additional statistics
            if cnts > 1:
                tss = self._calc_mean(ts)
                tavg = np.mean(tss)
                tstd = np.std(tss)
                print(
                    "---            avg: {:8.5f} s   /  {:3d}  ({})".format(
                        tavg, cnts, tavg
                    )
                )
                print(
                    "---            std: {:8.5f} s   /  {:3d}  ({})".format(
                        tstd, cnts, tstd
                    )
                )

    def comp_stats(self):
        if self.do_print:
            print("--------------------------------------------------")

        # calculate total duration
        ts = np.array(self.names["step"])
        cnts = len(ts) - 1
        t0tot = ts[-1] - ts[0]
        if t0tot == 0.0:
            for name in self.components:
                ts = np.array(self.components[name])
                t0tot += np.sum(ts)


        totper = 0.0

        for name in self.components:
            ts = np.array(self.components[name])

            tavg = np.mean(ts)
            cnts = len(ts)
            tstd = np.std(ts)

            t0 = np.sum(ts)
            relt = 100.0 * t0 / t0tot
            totper += relt

            if self.do_print:
                t0 = tavg
                if t0 < 1.0e-5:
                    print(
                        "--- {:.20}   {:6.3f}%  |  time: {:8.5f} mus/{:3d}  ({:.8})".format(
                            name.ljust(20), relt, t0 * 1.0e6, cnts, t0
                        )
                    )
                elif 1.0e-5 < t0 < 1.0e-2:
                    print(
                        "--- {:.20}   {:6.3f}%  |  time: {:8.5f} ms /{:3d}  ({:.8})".format(
                            name.ljust(20), relt, t0 * 1.0e3, cnts, t0
                        )
                    )
                elif 1.0e-2 < t0:
                    print(
                        "--- {:.20}   {:6.3f}%  |  time: {:8.5f} s  /{:3d}  ({:.8})".format(
                            name.ljust(20), relt, t0, cnts, t0
                        )
                    )

                # if cnts > 1:
                #    print("- - -          avg: {:8.5f} s   /  {:3d}  ({:.8})".format(tavg, cnts, tavg))
                #    print("- - -          std: {:8.5f} s   /  {:3d}  ({:.8})".format(tstd, cnts, tstd))

        if self.do_print:
            print("                        += {:6.3f}% ".format(totper))
            print("--------------------------------------------------")

    def dump(self):

        if self.do_print:
            print("--------------------------------------------------")
            for name in self.names:
                self.stats(name)
            print("--------------------------------------------------")
