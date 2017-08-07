import time
import numpy as np

class Timer:

    times = []

    def __init__(self, names):
        self.names = names
        self.n = len(names)

        for i in range(self.n):
            self.times.append( [] )

    def _look_name(self, name):
        for i, ref_name in enumerate(self.names):
            if name == ref_name:
                return i
        return 99

    def start(self, name):
        indx = self._look_name(name)
        self.times[indx] = [] #purge old list
        self.times[indx].append( time.time() )

    def stop(self, name):
        t0 = time.time() #stop before lookup
        indx = self._look_name(name)
        self.times[indx].append( t0 )

    def lap(self, name):
        t0 = time.time() #stop before lookup
        indx = self._look_name(name)
        self.times[indx].append( t0 )

    def _calc_mean(self, arr):
        x = np.zeros(len(arr)-1)
        for ai in range(0, len(arr)-1):
            x[ai] = arr[ai+1] - arr[ai]
        return x

    def stats(self, name):

        indx = self._look_name(name)
        ts = np.array( self.times[indx] )
        cnts = len(ts) - 1
        if cnts < 1:
            return

        t0 = ts[-1] - ts[0]

        if t0 < 1.0e-6:
            print "--- {:.6}   time: {:8.5f} mus /  {:3d}  ({})".format(name, t0/1.0e6,  cnts, t0)
        elif 1.0e-6 < t0 < 1.0e-3:
            print "--- {:.6}   time: {:8.5f} ms  /  {:3d}  ({})".format(name, t0/1.0e3, cnts, t0)
        elif 1.0e-3 < t0:
            print "--- {:.6}   time: {:8.5f} s   /  {:3d}  ({})".format(name, t0,       cnts, t0)

        #print additional statistics
        if cnts > 1:
            tss = self._calc_mean(ts)
            tavg = np.mean(tss)
            tstd = np.std(ts)
            print "---          avg: {:8.5f} s   /  {:3d}  ({})".format(tavg, cnts, tavg)
            print "---          std: {:8.5f} s   /  {:3d}  ({})".format(tstd, cnts, tstd)


    def dump(self):

        print "--------------------------------------------------"
        for name in self.names:
            self.stats(name)
        print "--------------------------------------------------"






