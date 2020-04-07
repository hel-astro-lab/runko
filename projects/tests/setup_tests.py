from __future__ import print_function

from pytools import Configuration

import numpy as np
from numpy import sqrt, pi



def lap2time(lap, conf, do_print=False):
    # from relativistic pair plasma oscillation freq to classical single species freq
    omp = conf.omp*np.sqrt(conf.gamma/2)
    return lap*omp


# extend default conf class with problem specific parameters
class Configuration_Test(Configuration):
    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        # problem specific initializations
        if do_print:
            print("Initializing problem setup...")

        # local variables just for easier/cleaner syntax
        me = np.abs(self.me)
        mi = np.abs(self.mi)
        c = self.cfl
        ppc = self.ppc * 2.0  # multiply x2 to account for 2 species/pair plasma

        # if gammas < 1, we interpret it as beta/c
        if self.gamma == 0:
            self.gamma = 1.0
            self.beta = 0.0
        else:
            if self.gamma < 1.0:
                self.gamma = sqrt(1.0 / (1.0 - self.gamma ** 2.0))
            self.beta = sqrt(1.0 - 1.0 / self.gamma ** 2.0)

        # plasma reaction & subsequent normalization
        self.omp = c / self.c_omp
        self.qe = -(self.omp ** 2.0 * self.gamma) / ((ppc * 0.5) * (1.0 + me / mi))
        self.qi = -self.qe

        # take given charge into account
        self.qe *= -np.sign(self.me)
        self.qi *=  np.sign(self.mi)

        me *= abs(self.qi)
        mi *= abs(self.qi)

        # temperatures
        self.delgam_e = self.delgam
        self.delgam_i = self.temp_ratio * self.delgam_e

        # ---------cold plasma-----------
        # parse external magnetic field strength from sigma_ext

        # determine initial magnetic field based on magnetization sigma which
        # is magnetic energy density/ kinetic energy density
        # this definition works even for nonrelativistic flows.
        self.binit = sqrt(
            (self.gamma) * ppc * 0.5 * c ** 2.0 * (me * (1.0 + me / mi)) * self.sigma
        )



