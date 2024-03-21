from __future__ import print_function

import pytools
from pytools import Configuration

import numpy as np
from numpy import sqrt, pi


# Prtcl velocity (and location modulation inside cell)
#
# NOTE: Cell extents are from xloc to xloc + 1, i.e., dx = 1 in all dims.
#       Therefore, typically we use position as x0 + RUnif[0,1).
#
# NOTE: Default injector changes odd ispcs's loc to (ispcs-1)'s prtcl loc.
#       This means that positrons/ions are on top of electrons to guarantee
#       charge conservation (because zero charge density initially).
#
def velocity_profile(xloc, ispcs, conf):

    # electrons
    if ispcs == 0:
        delgam = conf.delgam_e

    # positrons/ions/second species
    elif ispcs == 1:
        delgam = conf.delgam_i

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    zz = xloc[2] + np.random.rand(1)

    # velocity sampling
    gamma = 0.0
    direction = +1
    ux, uy, uz, uu = pytools.sample_boosted_maxwellian(
        delgam, gamma, direction=direction, dims=3
    )

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0


# number of prtcls of species 'ispcs' to be added to cell at location 'xloc'
#
# NOTE: Plasma frequency is adjusted based on conf.ppc (and prtcl charge conf.qe/qi
#       and mass conf.me/mi are computed accordingly) so modulate density in respect
#       to conf.ppc only.
#
def density_profile(xloc, ispcs, conf):

    # uniform plasma with default n_0 number density
    return conf.ppc


# Field initialization
def insert_em_fields(grid, conf):

    for tile in pytools.tiles_all(grid):
        yee = tile.get_grids(0)

        ii,jj,kk = tile.index if conf.threeD else (*tile.index, 0)

        # insert values into Yee lattices; includes halos from -3 to n+3
        for n in range(-3, conf.NzMesh + 3):
            for m in range(-3, conf.NyMesh + 3):
                for l in range(-3, conf.NxMesh + 3):

                    # get global coordinates
                    #iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                    #r = np.sqrt(iglob ** 2 + jglob ** 2 + kglob ** 2)

                    yee.bx[l, m, n] = 0.0
                    yee.by[l, m, n] = 0.0
                    yee.bz[l, m, n] = conf.binit 

                    yee.ex[l, m, n] = 0.0
                    yee.ey[l, m, n] = 0.0
                    yee.ez[l, m, n] = 0.0
    return






# extend default conf class with problem specific parameters
class Configuration_Turbulence(Configuration):

    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        # problem specific initializations
        if do_print:
            print("Initializing turbulence setup...")

        # local variables just for easier/cleaner syntax
        me = np.abs(self.me)
        mi = np.abs(self.mi)
        c = self.cfl
        ppc = self.ppc * 2.0  # multiply x2 to account for 2 species/pair plasma

        # plasma reaction & subsequent normalization
        self.gamma = 1.0
        self.omp = c / self.c_omp
        self.qe = -(self.omp ** 2.0 * self.gamma) / ((ppc * 0.5) * (1.0 + me / mi))
        self.qi = -self.qe

        me *= abs(self.qi)
        mi *= abs(self.qi)

        # temperatures
        self.delgam_e = self.delgam
        self.delgam_i = self.delgam


        # ---------cold plasma-----------
        # parse external magnetic field strength from sigma_ext

        # determine initial magnetic field based on magnetization sigma which
        # is magnetic energy density/ kinetic energy density
        # this definition works even for nonrelativistic flows.
        #self.binit = sqrt(
        #    (self.gamma) * ppc * 0.5 * c ** 2.0 * (me * (1.0 + me / mi)) * self.sigma
        #)

        # no corrections; cold sigma
        self.binit_nc = sqrt(ppc*c**2.*self.sigma*me)

        # approximative \gamma_th = 1 + 3\theta
        self.gammath = 1.0 + 3.0*self.delgam_e
        self.binit_approx = sqrt(self.gammath*ppc*me*c**2.*self.sigma)

        #manual value for theta=0.3
        #self.gammath = 1.55 

        self.gammath = 1.0 # cool plasma
        self.binit = sqrt(self.gammath*ppc*me*c**2.*self.sigma)

        self.sigma_perp = self.sigma*self.drive_ampl**2
        self.binit_perp = sqrt(self.gammath*ppc*me*c**2.*self.sigma_perp)

        if do_print:
            print("init: sigma: ", self.sigma)
            print("init: mass term: ", sqrt(mi+me))
            print("")
            print("init: B_guide (manual): ", self.binit)
            print("init: B_guide (no corrections): ", self.binit_nc)
            print("init: B_guide (approx): ", self.binit_approx)
            print("")


        #forcing scale in units of skin depth
        self.l0 = self.Nx*self.NxMesh/self.max_mode/self.c_omp  
        
        #thermal larmor radius
        lth = self.gammath / np.sqrt(self.sigma_perp)*np.sqrt(self.gammath) #gamma_0

        #reconnection radius; gam ~ sigma
        lsig = self.sigma_perp / np.sqrt(self.sigma_perp)*np.sqrt(self.gammath) #gamma_0

        self.g0 = self.l0*np.sqrt(self.sigma_perp)*np.sqrt(self.gammath) #gamma_0

        if do_print:
            print("init: l_0:", self.l0)
            print("init: l_th:", lth)
            print("init: l_sig:", lsig)
            print("init: gamma_0: ", self.g0)
            print("")

        #running time estimates
        lx = self.Nx*self.NxMesh/self.max_mode
        t0 = lx/self.cfl/self.drive_ampl**2

        if do_print:
            print("init: lap(t = 5  l_0/c):",  5*t0)
            print("init: lap(t = 10 l_0/c):", 10*t0)
            print("init: lap(t = 20 l_0/c):", 20*t0)
            print("init: sampling rate:",  self.interval/t0)


        # NOTE manually switch off rank memory optimization mode
        self.mpi_task_mode = False

