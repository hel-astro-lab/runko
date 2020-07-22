from __future__ import print_function
from pytools import Configuration
import pytools

from numpy import sqrt, pi
import numpy as np


class Configuration_Packets(Configuration):

    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        #-------------------------------------------------- 
        # problem specific initializations
        if do_print:
            print("Initializing reconnection setup...")
    
        #-------------------------------------------------- 
        # particle initialization

        #local variables just for easier/cleaner syntax
        me  = np.abs(self.me)
        mi  = np.abs(self.mi)
        c   = self.cfl
        ppc = self.ppc*2.0 #multiply x2 to account for 2 species/pair plasma

        # if gammas < 1, we interpret them as beta/c
        self.gamma = 1.0
        #if(self.gamma < 1.):
        #    self.gamma = sqrt(1./(1.-self.gamma**2.)) 
        #self.beta = sqrt(1.-1./self.gamma**2.)
        
        #plasma reaction & subsequent normalization
        omp=c/self.c_omp
        self.qe = -(omp**2.*self.gamma)/((ppc*.5)*(1.+me/mi))
        self.qi = -self.qe

        me *= abs(self.qi)
        mi *= abs(self.qi)

        #temperatures
        self.delgam_e = self.delgam
        self.delgam_i = self.delgam_e
        

        #-------------------------------------------------- 
        # field initialization

        # parse external magnetic field strength from sigma_ext
        #if self.external_fields:
        #    self.bz_ext = sqrt(
        #        (self.gamma-1.0)*.5*ppc*c**2*(mi+me)*self.sigma_ext)

        #determine initial magnetic field based on magnetization sigma which 
        #is magnetic energy density/ kinetic energy density
        #this definition works even for nonrelativistic flows. 
        self.binit = sqrt((self.gamma)*ppc*.5*c**2.*(mi+me)*self.sigma)

        # Alfven speed
        #NOTE: FFE Alfven speed is 1
        #self.beta = np.sqrt(self.sigma/(self.sigma + 1.0))
        self.beta = 1.0


        #-------------------------------------------------- 
        # special setups for colliding packets
        if True:
            Lx = self.Nx*self.NxMesh
            Ly = self.Ny*self.NyMesh
            Lz = self.Nz*self.NzMesh

            #b0 = self.binit   # equilibrium magnetic field
            #zeta = self.zeta  # perturbation amplitude

            # maximum perpendicular and parallel wavenumbers
            #kperp = 1.0*pi/Lx #half of the sine only
            #kpar  = 2.0*pi/Lz

            #phase shift
            #om0 = 0.0 

            #amplitude; TODO
            #self.A = b0*zeta/2.0

            #-------------------------------------------------- 
            # packet centers
            
            # middle of the box
            x0 = Lx*0.5
            y0 = Ly*0.5

            # distort packet centers
            x1 = x0 - 0.5*self.impact_param
            x2 = x0 + 0.5*self.impact_param

            # z center of both pkgs
            z1 = Lz*0.25
            z2 = Lz*0.75

            #print("x0 {} {} {}".format(x0, y0, z1))

            # position of the centers as Stagger objects
            self.pkg_loc1 = pytools.Stagger(x1, y0, z1)
            self.pkg_loc2 = pytools.Stagger(x2, y0, z2)


        if do_print:
            print("init: sigma:             ", self.sigma)
            print("init: mass term:         ", sqrt(mi+me))
            print("init: B_guide (manual) : ", self.binit)
            print("init: beta (Alfven vel): ", self.beta)


        #forcing scale in units of skin depth
        self.l0 = self.Nx*self.NxMesh/self.c_omp  
        self.l1 = self.ell/self.c_omp  

        #thermal larmor radius
        self.gammath = 1.0 # cool plasma
        lth = self.gammath / np.sqrt(self.sigma)*np.sqrt(self.gammath) #gamma_0

        #reconnection radius; gam ~ sigma
        lsig = self.sigma / np.sqrt(self.sigma)*np.sqrt(self.gammath) #gamma_0
        self.g0 = self.l0*np.sqrt(self.sigma)*np.sqrt(self.gammath) #gamma_0

        if do_print:
            print("init: ell:", self.ell, 'dx ', self.ell/self.c_omp, 'comp ', self.ell/(self.Nx*self.NxMesh), '%')
            print("init: l_0:", self.l0)
            print("init: l_1:", self.l1)
            print("init: l_th:", lth)
            print("init: l_sig:", lsig)
            print("init: gamma_0: ", self.g0)
            #print("init: A = {}".format(A))

