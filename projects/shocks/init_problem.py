from __future__ import print_function
from configSetup import Configuration

from numpy import sqrt, pi
import numpy as np


class Configuration_Shocks(Configuration):

    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        #-------------------------------------------------- 
        # problem specific initializations
        if do_print:
            print("Initializing shock setup...")
    
        #-------------------------------------------------- 
        # particle initialization

        #local variables just for easier/cleaner syntax
        me  = np.abs(self.me)
        mi  = np.abs(self.mi)
        c   = self.cfl
        ppc = self.ppc*2.0 #multiply x2 to account for 2 species/pair plasma

        # if gammas < 1, we interpret it as beta/c
        if(self.gamma < 1.):
            self.gamma = sqrt(1./(1.-self.gamma**2.)) 
        self.beta = sqrt(1.-1./self.gamma**2.)
        
	#plasma reaction & subsequent normalization
        omp=c/self.c_omp
        self.qe = -(omp**2.*self.gamma)/((ppc*.5)*(1.+me/mi)) 
        self.qi = -self.qe 

        me *= abs(self.qi)
        mi *= abs(self.qi)

        #temperatures
        self.delgam_e = self.delgam
        self.delgam_i = self.temp_ratio*self.delgam_e
        
        #--------------------------------------------------
        # printing 

        if do_print:
            print("init: Positron thermal spread: ", self.delgam_i)
            print("init: Electron thermal spread: ", self.delgam_e)

        sigmaeff = self.sigma #* temperature corrections
        if do_print:
            print("init: Alfven (outflow) three-velocity: ",sqrt(sigmaeff/(1.+sigmaeff)))
            print("init: Ion beta: ",     2.*self.delgam_i/(self.sigma*(mi/me+1.)/(mi/me)) )
            print("init: Electron beta: ",2.*self.delgam_e/(self.sigma*(mi/me+1.)) )

        #-------------------------------------------------- 
        # field initialization

        #---------cold plasma-----------
        # parse external magnetic field strength from sigma_ext
	#self.bz_ext = sqrt( (self.gamma-1.0)*.5*ppc*c**2*(mi+me)*self.sigma_ext)

	#determine initial magnetic field based on magnetization sigma which 
        #is magnetic energy density/ kinetic energy density
	#this definition works even for nonrelativistic flows. 
	#self.binit = sqrt((self.gamma)*ppc*.5*c**2.*(mi+me)*self.sigma)


        #----------hot plasma----------
        corrdelgam_qe = 1.0
        corrdelgam_sig = 1.0

        delgam_i = self.delgam_i
        delgam_e = self.delgam_e

        zeta=delgam_i/(0.24 + delgam_i)
        gad_i=1./3.*(5 - 1.21937*zeta + 0.18203*zeta**2 - 0.96583*zeta**3 + 2.32513*zeta**4 - 2.39332*zeta**5 + 1.07136*zeta**6)
        delgam_e=self.delgam*mi/me*self.temp_ratio 
        zeta=delgam_e/(0.24 + delgam_e)
        gad_e=1./3.*(5 - 1.21937*zeta + 0.18203*zeta**2 - 0.96583*zeta**3 + 2.32513*zeta**4 - 2.39332*zeta**5 + 1.07136*zeta**6)

        self.binit=1.*sqrt(ppc*.5*c**2.* \
                (mi*(1.+corrdelgam_sig*gad_i/(gad_i-1.)*self.delgam_i)+me*(1.+ \
                corrdelgam_sig*gad_e/(gad_e-1.)*self.delgam_e))*self.sigma)

        if do_print:
            print("init: sigma: ", self.sigma)
            print("init: B_guide: ", self.binit)




        #-------------------------------------------------- 
        # radiation drag, if any

        if "gammarad" in self.__dict__:
            if not(self.gammarad == 0.0):
                self.drag_amplitude = 0.1*self.binit/(self.gammarad**2.0)

                if do_print:
                    print("using radiation drag...")
                    print(" drag amplitude: {} with gamma_rad: {}".format(self.drag_amplitude, self.gammarad))
        else:
            self.gammarad = 0.0

        if "radtemp" in self.__dict__:
            if not(self.radtemp == 0.0):
                if do_print:
                    print("using radiation drag with temperature...")
                    print(" drag temperature: {}".format(self.radtemp))
        else:
            self.radtemp = 0.0

