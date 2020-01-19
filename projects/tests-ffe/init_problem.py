from __future__ import print_function
from configSetup import Configuration

from numpy import sqrt, pi
import numpy as np


class Configuration_Reconnection(Configuration):

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
        self.delgam_i = self.delgam_e
        
        #--------------------------------------------------
        # problem related
        #FIXME
        #self.dstripe  = 1.0/(self.dstripe*self.c_omp)
        #self.dvstripe = 1.0/(self.dvstripe*self.c_omp)

        self.sheet_thickness = self.sheet_thickness*self.c_omp #skind into cells
        self.pinch_width = self.pinch_width*self.c_omp         #skind into cells


        #box center
        self.mx = self.Nx*self.NxMesh
        self.my = self.Ny*self.NyMesh
        self.mz = self.Nz*self.NzMesh

        #middle of the box
        if not(self.periodicx):
            self.mxhalf  = int(0.5*self.mx)
            self.lstripe = 1.0 #FIXME delete
        else:
            self.mxhalf  = int(0.25*self.mx)
            self.lstripe = int(self.mx) #FIXME delete
        self.myhalf = int(0.5*self.my)
        self.mzhalf = int(0.5*self.mz)


        #--------------------------------------------------
        # printing 

        if do_print:
            sigmaeff = self.sigma #* temperature corrections
            print("init: Alfven (outflow) three-velocity: ",sqrt(sigmaeff/(1.+sigmaeff)))
            print("init: sheet thickness ",    self.sheet_thickness )
            print("init: sheet density ", self.sheet_density )
            print("init: trigger B", self.trigger )
            print("init: pinch width ", self.pinch_width )
            print("init: trigger Ex", self.trigger_field )

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


        #FIXME
        if False:

            #hot plasma version
            corrdelgam_qe = 1.0
            corrdelgam_sig = 1.0

            delgam_i = self.delgam_i
            delgam_e = self.delgam_e

            zeta=delgam_i/(0.24 + delgam_i)
            gad_i=1./3.*(5 - 1.21937*zeta + 0.18203*zeta**2 - 0.96583*zeta**3 + 2.32513*zeta**4 - 2.39332*zeta**5 + 1.07136*zeta**6)
            delgam_e=self.delgam*mi/me*self.temperature_ratio 
            zeta=delgam_e/(0.24 + delgam_e)
            gad_e=1./3.*(5 - 1.21937*zeta + 0.18203*zeta**2 - 0.96583*zeta**3 + 2.32513*zeta**4 - 2.39332*zeta**5 + 1.07136*zeta**6)

            self.binit=1.*sqrt(ppc*.5*c**2.* \
                    (mi*(1.+corrdelgam_sig*gad_i/(gad_i-1.)*self.delgam_i)+me*(1.+ \
                    corrdelgam_sig*gad_e/(gad_e-1.)*self.delgam_e))*self.sigma)

            self.bphi = self.bphi/180.*pi

            if do_print:
                print("sin(bphi:", np.sin(self.bphi))

