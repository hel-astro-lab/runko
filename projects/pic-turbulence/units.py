from __future__ import print_function
from problem import Configuration_Turbulence as Configuration


from numpy import sqrt, pi
import numpy as np
from scipy.special import kn


# Class to create and store simulation units
class Units:

    def __init__(self, conf):

        #--------------------------------------------------
        # standard things
        self.qe = np.abs(conf.qe)
        self.me = np.abs(conf.me)*self.qe
        self.c  = conf.cfl 
        self.theta = conf.delgam

        #thermal gamma in the beginning
        #self.gam_th = 1.0 + 3.0*self.theta #gamma_thermal
        self.gam_th = 1.55 #manual value for theta = 0.3

        #grid points
        self.Lx = conf.NxMesh*conf.Nx
        self.Ly = conf.NyMesh*conf.Ny
        self.Lz = conf.NzMesh*conf.Nz

        # initial magnetic field energy
        #Benerg = 0.5*Lx*Ly*conf.binit**2

        #forcing scale in units of skin depth
        self.l0 = conf.Nx*conf.NxMesh/conf.max_mode/conf.c_omp  


        # gamma_0 (maximum attainable lorentz factor)
        self.g0 = self.l0*np.sqrt(conf.sigma)*np.sqrt(self.gam_th) #gamma_0
        #self.A = (0.1*3.0/4.0)*self.g0/self.gammarad**2 #definition we use in Runko
        #print("A = {}".format(self.A))


        #--------------------------------------------------
        #emf normalization; OLD; 

        #conf.binit = 1.0
        self.norm_flds = 2.0/conf.binit**2 #mean initial energy density


        # every stride:th value stored 
        self.norm_flds /= (self.Lx/conf.stride)*(self.Ly/conf.stride) #divide by volume of simulation domain
        if self.Lz > 1:
            self.norm_flds /= (self.Lz/conf.stride)


        #--------------------------------------------------
        #particles
        self.n_prtcls = self.Lx*self.Ly*self.Lz*conf.ppc*2.0 #pair-plasma so multiply by 2
        try:
            self.n_test_prtcls = conf.n_test_prtcls
        except:
            print('Exception; no prtcl number found')
            self.n_test_prtcls = 999424.0 # FIXME this is manually deduced number of test particles; 


        #self.n_test_prtcls = conf.n_test_prtcls
        print("total number of particles: ", self.n_prtcls)

        #how many prtcls do test particles represent
        self.norm_prtcls = self.n_prtcls/self.n_test_prtcls 

        print("prtcl/test_prtcl ratio: ", self.norm_prtcls)
        self.norm_prtcls *= self.me*self.c**2 #E = mc^2


        # final re-normalization

        # set energy normalization to initial total magnetic field energy in the box
        # NOTE: problem dependent selection
        self.norm_energy = 0.5*self.Lx*self.Ly*self.Lz*conf.binit**2

        #normalize all to selected energy scale
        self.norm_flds    = self.norm_flds * 1.0 # already in units of total energy in B_0^2/2
        #self.norm_rad     = self.norm_rad / self.norm_energy
        self.norm_prtcls  = self.norm_prtcls / self.norm_energy


    def lap2time(self, lap):
        t0 = 0.5*self.Lz/self.c

        return lap/t0



