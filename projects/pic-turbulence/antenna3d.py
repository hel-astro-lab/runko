import numpy as np
import pyrunko.tools as pytools
import sys

from numpy import pi

# Magnetic antenna that perturbs B field with sinusoidal modes
class Antenna:
    def __init__(self, mode_perp, mode_par, conf):
        # print("initializing...")

        self.NxMesh = conf.NxMesh
        self.NyMesh = conf.NyMesh
        self.NzMesh = conf.NzMesh

        # temporary working arrays
        self.bxm = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))
        self.bym = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))

        # domain sizes
        self.Lx = conf.Nx * self.NxMesh / conf.c_omp
        self.Ly = conf.Ny * self.NyMesh / conf.c_omp
        self.Lz = conf.Nz * self.NzMesh / conf.c_omp
        self.L = self.Lx  # assuming L = Lx=Ly=Lz box
        # print("L = ({},{},{})".format(self.Lx, self.Ly, self.Lz))

        # absolute driving amplitude
        self.A0 = conf.binit * conf.drive_ampl

        # manually set antenna modes
        self.n_perp = mode_perp
        self.n_par = mode_par

        # initialize seed (so that this is deterministic)
        seed = int(mode_perp)
        np.random.seed(seed)

        # phi_mn and \varphi_mn
        self.ph1 = 2.0 * np.pi * np.random.rand(self.n_perp, self.n_perp, self.n_par)
        self.ph2 = 2.0 * np.pi * np.random.rand(self.n_perp, self.n_perp, self.n_par)
        self.ph3 = 2.0 * np.pi * np.random.rand(self.n_perp, self.n_perp, self.n_par)

        # location matrix
        # NOTE: numpy vs runko array conventions
        self.xloc, self.yloc, self.zloc = np.meshgrid(
            np.linspace(0.0, self.NyMesh-1, self.NyMesh),
            np.linspace(0.0, self.NxMesh-1, self.NxMesh),
            np.linspace(0.0, self.NzMesh-1, self.NzMesh),
        )

        self.kx = 2.0 * np.pi / (conf.c_omp * self.Lx)
        self.ky = 2.0 * np.pi / (conf.c_omp * self.Ly)
        self.kz = 2.0 * np.pi / (conf.c_omp * self.Lz)


    # calculate dB_x, and dB_y
    def get_dB(self, tile):

        # length of the wave mode vector
        def beta(n, m):
            return np.sqrt(8.0)/np.sqrt(n*n + m*m)/self.n_perp/np.sqrt(self.n_par)


        t_mins = np.array(tile.mins)
        t_maxs = np.array(tile.maxs)

        mins = np.zeros(3)
        maxs = np.zeros(3)
        for i in range(0, len(t_mins)):
            mins[i] = 1.0 * t_mins[i]
            maxs[i] = 1.0 * t_maxs[i]

        # NOTE: numpy vs runko array conventions
        xx = self.xloc + mins[1]
        yy = self.yloc + mins[0]
        zz = self.zloc + mins[2]

        # zero initialize arrays before vector operations
        self.bxm[:,:,:] = 0.0
        self.bym[:,:,:] = 0.0

        # loop over modes
        for n in range(1,self.n_perp+1):
            for m in range(1,self.n_perp+1):

                norm = beta(n, m)
                for o in range(1,self.n_par+1):

                    xmodx = np.sin(m*self.kx*xx + self.ph1[n-1,m-1,o-1]) 
                    xmody = np.cos(n*self.ky*yy + self.ph2[n-1,m-1,o-1])

                    ymodx = np.cos(m*self.kx*xx + self.ph1[n-1,m-1,o-1]) 
                    ymody = np.sin(n*self.ky*yy + self.ph2[n-1,m-1,o-1])

                    zmod  = np.sin(o*self.kz*zz +  self.ph3[n-1,m-1,o-1])
                    #zmod = 1.0

                    self.bxm += norm*n*xmodx*xmody*zmod
                    self.bym -= norm*m*ymodx*ymody*zmod

        # normalize from 0,1 to amplitude
        self.bxm *= self.A0
        self.bym *= self.A0


    def add_driving(self, tile):
        gs = tile.get_grids(0)

        self.get_dB(tile)

        #copy from numpy arrys into tile straight
        for r in range(self.NzMesh):
            for s in range(self.NyMesh):
                for q in range(self.NxMesh):
                    gs.bx[q,s,r] = self.bxm[q,s,r]
                    gs.by[q,s,r] = self.bym[q,s,r]


