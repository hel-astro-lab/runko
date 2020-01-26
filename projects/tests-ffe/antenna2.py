
import numpy as np
import pyrunko.tools as pytools
import sys

from numpy import pi 


# 
#
class Antenna:


    def __init__(self, min_mode, max_mode, conf):
        #print("initializing...")

        self.NxMesh = conf.NxMesh
        self.NyMesh = conf.NyMesh
        self.NzMesh = conf.NzMesh

        self.Bx_rms = pytools.Mesh_H3(self.NxMesh, self.NyMesh, self.NzMesh)
        self.By_rms = pytools.Mesh_H3(self.NxMesh, self.NyMesh, self.NzMesh)

        self.Lx = conf.Nx*self.NxMesh/conf.c_omp
        self.Ly = conf.Ny*self.NyMesh/conf.c_omp
        self.Lz = conf.Nz*self.NzMesh/conf.c_omp
        self.L  = self.Lx                        #assuming L = Lx=Ly=Lz box
        #print("L = ({},{},{})".format(self.Lx, self.Ly, self.Lz))

        self.min_mode = min_mode
        self.max_mode = max_mode

        # m,n \in {1,...,N}
        self.ms = np.arange(1,max_mode+1)
        self.ns = np.arange(1,max_mode+1)

        # absolute driving amplitude
        self.A0 = conf.binit*conf.drive_ampl


        # power in different modes
        self.beta = np.zeros((max_mode, max_mode))
        # TODO:
        # modify max_mode = Nmax-Nmin
        # and change limits to n = Nmin:Nmax etc.
        for n in range(1,max_mode+1):
            if n < min_mode:
                continue
            for m in range(1,max_mode+1):
                if m < min_mode:
                    continue
                self.beta[m-1,n-1] = 2.0*self.A0/np.sqrt(n*n + m*m)/max_mode


        #initialize seed (so that this is deterministic)
        seed = np.int(max_mode)
        np.random.seed(seed)

        # phi_mn and \varphi_mn
        self.phases1 = 2.0*np.pi*np.random.rand(max_mode, max_mode)
        self.phases2 = 2.0*np.pi*np.random.rand(max_mode, max_mode)

        # location matrix
        self.xloc, self.yloc, self.zloc = np.meshgrid(
                np.linspace(0.0, self.NxMesh-1.0, self.NxMesh),
                np.linspace(0.0, self.NyMesh-1.0, self.NyMesh),
                np.linspace(0.0, self.NzMesh-1.0, self.NzMesh))

        self.kx = 2.0*np.pi / (conf.c_omp*self.Lx)
        self.ky = 2.0*np.pi / (conf.c_omp*self.Ly)
        self.kz = 2.0*np.pi / (conf.c_omp*self.Lz)





    # step driver coefficient according to eq. 14
    def step_driver(self, lap, dt):
        #print("stepping...")

        np.random.seed(np.int(lap))

        #delta-correlated uniform random complex number 
        # between Re(un) in [-1/2, 1/2] and Im(un) in [-1/2, 1/2]
        #un = (np.random.rand() - 0.5) + 1.0j*(np.random.rand()-0.5)

        # finite-difference stepping (remembering the dt that was missing from sigma)
        #self.an += dt*(-1.0j*self.wa*self.an + self.sigma*un/np.sqrt(dt) )

        # TODO: need to step phases here


    # calculate B_{rms,x}, and B_{rms,y}
    def get_Brms(self, tile):
        t_mins = np.array( tile.mins ) 
        t_maxs = np.array( tile.maxs )

        mins = np.zeros(3)
        maxs = np.zeros(3)
        for i in range(0,len(t_mins)):
            mins[i] = 1.0*t_mins[i]
            maxs[i] = 1.0*t_maxs[i]

        #expand axis to vectorize this
        xx = self.xloc[:, :, :, np.newaxis, np.newaxis] + mins[0]
        yy = self.yloc[:, :, :, np.newaxis, np.newaxis] + mins[1]

        Bx =  self.beta*self.ns  \
                       *np.sin(self.kx*self.ms*xx + self.phases1) \
                       *np.cos(self.ky*self.ns*yy + self.phases2)

        #By = -self.beta*self.ms \
        By = -(self.beta*self.ms).T \
                       *np.cos(self.kx*self.ms*xx + self.phases1) \
                       *np.sin(self.ky*self.ns*yy + self.phases2)

        Bx = np.sum(Bx, axis=(3,4))
        By = np.sum(By, axis=(3,4))

        #print("Bx")
        #print(self.Bx)
        #print(np.shape(self.xloc))
        #print(np.shape(self.Bx))

        #print("By")
        #print(self.By)
        #print(np.shape(self.xloc))
        #print(np.shape(self.By))

        # NOTE: flipping of x vs y indices
        for r in range(self.NzMesh):
            for s in range(self.NyMesh):
                for q in range(self.NxMesh):
                    #self.Jext[q,s,r] = kdotr[s,q,r]
                    self.Bx_rms[q,s,r] = Bx[s,q,r]
                    self.By_rms[q,s,r] = By[s,q,r]


    def get_Brms2(self, tile):
        t_mins = np.array( tile.mins ) 
        t_maxs = np.array( tile.maxs )

        mins = np.zeros(3)
        maxs = np.zeros(3)
        for i in range(0,len(t_mins)):
            mins[i] = 1.0*t_mins[i]
            maxs[i] = 1.0*t_maxs[i]

        #xx = self.xloc[:, :, :] + mins[0]
        #yy = self.yloc[:, :, :] + mins[1]

        #Bx =  self.beta*self.ns  \
        #               *np.sin(self.kx*self.ms*xx + self.phases1) \
        #               *np.cos(self.ky*self.ns*yy + self.phases2)

        #By = -(self.beta*self.ms).T \
        ##By = -self.beta*self.ms \
        #               *np.cos(self.kx*self.ms*xx + self.phases1) \
        #               *np.sin(self.ky*self.ns*yy + self.phases2)

        xloc = np.linspace(0.0, self.NxMesh-1.0, self.NxMesh) + mins[0]
        yloc = np.linspace(0.0, self.NyMesh-1.0, self.NyMesh) + mins[1]
        zloc = np.linspace(0.0, self.NzMesh-1.0, self.NzMesh) + mins[2]

        max_mode = self.max_mode
        for r in range(self.NzMesh):
            for s in range(self.NyMesh):
                for q in range(self.NxMesh):

                    self.Bx_rms[q,s,r] = 0.0
                    self.By_rms[q,s,r] = 0.0

                    for n in range(1,max_mode+1):
                        for m in range(1,max_mode+1):
                            beta = 2.0*self.A0/np.sqrt(n*n + m*m)/max_mode

                            Bx =  beta*n \
                                    *np.sin(self.kx*m*xloc[q] + self.phases1[n-1,m-1])  \
                                    *np.cos(self.ky*n*yloc[s] + self.phases2[n-1,m-1])

                            By = -beta*m \
                                    *np.cos(self.kx*m*xloc[q] + self.phases1[n-1,m-1]) \
                                    *np.sin(self.ky*n*yloc[s] + self.phases2[n-1,m-1])

                            self.Bx_rms[q,s,r] += Bx
                            self.By_rms[q,s,r] += By



    def add_driving(self, tile):
        yee = tile.get_yee(0)
        #self.get_Brms(tile)
        self.get_Brms2(tile)

        yee.bx += self.Bx_rms
        yee.by += self.By_rms
        

        #yee.jz += self.Jext

        #for mode in self.modes:
        #    self.get_current(mode, tile)
        #    yee.jz += self.Jext







