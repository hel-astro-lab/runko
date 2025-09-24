
import numpy as np
import pytools  # runko python tools
import sys


from numpy import pi, sqrt, exp
from mpi4py import MPI
from numpy import newaxis as nax


# 3D Langevin antenna that drives a current
class Antenna:
    def __init__(self, min_mode, max_mode, conf):
        # print("initializing...")

        self.cfl = conf.cfl

        self.Nx = conf.Nx
        self.Ny = conf.Ny
        self.Nz = conf.Nz

        self.NxMesh = conf.NxMesh
        self.NyMesh = conf.NyMesh
        self.NzMesh = conf.NzMesh

        # local magnetic field storage
        self.bx = np.zeros((self.NxMesh+2, self.NyMesh+2, self.NzMesh+2))
        self.by = np.zeros((self.NxMesh+2, self.NyMesh+2, self.NzMesh+2))
        self.bz = np.zeros((self.NxMesh+2, self.NyMesh+2, self.NzMesh+2))

        # local current storage
        self.jx = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))
        self.jy = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))
        self.jz = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))

        lx = conf.Nx*self.NxMesh #/ conf.c_omp
        ly = conf.Ny*self.NyMesh #/ conf.c_omp
        lz = conf.Nz*self.NzMesh #/ conf.c_omp

        self.b0 = conf.binit
        self.db = conf.drive_ampl*conf.binit # amplitude of dB perturbations in code units

        # normal indexing
        self.xloc, self.yloc, self.zloc = np.meshgrid(
            np.linspace(-1, self.NxMesh-1+1, self.NxMesh+2),
            np.linspace(-1, self.NyMesh-1+1, self.NyMesh+2),
            np.linspace(-1, self.NzMesh-1+1, self.NzMesh+2),
            indexing='ij',
        )

        #--------------------------------------------------
        # modes that we inject
        self.n_modes = 8
        self.kx = np.zeros(self.n_modes)
        self.ky = np.zeros(self.n_modes)
        self.kz = np.zeros(self.n_modes)

        # standard plane modes
        self.kx[0], self.ky[0], self.kz[0] =  1, 0, 1
        self.kx[1], self.ky[1], self.kz[1] =  1, 0,-1
        self.kx[2], self.ky[2], self.kz[2] =  0, 1, 1
        self.kx[3], self.ky[3], self.kz[3] =  0, 1,-1

        # opposite modes (= -k_m)
        for m in range(4,8):
            n = m - self.n_modes//2 # complementary mode
            self.kx[m], self.ky[m], self.kz[m] = -self.kx[n], -self.ky[n], -self.kz[n]

        # normalization
        self.kx[:] *= 2*pi*max_mode/lx
        self.ky[:] *= 2*pi*max_mode/ly
        self.kz[:] *= 2*pi*max_mode/lz

        #--------------------------------------------------
        self.t0 = lx/conf.cfl/max_mode # light-crossing time across the box (in units of dt)
        self.l0 = lx/max_mode # eddy length 
        valf = np.sqrt(conf.sigma/(1.0 + conf.sigma)) # Alfven velocity in units of c
        w0 = 2.0*pi*conf.cfl*valf/self.l0  # = k_z*v_A*dt

        self.omega = conf.drive_freq*w0  #0.01 # driving frequency
        self.gamma = conf.decorr_time*w0 #0.001 # decorrelation time

        #--------------------------------------------------
        # arrays for normalization of Fourier mode summations
        self.betax = np.zeros(self.n_modes)
        self.betay = np.zeros(self.n_modes)

        for m in range(self.n_modes):
            if np.sqrt(self.kx[m]**2 + self.ky[m]**2) > 0:
                self.betax[m] = +self.ky[m]/sqrt( self.kx[m]**2 + self.ky[m]**2 )
                self.betay[m] = -self.kx[m]/sqrt( self.kx[m]**2 + self.ky[m]**2 )

            else: # TODO experimental branch; not sure this is correct norm
                self.betax[m] = +self.ky[m] #/sqrt( self.n_modes ) #*( self.kx[m]**2 + self.ky[m]**2) )
                self.betay[m] = -self.kx[m] #/sqrt( self.n_modes ) #*( self.kx[m]**2 + self.ky[m]**2) )

        # normalize with amplitude and mode numbers 
        self.betax[:] *= self.db/np.sqrt(self.n_modes)
        self.betay[:] *= self.db/np.sqrt(self.n_modes)


        # initialize random phases for the modes
        self.a_n = np.zeros((self.n_modes), dtype=complex)

        #--------------------------------------------------
        # rank 0 calculates and broadcasts
        twopi_zeta = np.zeros(self.n_modes//2)
        if MPI.COMM_WORLD.Get_rank() == 0:
            #twopi_zeta[:] = 2.0*pi*np.random.rand(self.n_modes//2)

            # manual random seeds
            twopi_zeta    = np.array([5.18222009, 2.85776313, 5.12464465, 0.72490604], dtype='d')
            #twopi_zeta    = np.array([0., 2.85776313, 5.12464465, 0.72490604]) # TODO FIXME

        MPI.COMM_WORLD.Bcast(
                [twopi_zeta, self.n_modes//2, MPI.DOUBLE], 
                root = 0,)

        #--------------------------------------------------
        # initial phases for the driving coefficient
        for m in range(self.n_modes//2):
            self.a_n[m] = exp(1j*twopi_zeta[m] )
            self.a_n[m + self.n_modes//2] = np.conj( self.a_n[m] )
            #print('a_n', m, self.a_n[m], self.a_n[m + self.n_modes//2], 'zeta', twopi_zeta[m])
            #print('other mode', m, m + self.n_modes//2)


        # physical timer inside the antenna to track measured time
        self.tcur = 0.0

        return

    # update internal state of the antenna
    def update_rnd_phases(self):

        ur = np.zeros(self.n_modes, dtype='d')
        ui = np.zeros(self.n_modes, dtype='d')

        # rank0 calculates random number
        if MPI.COMM_WORLD.Get_rank() == 0:
            ur[:] = np.random.rand((self.n_modes))-0.5
            ui[:] = np.random.rand((self.n_modes))-0.5
        
        #--------------------------------------------------
        #MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write

        # update to others
        MPI.COMM_WORLD.Bcast(
                [ur, self.n_modes, MPI.DOUBLE], 
                root = 0,)
        MPI.COMM_WORLD.Bcast(
                [ui, self.n_modes, MPI.DOUBLE], 
                root = 0,)
        #MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write
        #--------------------------------------------------

        for m in range(self.n_modes//2):
            u = ur[m] + 1j*ui[m]

            # a_(n+1) = a_n * exp( -i*(omega + i*gamma)*dt )
            # note that exp(-x) \approx 1 - x which would give finite-difference formula for a_(n+1) - a_n/dt
            self.a_n[m] *= exp(-1j*self.omega - self.gamma)

            # a_n += sigma*u_n *dt
            # where sigma = sqrt(12*|gamma|/dt)
            # this is valid only for uniformly distributed random numbers u_n
            self.a_n[m] += u*sqrt(12.0*self.gamma)

            self.a_n[m + self.n_modes//2] = np.conj(self.a_n[m])

            #print('upd a_n',m, self.a_n[m])
        

        # light crossing time across the largest eddy
        self.tcur += 1.0/self.t0 

        #if MPI.COMM_WORLD.Get_rank() == 0:
        #    print('t ', self.tcur, ' a_0', self.a_n[0]/self.db)

        return


    # calculate next step
    def comp_cur(self, tile):

        self.bx[:,:,:] = 0.0
        self.by[:,:,:] = 0.0

        # normal ij indexing w/ meshgrid
        xx = self.xloc + tile.mins[0]
        yy = self.yloc + tile.mins[1]
        zz = self.zloc + tile.mins[2]

        xxp = xx + 0.5
        yyp = yy + 0.5
        zzp = zz + 0.5
          
        #--------------------------------------------------
        # regularize
        arms = np.sqrt( np.sum( self.a_n*np.conj(self.a_n) ) )
        self.norm_reg = 2*np.sqrt(2)/arms

        if False: # if regularize the power output
            a_n = self.norm_reg*self.a_n[:]
        else:
            a_n = self.a_n


        #--------------------------------------------------
        # fully vectorized version w/ implicit loops over physical dimensions and modes
        #--- add_antenna            56.625%  |  time:  3.83968 s  / 10  (3.8396844)

        #kdotrx =  self.kx[nax, nax, nax, :]*xx[ :,:,:,nax] \
        #        + self.ky[nax, nax, nax, :]*yyp[:,:,:,nax] \
        #        + self.kz[nax, nax, nax, :]*zzp[:,:,:,nax]
        #
        #kdotry =  self.kx[nax, nax, nax, :]*xxp[:,:,:,nax] \
        #        + self.ky[nax, nax, nax, :]*yy[ :,:,:,nax] \
        #        + self.kz[nax, nax, nax, :]*zzp[:,:,:,nax]

        #bxm = self.betax[nax, nax, nax, :]*np.real( 1j*self.a_n[nax, nax, nax, :]*exp(1j*kdotrx[:,:,:,:]))
        #bym = self.betay[nax, nax, nax, :]*np.real( 1j*self.a_n[nax, nax, nax, :]*exp(1j*kdotry[:,:,:,:]))

        #self.bx[:,:,:] = np.sum(bxm, axis=3)
        #self.by[:,:,:] = np.sum(bym, axis=3)


        #--------------------------------------------------
        # semivectorized loop version
        #--- add_antenna            57.802%  |  time:  3.80180 s  / 10  (3.8018026)

        #for m in range(4): # TODO FIXME
        for m in range(self.n_modes):
            kdotrx = self.kx[m]*xx  + self.ky[m]*yyp + self.kz[m]*zzp
            kdotry = self.kx[m]*xxp + self.ky[m]*yy  + self.kz[m]*zzp

            self.bx[:,:,:] += self.betax[m]*np.real( 1j*a_n[m]*exp(1j*kdotrx[:,:,:] ))
            self.by[:,:,:] += self.betay[m]*np.real( 1j*a_n[m]*exp(1j*kdotry[:,:,:] ))

            # check against NaNs
            #is_nanx = np.any(np.isnan(self.bx))
            #is_nany = np.any(np.isnan(self.by))

            #if is_nany or is_nany:
            #    print('nan detected')
            #    print('kx', kdotrx)
            #    print('ky', kdotry)
            #    print('bx', self.bx)
            #    print('by', self.by)
            #    sys.exit()

        #print('difference between loop and vec x:', np.sum(self.bx - bx2))
        #print('difference between loop and vec y:', np.sum(self.by - by2))


        #--------------------------------------------------
        #if False: # slow manual loop version

        #    #--------------------------------------------------
        #    #i0,j0,k0 = pytools.get_index(tile, conf)
        #    i0,j0,k0 = tile.index

        #    #iglob, jglob, kglob = pytools.ind2loc((i0, j0, k0), (0, 0, 0), conf)
        #    xmin = 0.0
        #    ymin = 0.0
        #    zmin = 0.0
        #    iglob = xmin + i0*self.NxMesh 
        #    jglob = ymin + j0*self.NyMesh 
        #    kglob = zmin + k0*self.NzMesh 
        #    #--------------------------------------------------

        #    for k in range(self.NzMesh):
        #        z = kglob + k # global z coord
        #        for j in range(self.NyMesh):
        #            y = jglob + j # global y coord
        #            for i in range(self.NxMesh):
        #                x = iglob + i  # global x coord

        #                for m in range(self.n_modes):
        #                    # dot product of k_m vector with staggered coordinate location
        #                    kdotr_x = self.kx[m]*x        + self.ky[m]*(y+0.5)  + self.kz[m]*(z+0.5)
        #                    kdotr_y = self.kx[m]*(x+0.5)  + self.ky[m]*y        + self.kz[m]*(z+0.5)

        #                    # +1 padding to counter the -1,N+1 limits
        #                    self.bx[i+1,j+1,k+1] += self.betax[m]*np.real( 1j*self.a_n[m]*exp(1j*kdotr_x) )
        #                    self.by[i+1,j+1,k+1] += self.betay[m]*np.real( 1j*self.a_n[m]*exp(1j*kdotr_y) )


        #--------------------------------------------------
        # get current from B field; j = c nabla \times B

        # NOTE this is actually j\Delta t = j*c
        # NOTE no minus sign here since it is included in current deposit routine in runko
        self.jx = +self.cfl*(+self.by[1:-1, 1:-1, 0:-2] - self.by[1:-1, 1:-1, 1:-1] ) # d_z b_y
        self.jy = +self.cfl*(+self.bx[1:-1, 1:-1, 1:-1] - self.bx[1:-1, 1:-1, 0:-2] ) #-d_z b_x
        self.jz = +self.cfl*(+self.bx[1:-1, 0:-2, 1:-1] - self.bx[1:-1, 1:-1, 1:-1]   # d_y b_x 
                             -self.by[0:-2, 1:-1, 1:-1] + self.by[1:-1, 1:-1, 1:-1])  #-d_x b_y

        #print('max j', np.max(np.abs(self.jx)), np.max(np.abs(self.jy)), np.max(np.abs(self.jz)))

        return


    def add_driving(self, tile):
        g = tile.get_grids(0)
        self.comp_cur(tile) # update internal array

        #copy from numpy arrys into tile straight
        for r in range(self.NzMesh):
            for s in range(self.NyMesh):
                for q in range(self.NxMesh):
                    g.bx[q,s,r] += self.bx[q+1,s+1,r+1]
                    g.by[q,s,r] += self.by[q+1,s+1,r+1]
        return


    def add_ext_cur(self, tile):
        g = tile.get_grids(0)
        self.comp_cur(tile) # update internal array

        for r in range(self.NzMesh):
            for s in range(self.NyMesh):
                for q in range(self.NxMesh):
                    g.jx[q,s,r] += self.jx[q,s,r]
                    g.jy[q,s,r] += self.jy[q,s,r]
                    g.jz[q,s,r] += self.jz[q,s,r]


        #max_cur = 0.0
        #for r in range(self.NzMesh):
        #    for s in range(self.NyMesh):
        #        for q in range(self.NxMesh):
        #            jz = np.abs(g.jz[q,s,r])
        #            max_cur = jz if jz > max_cur else max_cur
        #print('max cur', max_cur)

        return


    def get_brms(self, grid):
        vol = self.NxMesh*self.NyMesh*self.NzMesh # tile volume
        volG = self.Nx*self.Ny*self.Nz # tile volume

        # sum power over all local tiles
        b2 = 0.0
        check_x = 0.0
        check_y = 0.0
        check_z = 0.0

        check_x_curl = 0.0
        check_y_curl = 0.0
        check_z_curl = 0.0

        check_x_jdiv = 0.0
        check_y_jdiv = 0.0
        check_z_jdiv = 0.0
        check_jdiv   = 0.0

        check_x_bdiv = 0.0
        check_y_bdiv = 0.0
        check_z_bdiv = 0.0
        check_bdiv   = 0.0

        for tile in pytools.tiles_local(grid):
            self.comp_cur(tile)

            b2 += np.sum(
                      self.bx[1:-1, 1:-1, 1:-1]**2 
                    + self.by[1:-1, 1:-1, 1:-1]**2 )/vol

            # check that curl of J is taken properly and j = curl B
            curlBx = self.jx[:,:,:] - self.cfl*(self.by[1:-1, 1:-1, 0:-2] - self.by[1:-1, 1:-1, 1:-1] 
                                              - self.bz[1:-1, 0:-2, 1:-1] + self.bz[1:-1, 1:-1, 1:-1])
            curlBy = self.jy[:,:,:] - self.cfl*(self.bz[0:-2, 1:-1, 1:-1] - self.bz[1:-1, 1:-1, 1:-1] 
                                              - self.bx[1:-1, 1:-1, 0:-2] + self.bx[1:-1, 1:-1, 1:-1])
            curlBz = self.jz[: :,:] - self.cfl*(self.bx[1:-1, 0:-2, 1:-1] - self.bx[1:-1, 1:-1, 1:-1] 
                                              - self.by[0:-2, 1:-1, 1:-1] + self.by[1:-1, 1:-1, 1:-1])

            check_x_curl = np.max(np.abs(curlBx))
            check_y_curl = np.max(np.abs(curlBy))
            check_z_curl = np.max(np.abs(curlBz))

            # check that divJ = 0
            divJx = -self.jx[1: , :-1, :-1] + self.jx[:-1, :-1, :-1]
            divJy = -self.jy[:-1, 1: , :-1] + self.jy[:-1, :-1, :-1]
            divJz = -self.jz[:-1, :-1,  1:] + self.jz[:-1, :-1, :-1]

            divJ = divJx + divJy + divJz
            check_x_jdiv = np.max(np.abs(divJx))
            check_y_jdiv = np.max(np.abs(divJy))
            check_z_jdiv = np.max(np.abs(divJz))
            check_jdiv   = np.max(np.abs(divJ))


            # NOTE: div not implemented properly for staggered grid
            # check that divB = 0
            divBx = self.bx[1:, :-1, :-1] - self.bx[:-1, :-1, :-1]
            divBy = self.by[:-1, 1:, :-1] - self.by[:-1, :-1, :-1]
            divBz = self.bz[:-1,:-1,  1:] - self.bz[:-1, :-1, :-1]

            divB = divBx + divBy + divBz
            check_x_bdiv = np.max(np.abs(divBx))
            check_y_bdiv = np.max(np.abs(divBy))
            check_z_bdiv = np.max(np.abs(divBz))
            check_bdiv   = np.max(np.abs(divB))

        #--------------------------------------------------
        # reduce accross MPI ranks
        #rbuf = 0.0
        #MPI.COMM_WORLD.Reduce( 
        #        [b2,   MPI.DOUBLE],
        #        [rbuf, MPI.DOUBLE],
        #        op = MPI.SUM,
        #        root = 0,
        #        )
        #b2 = rbuf
        MPI.COMM_WORLD.reduce(b2, op=MPI.SUM, root=0)
        b2 *= 1.0/volG # normalize the tile numbers out

        dbrms = np.sqrt(b2)/self.b0**2

        #--------------------------------------------------
        check_x_curl = MPI.COMM_WORLD.reduce(check_x_curl, op=MPI.MAX, root=0)
        check_y_curl = MPI.COMM_WORLD.reduce(check_y_curl, op=MPI.MAX, root=0)
        check_z_curl = MPI.COMM_WORLD.reduce(check_z_curl, op=MPI.MAX, root=0)

        check_x_div = MPI.COMM_WORLD.reduce(check_x_jdiv, op=MPI.MAX, root=0)
        check_y_div = MPI.COMM_WORLD.reduce(check_y_jdiv, op=MPI.MAX, root=0)
        check_z_div = MPI.COMM_WORLD.reduce(check_z_jdiv, op=MPI.MAX, root=0)
        check_div   = MPI.COMM_WORLD.reduce(check_jdiv,   op=MPI.MAX, root=0)

        check_x_bdiv = MPI.COMM_WORLD.reduce(check_x_bdiv, op=MPI.MAX, root=0)
        check_y_bdiv = MPI.COMM_WORLD.reduce(check_y_bdiv, op=MPI.MAX, root=0)
        check_z_bdiv = MPI.COMM_WORLD.reduce(check_z_bdiv, op=MPI.MAX, root=0)
        check_bdiv   = MPI.COMM_WORLD.reduce(check_bdiv,   op=MPI.MAX, root=0)

        #--------------------------------------------------
        if MPI.COMM_WORLD.Get_rank() == 0:
            print(' antenna amplitude: dB_rms/B_0', dbrms)
            print(' curlB check', check_x_curl, check_y_curl, check_z_curl)
            print('  divJ check', check_jdiv, check_x_jdiv,  check_y_jdiv,  check_z_jdiv)
            print('  divB check', check_bdiv, check_x_bdiv,  check_y_bdiv,  check_z_bdiv)

        return


