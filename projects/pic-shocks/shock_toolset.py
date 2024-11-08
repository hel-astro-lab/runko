from __future__ import print_function
import numpy as np
import h5py as h5

from pytools import ind2loc
import pytools

from mpi4py import MPI


def MPI_OP_overlap(xmem, ymem, dt):
    x = np.frombuffer(xmem, dtype=np.float32)
    y = np.frombuffer(ymem, dtype=np.float32)

    # zero out x values that would overlap with y
    xnonz = np.where(x != 0.0)
    y[xnonz] = 0.0

    z = x + y
    y[:] = z

op_overlap = MPI.Op.Create(MPI_OP_overlap, commute=True)



def ifloor(x):
    return int(np.floor(x))

class ShockToolset:

    def __init__(self, conf):

        # global density profile along x
        self.density_rank = np.zeros((conf.Nx*conf.NxMesh), dtype='float32')
        self.density      = np.zeros((conf.Nx*conf.NxMesh), dtype='float32')

        self.nx = conf.Nx*conf.NxMesh
        self.ny = conf.Ny*conf.NyMesh
        self.nz = conf.Nz*conf.NzMesh

        self.outdir = conf.outdir

        self.comp = conf.c_omp

        #--------------------------------------------------
        # hi-res data cube size
        self.stride_box = conf.box_stride
        if conf.twoD and self.stride_box != 1:
            raise Exception('stride_box needs to be 1 for 2D setup')

        self.nx_box = conf.box_nx
        self.ny_box = conf.Ny*conf.NyMesh
        self.nz_box = conf.Nz*conf.NzMesh
        self.fld_data = np.zeros((
            int(self.nx_box/self.stride_box), 
            int(self.ny_box/self.stride_box), 
            int(self.nz_box/self.stride_box), 
            10), 
            dtype='float32',
            )

        #--------------------------------------------------
        # prtcl momentum grid
        self.n_mom_grid = 400
        self.prtcl_data = np.zeros((self.nx_box, self.n_mom_grid, 8), dtype='float32',)

        # energy pspectra dN/dlog\gam
        self.uvel_max = 1.0e3
        self.uvel_min = 1.0e-3
        self.bins  = np.logspace( np.log10(self.uvel_min), np.log10(self.uvel_max), self.n_mom_grid+1)

        # component-wise four velocity pspectras dN/d\beta_i \gam
        self.uveli_max = 500.0
        self.binsd = np.linspace( -self.uveli_max, self.uveli_max, self.n_mom_grid+1)

        # x grid inside the box
        self.xbin = np.linspace(0.0, conf.box_nx, conf.box_nx)


        #--------------------------------------------------
        #location ahead of the shock
        self.x_box = conf.box_shift #in units of cells

        # global shock front location index
        self.shock_loc = 15

        # min and max tiles of the box
        self.wloc0 = 0
        self.wloc0 = 1

        # current box limits
        self.bloc0 = 0
        self.bloc1 = 1

    # integrate density along the y and z axis projecting it to x axis
    def get_density_profile(self, grid, lap, conf):

        # initialize arrays
        self.density_rank[:] = 0
        self.density[:] = 0

        for tile in pytools.tiles_local(grid):
            i,j,k = pytools.get_index(tile, conf)

            mins = tile.mins
            maxs = tile.maxs

            # discard tiles behind/infront the box
            #if not( mins[0] > self.wloc0-2*conf.NxMesh and maxs[0] < self.wloc1+2*conf.NxMesh):
            #    continue

            # now get g grid and project to x axis
            g = tile.get_grids(0)

            density_slice_tile = np.zeros((conf.NxMesh), dtype='float32')
            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):
                        density_slice_tile[l] += g.rho[l,m,n]

            #normalize
            density_slice_tile[:] /= self.ny*self.nz

            # add to the global array
            iglob, jglob, kglob = pytools.ind2loc((i, j, k), (0, 0, 0), conf)

            i0 = int(iglob)
            i1 = int(iglob+conf.NxMesh)
            self.density_rank[i0:i1] += density_slice_tile[:]

        # reduce globally
        MPI.COMM_WORLD.Reduce(
                [self.density_rank, MPI.FLOAT],
                [self.density,      MPI.FLOAT],
                op = MPI.SUM,
                root = 0,
                #tag = 1,
                )
        #MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write


    # locate shock front after calling get_density_profile
    def find_shock_front(self, grid, lap, conf, xmin=0):

        #unpack where output; original return type is unstable
        def where_last(vec):
            a = np.argwhere(vec)
            return a[-1][0] if len(a) else 0

        #broadcast value to everybody
        buf = np.zeros((1),'i')

        if grid.rank() == 0:
            #find shock front location

            # upstream density units
            n0 = conf.ppc*2
            qe = np.abs(conf.qe)
            norm = n0 #*qe

            #v1: shock front tracking
            ind = where_last( self.density/norm > conf.shock_density_jump_thr)

            #v2 shock front tracking where we test also that values left of the index have n_d/n_u > thr
            #nu_nd = self.density/norm 
            #while True:
            #    ind = where_last( nu_nd > conf.shock_density_jump_thr )
            #    mean_dens_behind = np.mean(nu_nd[ind-3*conf.c_omp:ind]) # mean density value behind found location

            #    #print( 'ind at', ind, 'nu/nd:', nu_nd[ind], 'mean val behind:', mean_dens_behind)

            #    if mean_dens_behind > conf.shock_density_jump_thr:
            #        break
            #    elif ind == 0:
            #        break # unphysical situation
            #    else:
            #        nu_nd = nu_nd[:ind] # cut point out and analyze the next

            # try to predic location roughly if the above test does not work
            if ind == 0:
                ind = 1*self.shock_loc + conf.cfl*conf.interval*0.89
                #ind = 1*self.shock_loc + conf.cfl*conf.interval*0.50

            # cap to a minimum of reflector location
            if ind < xmin:
                ind = xmin

            self.shock_loc = ind


            #print('shock at ind {} | {} c/omp'.format(ind, ind/self.comp))

            buf[0] = self.shock_loc #send message for broadcasting

        MPI.COMM_WORLD.Bcast(
                [buf, MPI.INT], 
                root = 0,
                #tag = 2,
                )

        self.shock_loc = int( buf[0] )
        #print('shock loc is', self.shock_loc, grid.rank())

        #MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write



    # read hi-res field data into localized cube; then reduce it so that root has all the data
    def read_fields(self, grid, lap, conf):

        #current box limits
        self.bloc0 = self.shock_loc + self.x_box
        self.bloc1 = self.shock_loc + self.x_box + self.nx_box
        #print(self.bloc0, self.bloc1)
        num_elems = 0


        #initialize arrays
        self.fld_data[:,:,:,:] = 0.0

        for tile in pytools.tiles_local(grid):
            i,j,k = pytools.get_index(tile, conf)

            mins = tile.mins
            maxs = tile.maxs

            # discard tiles behind/infront the box
            if not( mins[0] > self.wloc0-2*conf.NxMesh and maxs[0] < self.wloc1+2*conf.NxMesh):
                continue

            # discard tiles not inside the hi-res box stripe
            #if not( mins[0] > self.bloc0-1 and maxs[0] < self.bloc1+1):
            #    continue

            # now get g grid and project to x axis
            g = tile.get_grids(0)

            for l in range(conf.NxMesh):
                iglob, jglob, kglob = pytools.ind2loc((i, j, k), (l, 0, 0), conf)
                if self.bloc0 <= iglob < self.bloc1:
                    #ir = int(iglob) #- int(self.bloc0)
                    ir = iglob-self.bloc0
                    jr = jglob
                    kr = kglob

                    num_elems += 1

                    # NOTE: i index already has l subindex baked into it via iglob
                    ist = ifloor(ir/self.stride_box) #real index with stride hopping
                    #ist = ifloor(ir/self.stride_box) - ifloor(self.bloc0/self.stride_box) #real index with stride hopping

                    #print('inserting iglob to ', iglob, ir, l)
                    for m in range(conf.NyMesh):
                        jst = ifloor((jr + m)/self.stride_box) #real index with stride hopping

                        for n in range(conf.NzMesh):
                            kst = ifloor((kr + n)/self.stride_box) #real index with stride hopping

                            self.fld_data[ist, jst, kst, 0] = g.rho[l,m,n]

                            self.fld_data[ist, jst, kst, 1] = g.ex[l,m,n]
                            self.fld_data[ist, jst, kst, 2] = g.ey[l,m,n]
                            self.fld_data[ist, jst, kst, 3] = g.ez[l,m,n]

                            self.fld_data[ist, jst, kst, 4] = g.bx[l,m,n]
                            self.fld_data[ist, jst, kst, 5] = g.by[l,m,n]
                            self.fld_data[ist, jst, kst, 6] = g.bz[l,m,n]

                            self.fld_data[ist, jst, kst, 7] = g.jx[l,m,n]
                            self.fld_data[ist, jst, kst, 8] = g.jy[l,m,n]
                            self.fld_data[ist, jst, kst, 9] = g.jz[l,m,n]
        

        # normalize stride summation away
        #self.fld_data[:,:,:, 0] /= self.stride_box**3 #rho
        #self.fld_data[:,:,:, 7] /= self.stride_box**3 #jx
        #self.fld_data[:,:,:, 8] /= self.stride_box**3 #jy
        #self.fld_data[:,:,:, 9] /= self.stride_box**3 #jz


        # sum reduce to root
        n = int( int(self.nx_box*self.ny_box*self.nz_box*10/(self.stride_box**3)) )
        rbuf = np.zeros((
            int(self.nx_box/self.stride_box), 
            int(self.ny_box/self.stride_box), 
            int(self.nz_box/self.stride_box), 
                10),
                dtype='float32')
        
        #print('collected ', num_elems, ' vs ', self.nx_box)

        MPI.COMM_WORLD.Reduce( 
                [self.fld_data, n, MPI.FLOAT],
                [rbuf, n, MPI.FLOAT],
                #op = MPI.SUM, 
                op = op_overlap,
                root = 0,
                #tag = 3,
                )

        self.fld_data[:] = rbuf[:] #copy reduced buffer back to original data array

        #MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write


    # read hi-res prtcl momentum data into localized cube; then reduce it so that root has all the data
    def read_prtcls(self, grid, lap, conf):

        #current box limits
        self.bloc0 = self.shock_loc + self.x_box
        self.bloc1 = self.shock_loc + self.x_box + self.nx_box
        #print(self.bloc0, self.bloc1)

        xbin = np.linspace(self.bloc0, self.bloc1+1, conf.box_nx+1)
        #print('xgrid:', xbin)

        num_elems = 0

        #initialize arrays
        self.prtcl_data[:,:,:] = 0.0

        for tile in pytools.tiles_local(grid):
            i,j,k = pytools.get_index(tile, conf)

            mins = tile.mins
            maxs = tile.maxs

            # discard tiles behind/infront the box
            if not( mins[0] > self.wloc0 - 2*conf.NxMesh and maxs[0] < self.wloc1 + 2*conf.NxMesh):
                continue

            # discard tiles not inside the hi-res box stripe
            #if not( mins[0] > self.bloc0-1 and maxs[0] < self.bloc1+1):
            #    continue

            for isp in [0,1]:
                cont = tile.get_container(isp)

                npp = cont.size
                #print('total number of prtcls in this tile:', npp)

                # get locations
                xp = np.array(cont.loc(0))

                # get four velocities
                uxp = np.array(cont.vel(0))
                uyp = np.array(cont.vel(1))
                uzp = np.array(cont.vel(2))

                gam = np.sqrt(1.0 + uxp**2 + uyp**2 + uzp**2)
                #betax = uxp/gam
                #betay = uyp/gam
                #betaz = uzp/gam

                #print(len(xp))
                #print(len(uxp))

                # bin prtcls to dlog(gam-1) grid 
                hist, bin1, bin2 = np.histogram2d(
                                    xp,
                                    gam-1.0,
                                    bins=(xbin, self.bins),
                                    )
                
                # bin prtcls to d \beta_x\gam grid 
                histx,bin1x,bin2x= np.histogram2d(
                                    xp,
                                    uxp,
                                    bins=(xbin, self.binsd),
                                    )

                # bin prtcls to d \beta_y\gam grid 
                histy,bin1y,bin2y= np.histogram2d(
                                    xp,
                                    uyp,
                                    bins=(xbin, self.binsd),
                                    )

                # bin prtcls to d \beta_y\gam grid 
                histz,bin1z,bin2z= np.histogram2d(
                                    xp,
                                    uzp,
                                    bins=(xbin, self.binsd),
                                    )

                #and append 
                if isp == 0: #electrons
                    self.prtcl_data[:,:,0] += hist[:,:]
                    self.prtcl_data[:,:,1] += histx[:,:]
                    self.prtcl_data[:,:,2] += histy[:,:]
                    self.prtcl_data[:,:,3] += histz[:,:]
                elif isp == 1: #positrons          [:,:]
                    self.prtcl_data[:,:,4] += hist[:,:]
                    self.prtcl_data[:,:,5] += histx[:,:]
                    self.prtcl_data[:,:,6] += histy[:,:]
                    self.prtcl_data[:,:,7] += histz[:,:]


        #--------------------------------------------------
        # end of loop over tiles

        # sum reduce to root
        n = int(self.nx_box*self.n_mom_grid*8)
        rbuf = np.zeros((self.nx_box, self.n_mom_grid, 8), dtype='float32')
        
        #print('collected ', num_elems, ' vs ', self.nx_box)

        MPI.COMM_WORLD.Reduce( 
                [self.prtcl_data, n, MPI.FLOAT],
                [rbuf, n, MPI.FLOAT],
                op = MPI.SUM, 
                root = 0,
                #tag = 4,
                )

        self.prtcl_data[:] = rbuf[:] #copy reduced buffer back to original data array

        #MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write


    # save output to disk
    def save(self, grid, lap, conf):

        if grid.rank() == 0: 
            fname = conf.outdir + '/shock_{}.h5'.format(str(lap)) #.rjust(4,'0'))
            #print(fname)
            f5 = h5.File(fname,'w')
            dset = f5.create_dataset('dens_profile', data=self.density)
            f5.attrs['shock_loc'] = self.shock_loc

            #--------------------------------------------------
            # save hi-res data cube
            dset0 = f5.create_dataset('rho', data=self.fld_data[:,:,:,0])

            dset1 = f5.create_dataset('ex', data=self.fld_data[:,:,:,1])
            dset2 = f5.create_dataset('ey', data=self.fld_data[:,:,:,2])
            dset3 = f5.create_dataset('ez', data=self.fld_data[:,:,:,3])

            dset4 = f5.create_dataset('bx', data=self.fld_data[:,:,:,4])
            dset5 = f5.create_dataset('by', data=self.fld_data[:,:,:,5])
            dset6 = f5.create_dataset('bz', data=self.fld_data[:,:,:,6])

            dset7 = f5.create_dataset('jx', data=self.fld_data[:,:,:,7])
            dset8 = f5.create_dataset('jy', data=self.fld_data[:,:,:,8])
            dset9 = f5.create_dataset('jz', data=self.fld_data[:,:,:,9])


            #--------------------------------------------------
            #prtcl momentum space
            dsetp1 = f5.create_dataset('gam_0', data=self.prtcl_data[:,:,0])
            dsetp2 = f5.create_dataset('ux_0',  data=self.prtcl_data[:,:,1])
            dsetp3 = f5.create_dataset('uy_0',  data=self.prtcl_data[:,:,2])
            dsetp4 = f5.create_dataset('uz_0',  data=self.prtcl_data[:,:,3])

            dsetp5 = f5.create_dataset('gam_1', data=self.prtcl_data[:,:,4])
            dsetp6 = f5.create_dataset('ux_1',  data=self.prtcl_data[:,:,5])
            dsetp7 = f5.create_dataset('uy_1',  data=self.prtcl_data[:,:,6])
            dsetp8 = f5.create_dataset('uz_1',  data=self.prtcl_data[:,:,7])

            dsetu1 = f5.create_dataset('bin_gams',    data=self.bins)
            dsetu2 = f5.create_dataset('bin_betagam', data=self.binsd)


            #--------------------------------------------------
            # cube metadata
            f5.attrs['bloc0'] = self.bloc0
            f5.attrs['bloc1'] = self.bloc1

            f5.attrs['Nx'] = self.nx_box
            f5.attrs['Ny'] = self.ny_box
            f5.attrs['Nz'] = self.nz_box

            f5.attrs['stride'] = self.stride_box

            f5.close()

        return 


