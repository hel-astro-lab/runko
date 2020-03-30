import numpy as np
import h5py as h5

class TestParticles:

    def __init__(self):
        self.prtcls = {}
        self.times = []


    def add_file(self, info):

        fname = info['particle_file'] + "_" + str(info['lap']) + ".h5"

        try:
            f5F = h5.File(fname,'r')
        except:
            return


        nx = f5F['Nx'].value
        ny = f5F['Ny'].value
        nz = f5F['Nz'].value
        print("reshaping 1D array into multiD with ({} {} {}) {}".format(nx,ny,nz, info['lap']))

        xloc = np.reshape( f5F['x'][:]    ,(nx,-1))
        yloc = np.reshape( f5F['y'][:]    ,(nx,-1)) 
        zloc = np.reshape( f5F['z'][:]    ,(nx,-1)) 
        ux   = np.reshape( f5F['vx'][:]   ,(nx,-1)) 
        uy   = np.reshape( f5F['vy'][:]   ,(nx,-1)) 
        uz   = np.reshape( f5F['vz'][:]   ,(nx,-1)) 
        wgt  = np.reshape( f5F['wgt'][:]  ,(nx,-1)) 
        ids  = np.reshape( f5F['id'][:]   ,(nx,-1)) 
        procs= np.reshape( f5F['proc'][:] ,(nx,-1))

        print("non zero prtcls:", np.count_nonzero(xloc))

        nt, npp = np.shape(xloc)
        for i in range(npp):
            key = (ids[0,i], procs[0,i])

            if not(key in self.prtcls):
                self.prtcls[key] = []

            for t in range(nt):
                gam = np.sqrt(1.0 + ux[t,i]**2 + uy[t,i]**2 + uz[t,i]**2)

                self.prtcls[key].append([
                        xloc[t,i],
                        yloc[t,i],
                        zloc[t,i],
                        ux[t,i],
                        uy[t,i],
                        uz[t,i],
                        wgt[t,i],
                        gam])

        print("adding {}".format(info['lap']))
        self.times.append(info['lap'])
