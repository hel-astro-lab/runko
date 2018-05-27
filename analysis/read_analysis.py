from __future__ import print_function

import numpy as np
import h5py

import corgi
import pyplasma as plasma

from read_mesh import TileInfo

# Read analysis lattice objects from h5 files
def get_analysis(f5, conf):

    # parse tile location
    i = conf.i
    j = conf.j
    k = conf.k
    ispcs = conf.ispcs


    fname = "analysis_"+str(i)+"_"+str(j)+"_"+str(k)+"-"+str(ispcs)
    dset = f5[fname]

    # initialize empty Yee lattice
    NxMesh = dset["Nx"].value
    NyMesh = dset["Ny"].value
    NzMesh = dset["Nz"].value
    aa = plasma.PlasmaMomentLattice(NxMesh, NyMesh, NzMesh)

    # read serialized 1D arrays and reshape into correct cubes
    rho    = dset["rho"].value.reshape((NxMesh, NyMesh, NzMesh))
    mgamma = dset["mgamma"].value.reshape((NxMesh, NyMesh, NzMesh))
    Vx     = dset["Vx"].value.reshape((NxMesh, NyMesh, NzMesh))
    Tx     = dset["Tx"].value.reshape((NxMesh, NyMesh, NzMesh))
    ekin   = dset["ekin"].value.reshape((NxMesh, NyMesh, NzMesh))

    # insert values into YeeLattice
    for s in range(NzMesh):
        for r in range(NyMesh):
            for q in range(NxMesh):
                aa.rho[q,r,s]    = rho[q,r,s]
                aa.mgamma[q,r,s] = mgamma[q,r,s]
                aa.Vx[q,r,s]     = Vx[q,r,s]
                aa.Tx[q,r,s]     = Tx[q,r,s]
                aa.ekin[q,r,s]   = ekin[q,r,s]

    return aa


# read full analysis data of the simulation snapshot into dictionary
def get_moments(prefix, ispcs, conf):

    # TODO fix hard coded limits
    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh

    #data
    data = {'x' : np.linspace(xmin, xmax, conf.Nx*conf.NxMesh),
           'rho':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'mgamma':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Vx':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Vy':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Vz':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Tx':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Ty':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Tz':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'ekin':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           }

    # initialize tile info
    tinfo = TileInfo()
    #tinfo.j = 0
    tinfo.j = 0
    tinfo.k = 0
    tinfo.ispcs = ispcs


    rank = 0 #TODO remove hard coded rank
    fname = prefix + "analysis-" + str(rank) + "_" + str(lap) + ".h5"
    f5 = h5py.File(fname,'r')
    print(fname)

    for i in range(conf.Nx):
        tinfo.i = i

        analysis = get_analysis(f5, tinfo)

        for s in range(conf.NxMesh):
            indx = i*conf.NxMesh + s

            data['rho'][indx] = analysis.rho[s, 0, 0]

            data['mgamma'][indx] = analysis.mgamma[s, 0, 0]

            data['Vx'][indx] = analysis.Vx[s, 0, 0]
            data['Vy'][indx] = analysis.Vy[s, 0, 0]
            data['Vz'][indx] = analysis.Vz[s, 0, 0]

            data['Tx'][indx] = analysis.Tx[s, 0, 0]
            data['Ty'][indx] = analysis.Ty[s, 0, 0]
            data['Tz'][indx] = analysis.Tz[s, 0, 0]

            data['ekin'][indx] = analysis.ekin[s, 0, 0]

    return data



# testing if run as main
if __name__ == "__main__":

    if False:
        lap = 0
        rank = 0
        fname = "../projects/tests/analysis-" + str(rank) + "_" + str(lap) + ".h5"

        f = h5py.File(fname,'r')

        tinfo = TileInfo()
        tinfo.i = 0
        tinfo.j = 0
        tinfo.k = 0
        tinfo.ispcs = 0

        aa = get_analysis(f, tinfo)


    if True:
        from configSetup import Configuration
        import matplotlib.pyplot as plt

        plt.fig = plt.figure(1, figsize=(8,3))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(1, 1)
        gs.update(hspace = 0.5)

        axs = []
        axs.append( plt.subplot(gs[0]) )

        lap = 0
        ispcs = 0
        prefix = "../projects/tests/"
        conf = Configuration(prefix + 'config-landau.ini') 

        analysis = get_moments(prefix, ispcs, conf)


        ###################################################
        # visualize

        axs[0].minorticks_on()
        axs[0].plot(analysis['x'], analysis['mgamma'], "b-")


        slap = str(lap).rjust(4, '0')
        fname = 'analysis-{}_{}.png'.format(0, slap)
        plt.savefig(fname)
