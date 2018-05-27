from __future__ import print_function

import numpy as np
import h5py

import corgi
import pyplasma as plasma

from read_mesh import TileInfo



# Read Yee lattice objects from fields h5 files
def get_fields(f5, conf):

    # parse tile location
    i = conf.i
    j = conf.j
    k = conf.k

    fname = "yee_"+str(i)+"_"+str(j)+"_"+str(k)
    dset = f5[fname]

    # initialize empty Yee lattice
    NxMesh = dset["Nx"].value
    NyMesh = dset["Ny"].value
    NzMesh = dset["Nz"].value
    yee = plasma.YeeLattice(NxMesh, NyMesh, NzMesh)

    # read serialized 1D arrays and reshape into correct cubes
    exm  = dset["ex"].value.reshape((NxMesh, NyMesh, NzMesh))
    bxm  = dset["bx"].value.reshape((NxMesh, NyMesh, NzMesh))
    jxm  = dset["jx"].value.reshape((NxMesh, NyMesh, NzMesh))
    rhom = dset["rho"].value.reshape((NxMesh, NyMesh, NzMesh))

    # insert values into YeeLattice
    for s in range(NzMesh):
        for r in range(NyMesh):
            for q in range(NxMesh):
                yee.ex[q,r,s]  = exm[q,r,s]
                #yee.ey[q,r,s]  = eym[q,r,s]
                #yee.ez[q,r,s]  = ezm[q,r,s]

                yee.bx[q,r,s]  = bxm[q,r,s]
                #yee.by[q,r,s]  = bym[q,r,s]
                #yee.bz[q,r,s]  = bzm[q,r,s]

                yee.jx[q,r,s]  = jxm[q,r,s]
                #yee.jy[q,r,s]  = jym[q,r,s]
                #yee.jz[q,r,s]  = jzm[q,r,s]

                #yee.jx1[q,r,s] = jx1m[q,r,s]

                yee.rho[q,r,s] = rhom[q,r,s]

    return yee



# read full Yee data of the simulation snapshot into dictionary
def get_yee(prefix, conf):

    # TODO fix hard coded limits
    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh

    # initialize Yee data dictionary
    data = {'x' : np.linspace(xmin, xmax, conf.Nx*conf.NxMesh),
            'ex':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ey':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ez':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ez':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'bx':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'by':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'bz':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jx':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jy':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jz':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jx1':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'rho':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           }


    # initialize tile info
    tinfo = TileInfo()
    #tinfo.j = 0
    tinfo.j = 0
    tinfo.k = 0


    rank = 0 #TODO remove hard coded rank
    fname = prefix + "fields-" + str(rank) + "_" + str(lap) + ".h5"
    f5 = h5py.File(fname,'r')
    print(fname)

    for i in range(conf.Nx):
        tinfo.i = i
        yee = get_fields(f5, tinfo)
        for s in range(conf.NxMesh):
            indx = i*conf.NxMesh + s

            data['ex'][indx] = yee.ex[s, 0, 0]
            data['ey'][indx] = yee.ey[s, 0, 0]
            data['ez'][indx] = yee.ez[s, 0, 0]

            data['bx'][indx] = yee.bx[s, 0, 0]
            data['by'][indx] = yee.by[s, 0, 0]
            data['bz'][indx] = yee.bz[s, 0, 0]

            data['jx'][indx] = yee.jx[s, 0, 0]
            data['jy'][indx] = yee.jy[s, 0, 0]
            data['jz'][indx] = yee.jz[s, 0, 0]

            data['jx1'][indx] = yee.jx1[s, 0, 0]

            data['rho'][indx] = yee.rho[s, 0, 0]

    return data




# testing if run as main
if __name__ == "__main__":

    if False:
        lap = 0
        rank = 0
        fname = "../projects/tests/fields-" + str(rank) + "_" + str(lap) + ".h5"

        f = h5py.File(fname,'r')

        tinfo = TileInfo()
        tinfo.i = 0
        tinfo.j = 0
        tinfo.k = 0

        pm = get_fields(f, tinfo)



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
        prefix = "../projects/tests/"
        conf = Configuration(prefix + 'config-landau.ini') 

        yee = get_yee(prefix, conf)


        ###################################################
        # visualize

        axs[0].minorticks_on()
        axs[0].plot(yee['x'], yee['rho'], "b-")


        slap = str(lap).rjust(4, '0')
        fname = 'fields-{}_{}.png'.format(0, slap)
        plt.savefig(fname)












