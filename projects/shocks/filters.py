import numpy as np

import pycorgi.twoD as corgi
import pyrunko.tools.twoD as pytools
import pyrunko.vlv.twoD as pyvlv
import pyrunko.pic.twoD as pypic
import pyrunko.fields.twoD as pyfld

import matplotlib.pyplot as plt

from visualize import imshow



np.random.seed(1)


class Conf:

    NxMesh = 100
    NyMesh = 100
    NzMesh = 1

    npasses = 20



def read_tile(tile, conf):

    nx,ny,nz = conf.NxMesh, conf.NyMesh, conf.NzMesh
    arr = np.zeros((nx,ny,nz))

    yee = tile.get_yee(0)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                arr[i,j,k] = yee.jx[i,j,k]

    return arr



if __name__ == "__main__":


    plt.fig = plt.figure(1, figsize=(5.0, 5.0), dpi=200)
    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(2, 2)
    #gs.update(hspace = 0.0)
    #gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    axs.append( plt.subplot(gs[1,0]) )
    axs.append( plt.subplot(gs[0,1]) )
    axs.append( plt.subplot(gs[1,1]) )



    conf = Conf()

    do_compensator = False

    #binomial + compensator filter
    if do_compensator:
        flt      = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        fltC     = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt.alpha  = 0.5 #make binomial
        fltC.alpha = conf.npasses/2 + 1 #binomial compensation is n/2 + 1
    else:
        #strided filters
        flt      = pyfld.General3pStrided(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt.alpha  = 0.5 #make binomial
        flt.stride = 3   #make binomial


    tile = pypic.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
    yee = tile.get_yee(0)
    yee.jx[50,50,0] = 1.0


    #--------------------------------------------------
    #sweep over npasses times
    for fj in range(conf.npasses):
        print("pass ", fj)
        flt.solve(tile)

    if do_compensator:
        print("compensator ")
        fltC.solve(tile)

    #--------------------------------------------------

    arr = read_tile(tile, conf)

    print(arr)

    x = np.arange(conf.NxMesh)
    
    axs[0].plot(x, arr[50,:,0], "k-", alpha=0.6)
    axs[0].plot(x, arr[:,50,0], "r--", alpha=0.6)


    xmin = 0
    ymin = 0
    xmax = conf.NxMesh
    ymax = conf.NyMesh

    im = imshow(axs[1], 
           arr[:,:,0], 
           xmin, xmax, ymin, ymax,
           cmap = 'viridis',
           vmin = 0,
           vmax = np.max(arr),
           clip = False,
           aspect=1,
           )

    #--------------------------------------------------
    # fft
    nx = conf.NxMesh//2

    freq = np.fft.fftfreq(conf.NxMesh, d = 1.0)
    farr = np.fft.fft2(arr[:,:,0])
    #farr = np.abs(farr)*np.abs(farr) #power
    farr = np.real(farr)**2

    im = imshow(axs[3], 
           farr[:nx,:nx],
           freq[0], freq[nx], freq[0], freq[nx],
           cmap = 'viridis',
           vmin = np.min(farr),
           vmax = np.max(farr),
           clip = False,
           aspect=1,
           )

    gainx = farr[:nx,0]
    gainy = farr[0,:nx]
    freqs = freq[:nx]

    axs[2].plot(freqs, gainx, "k-", alpha=0.6)
    axs[2].plot(freqs, gainy, "r--", alpha=0.6)

    axs[2].set_xscale('log')
    axs[2].set_yscale('log')

    axs[2].set_ylim((1.0e-11, 10))


    fdir = 'filters/'
    fname = fdir + 'test.pdf' 
    plt.savefig(fname)
