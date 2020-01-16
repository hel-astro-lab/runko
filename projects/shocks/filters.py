import numpy as np

import pycorgi.twoD as corgi
import pyrunko.tools.twoD as pytools
import pyrunko.vlv.twoD as pyvlv
import pyrunko.pic.twoD as pypic
import pyrunko.fields.twoD as pyfld

import matplotlib.pyplot as plt

from visualize import imshow

from timer import Timer


np.random.seed(1)


class Conf:

    NxMesh = 100
    NyMesh = 100
    NzMesh = 1

    npasses = 60


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
    tile = pypic.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
    yee = tile.get_yee(0)
    yee.jx[50,50,0] = 1.0

    timer = Timer()
    timer.start("filter")
    timer.do_print = True


    #binomial + compensator filter
    if False:
        conf.npasses = 60
        print("Binomial+compensator filter")
        flt      = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        fltC     = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt.alpha  = 0.5 #make binomial
        fltC.alpha = conf.npasses/2 + 1 #binomial compensation is n/2 + 1

        for fj in range(conf.npasses):
            flt.solve(tile)
        fltC.solve(tile)

    #strided filters
    elif True:
        conf.npasses = 5
        print("Strided filter")
        flt      = pyfld.General3pStrided(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt.alpha  = 0.5 #make binomial
        flt.stride = 3   #make binomial

        for fj in range(conf.npasses):
            flt.solve(tile)

    #strided filters + binomial
    elif False:
        print("Strided + binomial")
        flt      = pyfld.General3pStrided(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt.alpha  = 0.5 #make binomial
        flt.stride = 3   #make binomial

        flt2     = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt2.alpha  = 0.5 #make binomial

        fltC     = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        fltC.alpha = conf.npasses/2 + 1 #binomial compensation is n/2 + 1

        #strided passes
        for fj in range(conf.npasses):
            flt.solve(tile)
        #extra binomial pass
        flt2.solve(tile)

    timer.stop("filter") 
    timer.stats("filter") 


    #--------------------------------------------------
    arr = read_tile(tile, conf)

    x = np.arange(conf.NxMesh)
    xd = np.sqrt(2.0*x*x)
    xmid = 50.0

    axs[0].plot(x - xmid, arr[50,:,0], "k-", alpha=0.6)
    axs[0].plot(x - xmid, arr[:,50,0], "r--", alpha=0.6)

    arrd = arr.diagonal()[0]
    axs[0].plot(xd-np.sqrt(2.0)*xmid, arrd, "b:", alpha=0.6)


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
    #farr = np.real(farr)**2
    farr = np.abs(farr)

    #im = imshow(axs[3], 
    #       farr[:nx,:nx],
    #       freq[0], freq[nx], freq[0], freq[nx],
    #       cmap = 'viridis',
    #       vmin = np.min(farr),
    #       vmax = np.max(farr),
    #       clip = False,
    #       aspect=1,
    #       )
    im = imshow(axs[3], 
           np.log10(farr[:nx,:nx]),
           freq[0], freq[nx], freq[0], freq[nx],
           cmap = 'viridis',
           vmin = -11.0,
           vmax = 1.0,
           clip = None,
           aspect=1,
           )

    gainx = farr[:nx,0]
    gainy = farr[0,:nx]
    freqs = freq[:nx]
    ks = 1.0/freqs

    ksd = np.sqrt(2.0*ks*ks)
    gaind = farr[:nx,:].diagonal()

    axs[2].plot(ks, gainx, "k-", alpha=0.6)
    axs[2].plot(ks, gainy, "r--", alpha=0.6)
    axs[2].plot(ksd, gaind, "b:", alpha=0.6)

    axs[2].set_xscale('log')
    axs[2].set_yscale('log')

    axs[2].set_ylim((1.0e-11, 10))


    fdir = 'filters/'
    fname = fdir + 'test.pdf' 
    plt.savefig(fname)
