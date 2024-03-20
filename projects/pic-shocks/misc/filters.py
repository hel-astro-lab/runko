import numpy as np

import pycorgi.twoD as corgi
import pyrunko.tools.twoD as pytools
import pyrunko.vlv.twoD as pyvlv
import pyrunko.pic.twoD as pypic
import pyrunko.emf.twoD as pyfld

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

    if True:
        yee = tile.get_yee(0)
        yee.jx[50,50,0] = 1.0

    #inject cherenkov like noise by hand
    elif False:
        yee = tile.get_yee(0)
        #yee.jx[50,50,0] = 1.0
        
        N = 200
        iss1 = np.random.randint(10, 89, size=N)
        jss1 = np.random.randint(10, 89, size=N)
        rss1 = np.random.randint(1, 10, size=N)

        N = 200
        iss2 = np.random.randint(1, 100, size=N)
        jss2 = np.random.randint(1, 100, size=N)

        #insert diagonal noise
        if True:
            for (i,j,r) in zip(iss1,jss1,rss1):
                yee.jx[i,j,0    ] += 1.0
                for ri in range(r):
                    yee.jx[i+ri,j+ri,0] += 1.0

        if False:
            for (i,j) in zip(iss1,jss1):
                yee.jx[i,j,0    ] = 1.0
                yee.jx[i+1,j+1,0] = 1.0

        #insert diagonal noise
        if False:
            for (i,j) in zip(iss2,jss2):
                yee.jx[i,j,0    ] = 1.0
                yee.jx[i-1,j+1,0] = 1.0

    timer = Timer()
    timer.start("filter")
    timer.do_print = True

    mode = 2

    #binomial + compensator filter
    if mode == 1:
        conf.npasses = 32
        #conf.npasses = 1

        print("Binomial+compensator filter")
        flt      = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        fltC     = pyfld.General3p(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        fltC2    = pyfld.Compensator2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        flt.alpha  = 0.5 #make binomial
        fltC.alpha = conf.npasses/2 + 1 #binomial compensation is n/2 + 1

        for fj in range(conf.npasses):
            flt.solve(tile)
        #fltC.solve(tile)
        #fltC2.solve(tile) #Modified compensator

    #strided filters
    elif mode == 2:
        conf.npasses = 1
        print("Strided filter")
        flt      = pyfld.General3pStrided(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        flt.alpha  = 0.5 #make binomial
        flt.stride = 2   #make binomial

        fltC2 = pyfld.Compensator2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

        for fj in range(conf.npasses):
            flt.solve(tile)
        fltC2.solve(tile) #Modified compensator


    #strided filters
    elif mode == 3:
        conf.npasses = 1
        print("outrolled strided filter")
        flt      = pyfld.Binomial2Strided2(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        for fj in range(conf.npasses):
            flt.solve(tile)

    #strided filters + binomial
    elif mode == 4:
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
    xmid = 50

    #print outrolled filter
    if False:
        sk = 5
        normf = 16*16
        print(normf*arr[(xmid-sk+1):(xmid+sk), (xmid-sk+1):(xmid+sk),0])

        for j in range(-3, 4):
            for i in range(-3, 4):
                val = arr[xmid+i, xmid+j,0]*16*16
                print("mesh.jx(i{:+d}, j{:+d}, 0)*{}*wn +".format(i,j,val))


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
           np.log10(farr[:nx+1,:nx+1]),
           freq[0], freq[nx], freq[0], freq[nx],
           cmap = 'viridis',
           vmin = -11.0,
           vmax = 1.0,
           clip = None,
           aspect=1,
           )

    gainx = farr[:nx+1,0]
    gainy = farr[0,:nx+1]
    freqs = freq[:nx+1]
    ks = 1.0/freqs

    ksd = np.sqrt(2.0*ks*ks)
    gaind = farr[:nx+1,:].diagonal()

    axs[2].plot(ks, gainx, "k-", alpha=0.6)
    axs[2].plot(ks, gainy, "r--", alpha=0.6)
    axs[2].plot(ksd, gaind, "b:", alpha=0.6)

    axs[2].set_xscale('log')
    axs[2].set_yscale('log')

    axs[2].set_ylim((1.0e-11, 10))


    fdir = 'filters/'
    fname = fdir + 'test.pdf' 
    plt.savefig(fname)
