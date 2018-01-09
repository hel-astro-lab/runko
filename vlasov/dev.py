from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys, os

import pyplasmaDev as pdev


class conf:

    outdir = "out"

    Nxv = 3
    Nyv = 3
    Nzv = 3

    xmin = -1.0
    ymin = -1.0
    zmin = -1.0

    xmax =  1.0
    ymax =  1.0
    zmax =  1.0



def gauss(ux,uy,uz):

    return ux + uy + uz

    delgam = np.sqrt(2.0)
    mux = 0.0
    muy = 0.0
    muz = 0.0

    #f  = 1.0/np.sqrt(2.0*np.pi*delgam)
    f = 1.0
    f *= np.exp(-0.5*((ux - mux)**2)/delgam)
    #f *= np.exp(-0.5*((uy - muy)**2)/delgam)
    #f *= np.exp(-0.5*((uz - muz)**2)/delgam)

    return f
    #return 1.0



def fill(m):

    rfl = 0 #refinement level
    nx, ny, nz = m.get_length(rfl)

    print("setting:")
    #for i in range(nx-1):
    #    for j in range(ny-1):
    #        for k in range(nz-1):
    #            x,y,z = m.get_center([i,j,k], rfl)
    #            m[rfl,i,j,k] = gauss(x,y,z)
    #            #m[rfl,i,j,k] = np.random.rand()

    #            print(m[0,i,j,k])
    m[0,1,1,1] = -1.0


    print("getting:")
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                print(m[0,i,j,k])



def plotXslice(ax, m, j, k):

    rfl = 0 #refinement level
    nx, ny, nz = m.get_length(rfl)

    xx = np.zeros((nx))
    yy = np.zeros((nx))

    for i in range(nx):
        x,y,z = m.get_center([i,0,0], rfl)
        val   = m[rfl, i, 0, 0]
        #val = gauss(x,y,z)

        xx[i] = x
        yy[i] = val

    print("xx:")
    print(xx)

    print("yy:")
    print(yy)

    ax.plot(xx, yy)


    
def saveVisz(lap, conf):
    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/amr_{}.png'.format(slap)
    plt.savefig(fname)




if __name__ == "__main__":

    # set up the grid
    ################################################## 
    m = pdev.AdaptiveMesh3D()
    m.resize( [conf.Nxv,  conf.Nyv,  conf.Nzv])
    m.set_min([conf.xmin, conf.ymin, conf.zmin])
    m.set_max([conf.xmax, conf.ymax, conf.zmax])

    fill(m)

    print("get_level_0_cell_length")
    print(m.get_level_0_cell_length())

    print("beginning of the grid:")
    print(m.get_center([0,0,0], 0))
    print(m.get_center([1,1,1], 0))
    print(m.get_center([2,2,2], 0))


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(5, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )

    plotXslice(axs[0], m, 0, 0)
    saveVisz(0, conf)









