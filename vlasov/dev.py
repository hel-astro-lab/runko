from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys, os

import pyplasmaDev as pdev


class conf:

    outdir = "out"

    Nxv = 11
    Nyv = 21
    Nzv = 31

    xmin = -1.0
    ymin = -2.0
    zmin = -3.0

    xmax =  1.0
    ymax =  2.0
    zmax =  3.0



def gauss(ux,uy,uz):

    #return ux + uy + uz

    delgam = np.sqrt(1.0)
    mux = 0.0
    muy = 0.0
    muz = 0.0

    #f  = 1.0/np.sqrt(2.0*np.pi*delgam)
    f = 1.0
    f *= np.exp(-0.5*((ux - mux)**2)/delgam)
    f *= np.exp(-0.5*((uy - muy)**2)/delgam)
    f *= np.exp(-0.5*((uz - muz)**2)/delgam)

    return f
    #return 1.0



def fill(m):

    #rfl = 1 #refinement level
    #for rfl in range(m.maximum_refinement_level):
    for rfl in range(3):
        nx, ny, nz = m.get_length(rfl)

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x,y,z = m.get_center([i,j,k], rfl)
                    val = gauss(x,y,z)
                    m[i,j,k, rfl] =  val



def plotSlices(axs, m, rfl):


    nx, ny, nz = m.get_length(rfl)

    xmid = nx/2
    ymid = ny/2
    zmid = nz/2

    print("xmid: ", xmid)
    print("ymid: ", ymid)
    print("zmid: ", zmid)

    xx = np.zeros((nx))
    yy = np.zeros((ny))
    zz = np.zeros((nz))

    fx = np.zeros((nx))
    fy = np.zeros((ny))
    fz = np.zeros((nz))

    for i in range(nx):
        x,y,z = m.get_center([i,ymid,zmid], rfl)
        val   = m[i, ymid, zmid, rfl]

        xx[i] = x
        fx[i] = val

    for j in range(ny):
        x,y,z = m.get_center([xmid,j,zmid], rfl)
        val   = m[xmid, j, zmid, rfl]

        yy[j] = y
        fy[j] = val

    for k in range(nz):
        x,y,z = m.get_center([xmid,ymid,k], rfl)
        val   = m[xmid, ymid, k, rfl]

        zz[k] = z
        fz[k] = val

    cols = ["k", "b", "r", "g"]

    #axs[0].step(xx, fx, ".-", color=cols[rfl])

    axs[0].step(xx, fx, where="mid", color=cols[rfl])
    axs[1].step(yy, fy, where="mid", color=cols[rfl])
    axs[2].step(zz, fz, where="mid", color=cols[rfl])


    print(np.diff(xx))
    print(np.diff(yy))
    print(np.diff(zz))



    
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

    print("max. possible refinmenet:", m.get_maximum_possible_refinement_level())

    fill(m)


    #print(m.length)

    print("getting cids")
    print(m.get_length(0))
    print(m.get_cell([0,0,0], 0))
    print(m.get_cell([1,0,0], 0))
    print(m.get_cell([2,0,0], 0))

    print(m.get_length(1))
    print(m.get_cell([0,0,0], 1))
    print(m.get_cell([1,0,0], 1))
    print(m.get_cell([2,0,0], 1))


    print(m.get_length(2))
    print(m.get_cell([0,0,0], 2))
    print(m.get_cell([1,0,0], 2))
    print(m.get_cell([2,0,0], 2))
    print(m.get_cell([3,0,0], 2))
    print(m.get_cell([4,0,0], 2))

    #sys.exit()

    #m[0,0,0] = 1.0
    #print("000 = 1", m[0,0,0])
    #m[1,2,3] = 123.0
    #print("123 = 123", m[1,2,3])
    #print("000 = 1",   m[0,0,0])

    #print("get_level_0_cell_length")
    #print(m.get_level_0_cell_length())

    #print("beginning of the grid:")
    #print(m.get_center([0,0,0], 0))
    #print(m.get_center([1,1,1], 0))
    #print(m.get_center([2,2,2], 0))

    #print("end of the grid:")
    #print(m.get_center([conf.Nxv-1,conf.Nyv-1,conf.Nzv-1], 0))


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(3, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )

    plotSlices(axs, m, 0)
    plotSlices(axs, m, 1)
    plotSlices(axs, m, 2)
    #plotSlices(axs, m, 3)


    saveVisz(0, conf)









