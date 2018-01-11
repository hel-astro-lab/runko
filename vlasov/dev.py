from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys, os

import pyplasmaDev as pdev


class conf:

    outdir = "out"

    Nxv = 10
    Nyv = 10
    Nzv = 23

    xmin = -2.0
    ymin = -3.0
    zmin = -4.0

    xmax =  2.0
    ymax =  3.0
    zmax =  4.0



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



def level_fill(m, rfl):

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

    #print("xmid: ", xmid)
    #print("ymid: ", ymid)
    #print("zmid: ", zmid)

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


    #print(np.diff(xx))
    #print(np.diff(yy))
    #print(np.diff(zz))


def plotHierarchy(ax, m):

    cols = ["black", "blue", "green"]
    for rfl in range(0,3):
        nx, ny, nz = m.get_length(rfl)
        xmid = nx/2
        ymid = ny/2
        zmid = nz/2
        xx = np.zeros((nx))
        fx = np.zeros((nx))
        for i in range(nx):
            indc = [i,ymid,zmid]
            x,y,z = m.get_center(indc, rfl)
            val   = m[i, ymid, zmid, rfl]

            xx[i] = x
            fx[i] = val
        ax.step(xx, fx, where="mid", color=cols[rfl])


    #next level connection
    
    q = 0
    for rfl in range(3,4):
        nx, ny, nz = m.get_length(rfl)
        xmid = nx/2
        ymid = ny/2
        zmid = nz/2

        for i in range(nx):
            indc  = [i,ymid,zmid]
            x,y,z = m.get_center(indc, rfl)
            f     = m[i, ymid, zmid, rfl]

            #indp    = m.get_parent_indices(indc)
            #xp,yp,zp = m.get_center(indp, rfl-1)
            #fp       = m[ indp[0], indp[1], indp[2], rfl-1]

            # zero level parent
            indp    = m.get_level_0_parent_indices(indc, rfl)
            xp,yp,zp = m.get_center(indp, 0)
            fp       = m[ indp[0], indp[1], indp[2], 0]

            ax.plot( [x, xp], [f, 0.0], "r-" )




    
def saveVisz(lap, conf):
    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/amr_{}.png'.format(slap)
    plt.savefig(fname)




if __name__ == "__main__":

    # set up the grid
    ################################################## 
    m = pdev.AdaptiveMesh3D()
    m.resize( [conf.Nxv,  conf.Nyv,  conf.Nzv ])
    m.set_min([conf.xmin, conf.ymin, conf.zmin])
    m.set_max([conf.xmax, conf.ymax, conf.zmax])

    cidc = m.get_cell([2,1,1],1)
    indc = m.get_indices(cidc)
    indp = m.get_parent_indices(indc)
    cidp = m.get_parent(cidc)

    print(cidc, " ", indc)
    print(cidp, " ", indp)

    
    print("max. possible refinement:", m.get_maximum_possible_refinement_level())

    level_fill(m, 0)
    level_fill(m, 1)
    level_fill(m, 2)
    level_fill(m, 3)



    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(4, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )
    axs.append( plt.subplot(gs[3]) )

    plotSlices(axs, m, 0)
    #plotSlices(axs, m, 1)
    plotSlices(axs, m, 3)

    plotHierarchy(axs[3], m)


    saveVisz(0, conf)




