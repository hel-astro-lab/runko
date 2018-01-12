from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys, os

import pyplasmaDev as pdev

from visualize import imshow



class conf:

    outdir = "out"

    Nxv = 10
    Nyv = 10
    Nzv = 10

    xmin = -2.0
    ymin = -3.0
    zmin = -4.0

    xmax =  2.0
    ymax =  3.0
    zmax =  4.0



class Mesh:
    xx = None
    yy = None
    ff = None




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



def level_fill(m, rfl):

    nx, ny, nz = m.get_length(rfl)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x,y,z = m.get_center([i,j,k], rfl)
                val = gauss(x,y,z)
                m[i,j,k, rfl] =  val


def get_mesh(m, args):
    rfl = args["rfl"]

    n1, n2, n3 = m.get_length(rfl)

    #flip dimensions to right-handed coordinate system
    if   args["dir"] == "xy" :
        nx = n1
        ny = n2
    elif args["dir"] == "xz":
        nx = n1
        ny = n3
    elif args["dir"] == "yz":
        nx = n2
        ny = n3

    #empty arrays ready
    xx = np.zeros((nx))
    yy = np.zeros((ny))
    ff = np.zeros((nx, ny))

    if  args["dir"] == "xy" :
        for i in range(nx):
            x,y,z = m.get_center([i,0,0], rfl)
            xx[i] = x
        for j in range(ny):
            x,y,z = m.get_center([0,j,0], rfl)
            yy[j] = y
    elif args["dir"] == "xz":
        for i in range(nx):
            x,y,z = m.get_center([i,0,0], rfl)
            xx[i] = x
        for j in range(ny):
            x,y,z = m.get_center([0,0,j], rfl)
            yy[j] = z
    elif args["dir"] == "yz":
        for i in range(nx):
            x,y,z = m.get_center([0,i,0], rfl)
            xx[i] = y
        for j in range(ny):
            x,y,z = m.get_center([0,0,j], rfl)
            yy[j] = z


    # collect values from mesh
    q = args["q"]
    for i in range(nx):
        for j in range(ny):
            if  args["dir"] == "xy" :
                val   = m[i, j, q, rfl]
            elif args["dir"] == "xz":
                val   = m[i, q, j, rfl]
            elif args["dir"] == "yz":
                val   = m[q, i, j, rfl]
            ff[i,j] = val

    m = Mesh()
    m.xx = xx
    m.yy = yy
    m.ff = ff

    return m


def plot2DSlice(ax, m, args):

    mesh = get_mesh(m, args)

    imshow(ax,
           mesh.ff,
           mesh.xx[0], mesh.xx[-1],
           mesh.yy[0], mesh.yy[-1],
           vmin = 0.0,
           vmax = 1.0,
           cmap = "plasma",
           clip = 0.0
           )

    return





    
def saveVisz(lap, conf):
    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/amr2_{}.png'.format(slap)
    plt.savefig(fname)




if __name__ == "__main__":

    # set up the grid
    ################################################## 
    m = pdev.AdaptiveMesh3D()
    m.resize( [conf.Nxv,  conf.Nyv,  conf.Nzv ])
    m.set_min([conf.xmin, conf.ymin, conf.zmin])
    m.set_max([conf.xmax, conf.ymax, conf.zmax])

    print("max. possible refinement:", m.get_maximum_possible_refinement_level())

    level_fill(m, 0)
    level_fill(m, 1)
    level_fill(m, 2)
    #level_fill(m, 3)



    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,12))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(3, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )


    args = {"dir":"xy", 
            "q":   20,
            "rfl": 2 }
    plot2DSlice(axs[0], m, args)


    args = {"dir":"xz", 
            "q":   20,
            "rfl": 2 }
    plot2DSlice(axs[1], m, args)


    args = {"dir":"yz", 
            "q":   20,
            "rfl": 2 }
    plot2DSlice(axs[2], m, args)



    saveVisz(0, conf)




