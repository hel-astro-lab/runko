from __future__ import print_function

import numpy as np
from scipy.stats import multivariate_normal
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

    xmin = -4.0
    ymin = -5.0
    zmin = -4.0

    xmax =  4.0
    ymax =  5.0
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
    #f = 1.0
    #f *= np.exp(-0.5*((ux - mux)**2)/delgam)
    #f *= np.exp(-0.5*((uy - muy)**2)/delgam)
    #f *= np.exp(-0.5*((uz - muz)**2)/delgam)

    mean = [0.0, 0.0, 0.0]
    cov  = np.zeros((3,3))
    cov[0,0] = 1.0
    cov[1,1] = 2.0
    cov[2,2] = 3.0

    xxx = [ux, uy, uz]

    f = multivariate_normal.pdf(xxx, mean, cov)

    return f



def level_fill(m, rfl):

    nx, ny, nz = m.get_length(rfl)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x,y,z = m.get_center([i,j,k], rfl)
                val = gauss(x,y,z)
                m[i,j,k, rfl] =  val


def adaptive_fill(m, tol=0.001):

    nx, ny, nz = m.get_length(0)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x,y,z = m.get_center([i,j,k], 0)
                val = gauss(x,y,z)
                m[i,j,k, 0] =  val

    adapter = pdev.Adapter();


    adapter.check(m)
    adapter.refine(m)
    for cid in adapter.cells_created:
        rfl = m.get_refinement_level(cid)
        indx = m.get_indices(cid)
        x,y,z = m.get_center(indx, rfl)

        print("new cell {} at ({},{},{})/{}".format(cid,
            indx[0],
            indx[1],
            indx[2],
            rfl))

        val = gauss(x,y,z)
        m[indx[0], indx[1], indx[2], rfl] = val
        



    #next refine grid until tolerance
    #for sweep in [0]:
    if False:

        for cid in m.get_cells(True):
            rfl = m.get_refinement_level(cid)

            if(not m.is_leaf(cid) ):
                continue

            indx = m.get_indices(cid)
            grad = pdev.grad(m, indx, rfl)

            # our error indicator is | max{ Grad(f) } |
            err = np.max( np.abs(grad))
            #print("refining err: {}".format(err))
            if (err > 0.001):
            #if True:
                #m.split(cid) #let's not remove the parent

                #fill children with correct values
                for cidc in m.get_children(cid):
                    indxc    = m.get_indices(cidc)
                    xc,yc,zc = m.get_center(indxc, rfl+1)

                    val = 0.8*gauss(xc,yc,zc)
                    m[indxc[0], indxc[1], indxc[2], rfl+1] = val

                #m.update_from_children(cid)

        #induced refining goes here





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


    if args["q"] == "mid":
        if args["dir"] == "xy":
            q = n3/2
        elif args["dir"] == "xz":
            q = n2/2
        elif args["dir"] == "yz":
            q = n1/2
    else:
        q = args["q"]


    # collect values from mesh
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


    #normalize
    mesh.ff = mesh.ff / np.max(mesh.ff)

    imshow(ax,
           mesh.ff,
           mesh.xx[0], mesh.xx[-1],
           mesh.yy[0], mesh.yy[-1],
           vmin = 0.0,
           vmax = 1.0,
           cmap = "plasma_r",
           clip = 0.0
           )
    return


def plotAdaptiveSlice(ax, m, args):

    rfl_max = args["rfl"]
    nx, ny, nz = m.get_length(rfl_max)

    #empty arrays ready
    xx = np.zeros((nx))
    yy = np.zeros((ny))
    ff = np.zeros((nx, ny))

    for i in range(nx):
        x,y,z = m.get_center([i,0,0], rfl_max)
        xx[i] = x
    for j in range(ny):
        x,y,z = m.get_center([0,j,0], rfl_max)
        yy[j] = y


    if args["q"] == "mid":
            q = nz/2
    else:
        q = args["q"]


    #adapter = pdev.Adapter();
    #adapter.check(m)
    #cells = adapter.cells_to_refine
    #cells = m.get_cells(True)
    #for cid in cells:
    #        rfli = m.get_refinement_level(cid)
    #        indx  = m.get_indices(cid)
    #        i,j,k = indx
    #        if k != q:
    #            continue
    #        if rfli != rfl_max:
    #            continue
    #        val = m[i,j,k, rfli]
    #        ff[i,j] = val


    for i in range(nx):
        for j in range(ny):
            indx = [i,j,q]
            val = m[i,j,q,rfl_max]

            ff[i,j] = val





    imshow(ax,
           ff,
           xx[0], xx[-1],
           yy[0], yy[-1],
           vmin = 0.0,
           vmax = 0.02,
           cmap = "plasma_r",
           clip = 0.0
           )
    return





def get_guide_grid(m, rfl):
    nx, ny, nz = m.get_length(rfl)
    xx = np.zeros((nx))
    yy = np.zeros((ny))
    zz = np.zeros((nz))

    #get guide grids
    for i in range(nx):
        x,y,z = m.get_center([i,0,0], rfl)
        xx[i] = x
    for i in range(ny):
        x,y,z = m.get_center([0,i,0], rfl)
        yy[i] = y
    for i in range(nz):
        x,y,z = m.get_center([0,0,i], rfl)
        zz[i] = z

    return xx, yy, zz
    


def get_gradient(m, rfl):
    nx, ny, nz = m.get_length(rfl)

    #empty arrays ready
    ggg = np.zeros((nx, ny, nz, 3))

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                indx = [i,j,k]

                gr = pdev.grad(m, indx, rfl)
                ggg[i,j,k,:] = gr

    return ggg



def plotGradientSlice(ax, m, args):
    rfl = args["rfl"]

    nx, ny, nz = m.get_length(rfl)
    xx, yy, zz = get_guide_grid(m, rfl)
    gg = np.zeros((nx, ny, 2))

    #get 3D gradient cube and flatten into 2D
    ggg = get_gradient(m, rfl)

    if args["q"] == "mid":
        q = nz/2
        print(q)
    else:
        q = args["q"]
    for i in range(nx):
        for j in range(ny):
            gg[i,j,0] = ggg[i,j,q,0] #vx
            gg[i,j,1] = ggg[i,j,q,1] #vy


    X, Y  = np.meshgrid(xx, yy)
    U     = gg[:,:,0]
    V     = gg[:,:,1]
    #speed = np.sqrt(U*U + V*V)
    #speed = np.sqrt(U*U)
    #speed = np.sqrt(V*V)

    speed = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            speed[i,j] = np.maximum( np.abs(U[i,j]), np.abs(V[i,j]) )

        
    #errs = np.max( np.abs(speed) )
 
    imshow(ax,
           speed / speed.max(),
           xx[0], xx[-1],
           yy[0], yy[-1],
           vmin = 0.0,
           vmax = 1.0,
           cmap = "plasma_r",
           clip = 0.0
           )


    lw = 5*speed / speed.max()
    ax.streamplot(X, Y, U, V, density=0.9, color="k", linewidth=lw)

    #ax.quiver(X, Y, U, V, units='x', pivot='tip', width=0.022) #, scale=1/0.15)
    
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

    #level_fill(m, 0)
    #level_fill(m, 1)
    #level_fill(m, 2)
    #level_fill(m, 3)

    adaptive_fill(m)



    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,20))
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
            "q":  "mid",
            "rfl": 2 }
    #plot2DSlice(axs[0], m, args)
    #plotAdaptiveSlice(axs[0], m, args)
    plotGradientSlice(axs[0], m, args)

    args = {"dir":"xz", 
            "q":   "mid",
            "rfl": 2 }
    plot2DSlice(axs[1], m, args)


    args = {"dir":"yz", 
            "q":   "mid",
            "rfl": 2 }
    plot2DSlice(axs[2], m, args)





    indx = [5,29,3]
    print(indx) 
    pindx1 = m.get_parent_indices(indx)
    print(pindx1)
    pindx2 = m.get_parent_indices(pindx1)
    print("1:", pindx2)

    print("2:",m.get_level_0_parent_indices(indx, 2))
    

    print(indx) 
    indx= m.get_parent_indices(indx)
    print(indx)
    indx = m.get_parent_indices(indx)
    print("3:", indx)
    
    #m.cut()



    saveVisz(0, conf)




