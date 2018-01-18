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

    Nxv = 15
    Nyv = 15
    Nzv = 15

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

    return 10.0 + 1.0*ux + 2.0*uy + 3.0*uz
    #return np.abs( 0.1 + 0.1*ux*uy + 0.3*ux*uz + 0.9*uy*ux + 2.0*ux*uy*uz ) + 0.1
    #return 0.1 + 0.2*ux + 0.3*uy + 0.4*uz + 0.5*ux*uy + 0.6*ux*uz + 0.7*uy*ux + 0.8*ux*uy*uz

    delgam = np.sqrt(1.0)
    mux = 0.0
    muy = 0.0
    muz = 0.0

    #f  = 1.0/np.sqrt(2.0*np.pi*delgam)
    #f = 1.0
    #f *= np.exp(-0.5*((ux - mux)**2)/delgam)
    #f *= np.exp(-0.5*((uy - muy)**2)/delgam)
    #f *= np.exp(-0.5*((uz - muz)**2)/delgam)

    mean = [0.0, 0.0, 1.0]
    cov  = np.zeros((3,3))
    cov[0,0] = 1.0 
    cov[1,1] = 2.0 
    cov[2,2] = 5.0


    xxx = [ux, uy, uz]

    f = multivariate_normal.pdf(xxx, mean, cov)


    return f



def level_fill(m, rfl):

    nx, ny, nz = m.get_size(rfl)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x,y,z = m.get_center([i,j,k], rfl)
                val = gauss(x,y,z)
                m[i,j,k, rfl] =  val


def adaptive_fill(m, tol=0.001, sweeps=4):

    nx, ny, nz = m.get_size(0)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x,y,z = m.get_center([i,j,k], 0)
                val = gauss(x,y,z)
                m[i,j,k, 0] =  val

    adapter = pdev.Adapter();
    
    sweep = 1
    while(True):
        print("-------round {}-----------".format(sweep))

        adapter.check(m)
        adapter.refine(m)

        print("cells to refine: {}".format( len(adapter.cells_to_refine)))
        for cid in adapter.cells_created:
            rfl = m.get_refinement_level(cid)
            indx = m.get_indices(cid)
            x,y,z = m.get_center(indx, rfl)

            #print("new cell {} at ({},{},{})/{}".format(cid,
            #    indx[0],
            #    indx[1],
            #    indx[2],
            #    rfl))

            val = gauss(x,y,z)
            m[indx[0], indx[1], indx[2], rfl] = val
            

        adapter.unrefine(m)
        print("cells to be removed: {}".format( len(adapter.cells_removed)))

        sweep += 1
        if sweep > sweeps: break
    #m.clip_cells(1.0e-5)



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

                    val = gauss(xc,yc,zc)
                    m[indxc[0], indxc[1], indxc[2], rfl+1] = val

                #m.update_from_children(cid)

        #induced refining goes here





def get_mesh(m, args):
    rfl = args["rfl"]
    n1, n2, n3 = m.get_size(rfl)

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



def get_real_coords(m, indx, bvals, rfl):
    x,y,z = m.get_center(indx, rfl)
    xl,yl,zl = m.get_length(rfl)

    xval = (x - xl*0.5) + bvals[0]*xl
    yval = (y - yl*0.5) + bvals[1]*yl
    zval = (z - zl*0.5) + bvals[2]*zl

    return xval, yval, zval



def get_interpolated_mesh(m, args):

    rfl_max = args["rfl"]
    reso = 10

    nx, ny, nz = m.get_size(rfl_max)
    lvlm = 2**rfl_max

    #empty arrays ready
    xx = np.zeros((nx))
    yy = np.zeros((ny))
    zz = np.zeros((nz))


    x,y,z = m.get_center([0,0,0], rfl_max)
    xx[0] = x
    yy[0] = y
    zz[0] = z

    x,y,z = m.get_center([nx,0,0], rfl_max)
    xx[-1]= x
    x,y,z = m.get_center([0,ny,0], rfl_max)
    yy[-1]= y
    x,y,z = m.get_center([0,0,nz], rfl_max)
    zz[-1]= z


    if args["q"] == "mid":
        if args["dir"] == "xy":
            q = nz/2
        elif args["dir"] == "xz":
            q = ny/2
        elif args["dir"] == "yz":
            q = nx/2
    else:
        q = args["q"]

    if args["dir"] == "xy":
        ff = np.zeros((nx*reso, ny*reso))
    elif args["dir"] == "xz":
        ff = np.zeros((nx*reso, nz*reso))
    elif args["dir"] == "yz":
        ff = np.zeros((ny*reso, nz*reso))


    cells = m.get_cells(True)
    leafs = 0
    for cid in cells:

        if(not m.is_leaf(cid) ):
            continue

        leafs += 1

        rfl = m.get_refinement_level(cid)
        lvl = 2**rfl

        [i,j,k] = m.get_indices(cid)

        st = lvlm/lvl #stretch factor for ref lvl

        if args["dir"] == "xy" and k*st != q:
            continue
        if args["dir"] == "xz" and j*st != q:
            continue
        if args["dir"] == "yz" and i*st != q:
            continue

        val = m[i,j,k, rfl]

        i1 = i * st*reso
        j1 = j * st*reso
        k1 = k * st*reso

        arr = np.linspace(0.0, 1.0, st*reso)

        if args["dir"] == "xy":
            for ii in range(st*reso):
                for jj in range(st*reso):
                    val = pdev.trilinear_interp(m, [i,j,k], [arr[ii], arr[jj],0.5], rfl)

                    xval, yval, zval = get_real_coords(m, [i,j,k], [arr[ii], arr[jj],0.5], rfl)
                    val2 = gauss(xval, yval, zval)

                    ff[i1+ii, j1+jj] = val/val2


        if args["dir"] == "xz":
            for ii in range(st*reso):
                for jj in range(st*reso):
                    val = pdev.trilinear_interp(m, [i,j,k], [arr[ii], 0.5, arr[jj]], rfl)
                    ff[i1+ii, k1+jj] = val

        if args["dir"] == "yz":
            for ii in range(st*reso):
                for jj in range(st*reso):
                    val = pdev.trilinear_interp(m, [i,j,k], [0.5, arr[ii], arr[jj]], rfl)
                    ff[j1+ii, k1+jj] = val


    m = Mesh()
    m.xx = xx
    m.yy = yy
    m.ff = ff

    return m




def get_leaf_mesh(m, args):

    rfl_max = args["rfl"]
    nx, ny, nz = m.get_size(rfl_max)
    lvlm = 2**rfl_max

    #empty arrays ready
    xx = np.zeros((nx))
    yy = np.zeros((ny))
    zz = np.zeros((nz))

    for i in range(nx):
        x,y,z = m.get_center([i,0,0], rfl_max)
        xx[i] = x
    for j in range(ny):
        x,y,z = m.get_center([0,j,0], rfl_max)
        yy[j] = y
    for k in range(nz):
        x,y,z = m.get_center([0,0,k], rfl_max)
        zz[j] = z


    if args["q"] == "mid":
        if args["dir"] == "xy":
            q = nz/2
        elif args["dir"] == "xz":
            q = ny/2
        elif args["dir"] == "yz":
            q = nx/2
    else:
        q = args["q"]

    if args["dir"] == "xy":
        ff = np.zeros((nx, ny))
    elif args["dir"] == "xz":
        ff = np.zeros((nx, nz))
    elif args["dir"] == "yz":
        ff = np.zeros((ny, nz))


    cells = m.get_cells(True)
    leafs = 0
    for cid in cells:

        if(not m.is_leaf(cid) ):
            continue

        leafs += 1

        rfl = m.get_refinement_level(cid)
        lvl = 2**rfl

        [i,j,k] = m.get_indices(cid)

        st = lvlm/lvl #stretch factor for ref lvl

        if args["dir"] == "xy" and k*st != q:
            continue
        if args["dir"] == "xz" and j*st != q:
            continue
        if args["dir"] == "yz" and i*st != q:
            continue

        val = m[i,j,k, rfl]

        i *= st
        j *= st
        k *= st

        if args["dir"] == "xy":
            for ii in range(st):
                for jj in range(st):
                    ff[i+ii, j+jj] = val
        if args["dir"] == "xz":
            for ii in range(st):
                for jj in range(st):
                    ff[i+ii, k+jj] = val
        if args["dir"] == "yz":
            for ii in range(st):
                for jj in range(st):
                    ff[j+ii, k+jj] = val


    Nc = len(cells)
    print("cells: {} / leafs {} (ratio: {}) / full {} (compression {})".format(
        Nc,
        leafs, 
        1.0*Nc/(1.0*leafs),
        nx*ny*nz, 
        1.0* Nc /(1.0*nx*ny*nz) 
        ))


    m = Mesh()
    m.xx = xx
    m.yy = yy
    m.ff = ff

    return m




def plot2DSlice(ax, m, args):

    #mesh = get_mesh(m, args)
    #mesh = get_leaf_mesh(m, args)
    mesh = get_interpolated_mesh(m, args)

    if True:
        imshow(ax,
               mesh.ff,
               mesh.xx[0], mesh.xx[-1],
               mesh.yy[0], mesh.yy[-1],
               vmin = 0.8,
               vmax = 1.2,
               cmap = "RdYlBu",
               clip = 0.0
               )

        return

    else:

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
    nx, ny, nz = m.get_size(rfl_max)

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
    nx, ny, nz = m.get_size(rfl)
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
    nx, ny, nz = m.get_size(rfl)

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

    nx, ny, nz = m.get_size(rfl)
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

    speed = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            speed[i,j] = np.maximum( np.abs(U[i,j]), np.abs(V[i,j]) )
        
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
    





def get_indicator(m, rfl):
    nx, ny, nz = m.get_size(rfl)

    #empty arrays ready
    ggg = np.zeros((nx, ny, nz))

    adapter = pdev.Adapter()
    adapter.set_maximum_data_value(1.0)


    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                indx = [i,j,k]

                cid = m.get_cell_from_indices(indx, rfl)

                gr = adapter.maximum_gradient(m, cid)
                #gr = adapter.maximum_value(m, cid)

                #gr = 1.0e-4 / (gr + 1.0e-10)


                ggg[i,j,k] = gr

    return ggg



def plotIndicator(ax, m, args):
    rfl = args["rfl"]

    nx, ny, nz = m.get_size(rfl)
    xx, yy, zz = get_guide_grid(m, rfl)
    gg         = np.zeros((nx, ny))

    #get 3D gradient cube and flatten into 2D
    ggg = get_indicator(m, rfl)


    if args["q"] == "mid":
        q = nz/2
        print(q)
    else:
        q = args["q"]

    for i in range(nx):
        for j in range(ny):
            gg[i,j] = ggg[i,j,q]



    X, Y  = np.meshgrid(xx, yy)
        
    print("maximum of refining indicator: {}", gg.max() )

    gg = np.log10(gg)
    #gg /= gg.max()

    imshow(ax,
           gg,
           xx[0], xx[-1],
           yy[0], yy[-1],
           #vmin = 0.0,
           #vmax = 1.0,
           vmin =-8.0,
           vmax = 2.0,
           cmap = "plasma_r",
           #clip = 0.0
           clip =-9.0
           )





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

    if False:
        level_fill(m, 0)
        level_fill(m, 1)
        #level_fill(m, 2)
        #level_fill(m, 3)

        #m.clip_cells(1.0e-4)

        #adapter = pdev.Adapter();
        #adapter.check(m)
        #adapter.unrefine(m)
        #print("cells to be removed: {}".format( len(adapter.cells_removed)))

    else:
        adaptive_fill(m, sweeps=2)



    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(12,20))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(3, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    axs.append( plt.subplot(gs[1,0]) )
    axs.append( plt.subplot(gs[2,0]) )

    axsE = []
    axsE.append( plt.subplot(gs[0,1]) )
    axsE.append( plt.subplot(gs[1,1]) )
    axsE.append( plt.subplot(gs[2,1]) )


    args = {"dir":"xy", 
            "q":  "mid",
            "rfl": 2 }
    plot2DSlice(axs[0], m, args)
    #plotAdaptiveSlice(axs[0], m, args)
    #plotGradientSlice(axsE[0], m, args)

    args["rfl"] = 1
    plotIndicator(axsE[0], m, args)


    args = {"dir":"xz", 
            "q":   "mid",
            "rfl": 2 }
    plot2DSlice(axs[1], m, args)


    args = {"dir":"yz", 
            "q":   "mid",
            "rfl": 2 }
    plot2DSlice(axs[2], m, args)


    

    for ax in axs:
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)
    for ax in axs:
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)
    for ax in axsE:
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)

    for ax in axsE:
        ax.set_xlim(-4.0, 4.0)
        ax.set_ylim(-4.0, 4.0)


    saveVisz(0, conf)




