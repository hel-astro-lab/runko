from __future__ import print_function

import numpy as np
import h5py

from configSetup import Configuration

from visualize import imshow
from read_mesh import get_mesh
from read_mesh import TileInfo

from visualize_amr import xSliceMid, ySliceMid, zSliceMid
from visualize_amr import get_leaf_mesh




# Read velocity mesh from given location
def get_vmesh_as_cube(tinfo):

    rank = 0 #TODO remove hard coded rank
    fname = prefix + "meshes-" + str(rank) + "_" + str(lap) + ".h5"
    f5 = h5py.File(fname,'r')
    print(fname)

    vmesh = get_mesh(f5, tinfo)
    pym = get_leaf_mesh(vmesh, conf.refinement_level)

    return pym







if __name__ == "__main__":

    from matplotlib.colors import Normalize
    from matplotlib.cm import get_cmap
    from matplotlib import colorbar

    import matplotlib.pyplot as plt
    #from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig = plt.figure(1, figsize=(3.487, 2.5))
            
    plt.rc('font',  family='serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    #gs.update(hspace = 0.5)
    #gs.update(wspace = 0.05)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )

    axleft    = 0.18
    axbottom  = 0.16
    axright   = 0.98
    axtop     = 0.80

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    for ax in axs:
        ax.minorticks_on()

        ax.set_xlabel(r"$u_x$")
        ax.set_ylabel(r"$f(u_x)$")

        ax.set_xlim(-0.3, 0.3)
        ax.set_ylim( 1.e-22, 1.0e-3)

        ax.set_yscale('log')


    # what to plot
    ispcs = 0
    vdir = "x"
    prefix = "../projects/tests/"
    conf = Configuration(prefix + 'config-landau.ini') 


    #simulation box size
    xmin = 0.0
    ymin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymax = conf.dy*conf.Ny*conf.NyMesh

    # location of the mesh
    tinfo = TileInfo()
    tinfo.i = 0
    tinfo.j = 0
    tinfo.k = 0
        
    tinfo.q = 0
    tinfo.r = 0
    tinfo.s = 0

    tinfo.ispcs = 0


    # create panels from snapshots
    laps = [0,1,2]

    # create color map
    norm = Normalize(vmin=0, vmax=len(laps))
    cmap = get_cmap('Spectral')


    # and put colorbar on top
    pos1 = axs[0].get_position()
    
    axwidth  = axright - axleft
    axheight = (axtop - axbottom)*0.05
    axpad = 0.02
    cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])

    cb1 = colorbar.ColorbarBase(
            cax,
            cmap=cmap,
            norm=norm,
            orientation='horizontal',
            ticklocation='top')

    cb1.set_label('Time $t$ ($\omega_p^{-1}$)')


    for i, lap in enumerate(laps):

        col = cmap(norm(i))

        # read left
        ##################################################
        vslice = get_vmesh_as_cube(tinfo)

        sl = xSliceMid(vslice) #slice from middle
        #sl = ySliceMid(vslice) #slice from middle
        #sl = zSliceMid(vslice) #slice from middle

        #uxs = np.linspace(vslice.xx
        print(vslice.xx)
        print(sl)

        axs[0].plot(vslice.xx, sl, color=col)


    fname = 'spectra.pdf'
    plt.savefig(fname)


