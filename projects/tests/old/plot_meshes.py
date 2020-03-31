from __future__ import print_function

import numpy as np

from configSetup import Configuration

from visualize import imshow
from read_mesh import get_1d_meshes


# trick to make nice colorbars
# see http://joseph-long.com/writing/colorbars/
def colorbar(mappable, 
        loc="right", 
        orientation="vertical", 
        size="5%", 
        pad=0.05, 
        ticklocation='auto'):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(loc, size=size, pad=pad)
    return fig.colorbar(mappable, cax=cax, orientation=orientation, ticklocation=ticklocation)



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig = plt.figure(1, figsize=(6.974, 1.4))
    plt.rc('font', family='serif', size=8)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)
            
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    gs.update(wspace = 0.07)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )
    

    # what to plot
    ispcs = 0
    vdir = "x"

    prefix = 'twostream/out/'
    conf = Configuration('config-twostream-relativistic.ini') 

    #simulation box size
    xmin = 0.0
    ymin = 0.0
    xmax = conf.dx*conf.Nx*conf.NxMesh
    ymax = conf.dy*conf.Ny*conf.NyMesh
        

    # create panels from snapshots
    laps = range(0, 4001, 1000)
    for i, lap in enumerate(laps):
        print(lap)

        # check if this is top/bottom panel
        if i == 0:
            top = True
        else:
            top = False

        if i == len(laps)-1:
            bottom = True
        else:
            bottom = False

        # clear figure before we start
        for ax in axs:
            ax.clear()

        if top:
            axleft    = 0.10
            axbottom  = 0.24
            axright   = 0.96
            axtop     = 0.80
            fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

        useLog = True

        # read left
        ##################################################
        meshesL = get_1d_meshes(prefix, lap, conf, 0, vdir)


        print("atmosphere max dens: ", np.nanmax(meshesL['data']))
        norm_fac = np.nanmax(meshesL['data'])
        #print(norm_fac)
        meshesL['data'] = meshesL['data']/norm_fac

        if useLog:
            meshesL['data'] = np.log10(meshesL['data'])

        imL = imshow(axs[0], meshesL['data'],
               xmin, xmax,
               meshesL['vmin'], meshesL['vmax'],
               cmap = 'plasma_r',
               #vmin =   0.0,
               #vmax =   1.0,
               #clip =   1.0e-5,
               vmin =   -3.0,
               vmax =    0,
               clip = None
            )


        # read right
        ##################################################
        meshesR = get_1d_meshes(prefix, lap, conf, 1, vdir)


        print("beam max dens:", np.nanmax(meshesR['data']))
        #norm_fac = np.max(meshesR['data'])
        norm_fac = np.nanmax(meshesR['data'])
        meshesR['data'] = meshesR['data']/norm_fac
        print("relative beam max dens:", np.nanmax(meshesR['data']))

        if useLog:
            meshesR['data'] = np.log10(meshesR['data'])

        imR = imshow(axs[1], meshesR['data'],
               xmin, xmax,
               meshesR['vmin'], meshesR['vmax'],
               cmap = 'plasma_r',
               #vmin =   0.0,
               #vmax =   1.0,
               #clip =   1.0e-5,
               vmin =   -3.0,
               vmax =    0,
               clip = None
            )

        # remove ytick labels from second panel
        axs[1].set_yticklabels([])

        # set ylabel for leftmost panel
        axs[0].set_ylabel(r'$u_x$')


        # header panel
        if top:
            cbL = colorbar(imL, loc="top", orientation="horizontal", size="8%", pad=0.03, ticklocation='top')
            cbR = colorbar(imR, loc="top", orientation="horizontal", size="8%", pad=0.03, ticklocation='top')

        # middle panels
        if not(bottom):
            for ax in axs:
                ax.set_xticklabels([])

        # tail panel
        if bottom:
            for ax in axs:
                ax.set_xlabel(r'Location $x$ ($\lambda_p$)')

        # add timestamp
        
        tcur = lap*conf.dx*conf.cfl
        tstamp = "$t = $" + str(tcur) + " $\omega_p^{-1}$"
        txt = fig.text(axright+0.02, 0.5, tstamp, rotation=270, ha='center', va='center')



        slap = str(lap).rjust(5, '0')
        fname = prefix+'mesh_{}_{}.png'.format(0, slap)
        plt.savefig(fname, bbox_inches='tight')

        fname = prefix+'mesh_{}_{}.pdf'.format(0, slap)
        plt.savefig(fname, bbox_inches='tight')

        ##################################################
        # hacks to clean up the figure so we can re-use the same Fig object


        # finally we need to remove colorbar before next round
        if top:
            cbL.remove()
            cbR.remove()

        txt.remove()

