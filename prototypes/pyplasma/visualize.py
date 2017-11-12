import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
cmap = pal.wesanderson.Moonrise1_5.mpl_colormap

import numpy as np


Nrank = 4


##################################################
# plotting tools

# visualize matrix
def imshow(ax, grid, xmin, xmax, ymin, ymax):

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-3.0, 3.0)

    extent = [ xmin, xmax, ymin, ymax ]

    mgrid = np.ma.masked_where(grid == -1.0, grid)
    
    mgrid = mgrid.T
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = 0.0,
              vmax = Nrank-1,
              aspect='auto',
              #vmax = Nrank,
              #alpha=0.5
              )


# Visualize current cell ownership on node
def plot_node(ax, n, lap, conf):
    tmp_grid = np.ones( (n.getNx(), n.getNy()) ) * -1.0
    
    #for i in range(n.getNx()):
    #    for j in range(n.getNy()):
    #        cid = n.cell_id(i,j)
    #        if n.is_local(cid):
    #            tmp_grid[i,j] = 0.5


    for cid in n.getCellIds():
        c = n.getCell( cid )
        (i, j) = c.index()
        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in real cells at ({},{})".format(n.rank, i,j))
            sys.exit()
        tmp_grid[i,j] = c.owner

    #XXX add back
    #for cid in n.getVirtuals():
    #    c = n.getCell( cid )
    #    (i,j) = c.index()
    #    if tmp_grid[i,j] != -1.0:
    #        print("{}: ERROR in virtual cells at ({},{})".format(n.rank, i,j))
    #        sys.exit()
    #    tmp_grid[i,j] = c.owner

    imshow(ax, tmp_grid, n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax() )


    # add text label about number of neighbors
    for cid in n.getCellIds():
        c = n.getCell( cid )
        (i, j) = c.index()
        dx = n.getXmax() - n.getXmin()
        dy = n.getYmax() - n.getYmin()

        ix = n.getXmin() + dx*(i+0.5)/n.getNx()
        jy = n.getYmin() + dy*(j+0.5)/n.getNy()

        #Nv = n.number_of_virtual_neighbors(c)
        Nv = c.number_of_virtual_neighbors
        label = str(Nv)
        #label = "{} ({},{})/{}".format(cid,i,j,Nv)
        #label = "({},{})".format(i,j)
        ax.text(ix, jy, label, ha='center',va='center', size=8)

    #for cid in n.getVirtuals():
    #    c = n.getCell( cid )
    #    (i,j) = c.index()
    #    ix = n.getXmin() + n.getXmax()*(i+0.5)/n.getNx()
    #    jy = n.getYmin() + n.getYmin()*(j+0.5)/n.getNy()
    #    label = "Vir"
    #    ax.text(jy, ix, label, ha='center',va='center')

    #XXX add back
    #ax.set_title(str(len(n.getVirtuals() ))+"/"+str(len(n.getCellIds() )))

    #save
    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/node_{}_{}.png'.format(n.rank, slap)
    plt.savefig(fname)






