import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import palettable as pal
cmap = pal.wesanderson.Moonrise1_5.mpl_colormap

import logi


Nrank = 4



##################################################
# plotting tools

# visualize matrix
def imshow(ax, grid):
    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(logi.xmin, logi.xmax)
    ax.set_ylim(logi.xmin, logi.xmax)

    extent = [logi.xmin, logi.xmax, logi.ymin, logi.ymax]

    mgrid = np.ma.masked_where(grid == -1.0, grid)
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
def plot_node(ax, n, lap):
    tmp_grid = np.ones( (logi.Nx, logi.Ny) ) * -1.0

    for cid in n.get_cells():
        c = n.get_cell( cid )

        (i, j) = c.index()

        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print "{}: ERROR in real cells at ({},{})".format(n.rank, i,j)
            sys.exit()
        tmp_grid[i,j] = c.owner

    print "virs:", n.get_virtuals()
    for cid in n.get_virtuals():
        print "virtual cid:",cid
        c = n.get_cell( cid )
        (i,j) = c.index()
        if tmp_grid[i,j] != -1.0:
            print "{}: ERROR in virtual cells at ({},{})".format(n.rank, i,j)
            sys.exit()
        tmp_grid[i,j] = c.owner

    imshow(ax, tmp_grid)


    # add text label about number of neighbors
    for cid in n.get_cells():
        c = n.get_cell( cid )
        (i, j) = c.index()
        ix = logi.xmin + logi.xmax*(i+0.5)/logi.Nx
        jy = logi.ymin + logi.ymax*(j+0.5)/logi.Ny

        #Nv = n.number_of_virtual_neighbors(c)
        Nv = c.number_of_virtual_neighbors
        label = str(Nv)
        #label = "({},{})/{}".format(i,j,Nv)
        #label = "({},{})".format(i,j)
        ax.text(jy, ix, label, ha='center',va='center', size=8)


    #for c in n.virtuals:
    #    (i, j) = c.index()
    #    ix = conf.xmin + conf.xmax*(i+0.5)/conf.Nx
    #    jy = conf.ymin + conf.ymax*(j+0.5)/conf.Ny
    #    label = "Vir"
    #    ax.text(jy, ix, label, ha='center',va='center')

    ax.set_title(str(len(n.get_virtuals() ))+"/"+str(len(n.get_cells() )))


    #save
    slap = str(lap).rjust(4, '0')
    fname = fpath + '/node_{}_{}.png'.format(n.rank, slap)
    plt.savefig(fname)






if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )


    c = logi.Cell( 0, 0, 0 )
    print "cid:", c.index()
    print "owner", c.owner
    print "neigs", c.neighs(-1, -1)
    print "nhood:", c.nhood()

    print "--------------------------------------------------"
    print "testing node..."

    n = logi.Node()
    n.init_mpi()

    # Path to be created (avoid clutter by issuing this with master)
    fpath = "out/"
    if n.master:
        if not os.path.exists(fpath):
            os.makedirs(fpath)



    #make random starting order
    np.random.seed(4)
    if n.master:
        for i in range(logi.Nx):
            for j in range(logi.Ny):
                val = np.random.randint(n.Nrank)
                n.set_mpiGrid(i, j, val)
    n.bcast_mpiGrid()


    #load cells into local node
    for i in range(logi.Nx):
        for j in range(logi.Ny):
            if n.mpiGrid(i,j) == n.rank:
                c = logi.Cell(i, j, n.rank)

                #TODO load data to cell
                n.add_cell(c)



    plot_node(axs[0], n, 0)

    n.pack_all_virtuals()

    plot_node(axs[0], n, 1)



    n.finalize_mpi()
