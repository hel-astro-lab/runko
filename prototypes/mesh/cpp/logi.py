import numpy as np
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

import logi









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








    n.finalize_mpi()
