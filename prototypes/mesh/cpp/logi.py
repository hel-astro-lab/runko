import numpy as np
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

from logi import Cell
from logi import Node









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


    c = Cell( 0, 0, 0 )
    print "cid:", c.index()
    print "owner", c.owner
    print "neigs", c.neighs(-1, -1)
    print "nhood:", c.nhood()


    n = Node()

    n.initMPI()



    n.finalizeMPI()
