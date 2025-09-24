try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
except:
    pass

import numpy as np
import random


from visualize import plotNode
from visualize import plotTileBoundaries


class Particles:
    xs  = []
    ys  = []
    zs  = []

    uxs = []
    uys = []
    uzs = []

    wgt = []

    def clear(self):
        self.xs  = []
        self.ys  = []
        self.zs  = []

        self.uxs = []
        self.uys = []
        self.uzs = []

        self.wgt = []


def get_particles(grid, conf, ip):
    prtcl = Particles()
    prtcl.clear()

    for i in range(conf.Nx):
        for j in range(conf.Ny):
            for k in range(conf.Nz):
                if grid.get_mpi_grid(i,j) == grid.rank():
                    cid = grid.id(i,j)
                    c = grid.get_tile(cid)

                    x, y, z, ux, uy, uz, wgt = get_particles_from_tile(c, ip)

                    prtcl.xs.extend(x)
                    prtcl.ys.extend(y)
                    prtcl.zs.extend(z)

                    prtcl.uxs.extend(ux)
                    prtcl.uys.extend(uy)
                    prtcl.uzs.extend(uz)

                    prtcl.wgt.extend(wgt)

    return prtcl


def get_particles_from_tile(tile, ispcs):
    container = tile.get_container(ispcs)
    x  = container.loc(0)
    y  = container.loc(1)
    z  = container.loc(2)

    ux = container.vel(0)
    uy = container.vel(1)
    uz = container.vel(2)

    wgt = container.wgt()

    return x, y, z, ux, uy, uz, wgt



def plot2dParticles(ax, n, conf, downsample=0):

    #ax.clear()
    #ax.cla()
    plotNode(ax, n, conf)
    #plotTileBoundaries(ax, n, conf)

    prtcl = get_particles(n, conf, 0)
    Np = len(prtcl.xs)
    #print("particles to plot: {}".format(Np))

    if downsample > 0:
        rindxs = random.sample( range(0, Np-1), int(downsample*Np) )

        prtcl.xs = np.array( prtcl.xs )
        prtcl.ys = np.array( prtcl.ys ) 
        prtcl.zs = np.array( prtcl.zs )

        prtcl.xs = prtcl.xs[rindxs]
        prtcl.ys = prtcl.ys[rindxs]
        prtcl.zs = prtcl.zs[rindxs]

    ax.plot(prtcl.xs, prtcl.ys, ".", color='red')
    











