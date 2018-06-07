try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
except:
    pass

import numpy as np


from visualize import plotNode
from visualize import plotTileBoundaries


class Particles:
    xs  = []
    ys  = []
    zs  = []

    uxs = []
    uys = []
    uzs = []

    def clear(self):
        self.xs  = []
        self.ys  = []
        self.zs  = []

        self.uxs = []
        self.uys = []
        self.uzs = []

def get_particles(node, conf):
    prtcl = Particles()
    prtcl.clear()

    for i in range(conf.Nx):
        for j in range(conf.Ny):
            for k in range(conf.Nz):
                cid = node.cellId(i,j)
                c = node.getCellPtr(cid)

                x, y, z, ux, uy, uz = get_particles_from_cell(c)

                prtcl.xs.extend(x)
                prtcl.ys.extend(y)
                prtcl.zs.extend(z)

                prtcl.uxs.extend(ux)
                prtcl.uys.extend(uy)
                prtcl.uzs.extend(uz)

    return prtcl


def get_particles_from_cell(cell):
    x  = cell.container.loc(0)
    y  = cell.container.loc(1)
    z  = cell.container.loc(2)

    ux = cell.container.vel(0)
    uy = cell.container.vel(1)
    uz = cell.container.vel(2)

    return x, y, z, ux, uy, uz



def plot2dParticles(ax, n, conf):

    #ax.clear()
    ax.cla()
    plotNode(ax, n, conf)
    plotTileBoundaries(ax, n, conf)

    prtcl = get_particles(n, conf)
    Np = len(prtcl.xs)
    #print("particles to plot: {}".format(Np))

    ax.plot(prtcl.xs, prtcl.ys, ".", color='black')
    











