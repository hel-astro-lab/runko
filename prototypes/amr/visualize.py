import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm
import numpy as np


def plot_center(ax, vb):
    x = vb.loc[0]
    y = vb.loc[1]
    z = vb.loc[2]
    ax.plot(x, y, marker='.', color='black')


def get_vertex(vb):
    x = vb.loc[0]
    y = vb.loc[1]
    z = vb.loc[2]
    dx = vb.dls[0]
    dy = vb.dls[1]
    dz = vb.dls[2]
   
    vrtx = np.zeros((2))
    vrty = np.zeros((2))
    vrtz = np.zeros((2))

    vrtx[0] = ( x - dx/2.0 )
    vrtx[1] = ( x + dx/2.0 )

    vrty[0] = ( y - dy/2.0 )
    vrty[1] = ( y + dy/2.0 )

    vrtz[0] = ( z - dz/2.0 )
    vrtz[1] = ( z + dz/2.0 )

    return vrtx, vrty, vrtz


def get_xy_bounding_curve(vb):
    vrtx, vrty, vrtz = get_vertex(vb)

    vrtxs = np.zeros(4)
    vrtys = np.zeros(4)

    vrtxs[0] = vrtx[0]
    vrtxs[1] = vrtx[1]
    vrtxs[2] = vrtx[1]
    vrtxs[3] = vrtx[0]

    vrtys[0] = vrty[0]
    vrtys[1] = vrty[0]
    vrtys[2] = vrty[1]
    vrtys[3] = vrty[1]

    return vrtxs, vrtys


def plot_edges(ax, vb):
    vrtxs, vrtys = get_xy_bounding_curve(vb)

    if vb.refLevel == 0:
        col = 'black'

    ax.plot(vrtxs, vrtys, linestyle='solid', color=col)



def visualize_mesh(ax, mesh):

    ax.cla()
    ax.minorticks_on()
    ax.set_xlim(-10.0, 10.0)
    ax.set_ylim(-10.0, 10.0)

    for vb in mesh.mesh:
        plot_center(ax, vb)

        plot_edges(ax, vb)



