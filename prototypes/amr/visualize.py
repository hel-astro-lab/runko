import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm
import numpy as np

import vmesh


def plot_center(ax, mesh, cid):
    (x,y,z) = mesh.get_center(cid)
    ax.plot(x, y, marker='.', color='black')


def get_vertex(mesh, cid):
    (x,y,z) = mesh.get_center(cid)
    (dx, dy, dz) = mesh.get_size(cid)
   
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


def get_xy_bounding_curve(mesh, cid):
    vrtx, vrty, vrtz = get_vertex(mesh, cid)

    vrtxs = np.zeros(5)
    vrtys = np.zeros(5)

    vrtxs[0] = vrtx[0]
    vrtxs[1] = vrtx[1]
    vrtxs[2] = vrtx[1]
    vrtxs[3] = vrtx[0]
    vrtxs[4] = vrtx[0]

    vrtys[0] = vrty[0]
    vrtys[1] = vrty[0]
    vrtys[2] = vrty[1]
    vrtys[3] = vrty[1]
    vrtys[4] = vrty[0]

    return vrtxs, vrtys


def plot_edges(ax, mesh, cid, alpha=0.3):
    vrtxs, vrtys = get_xy_bounding_curve(mesh, cid)

    #if vb.refLevel == 0:
    #    col = 'black'
    col = 'black'

    ax.plot(vrtxs, vrtys, linestyle='solid', color=col, alpha=alpha)



def visualize_mesh(ax, mesh):

    ax.cla()
    ax.minorticks_on()
    ax.set_xlim(-11.0, 11.0)
    ax.set_ylim(-11.0, 11.0)

    for cellID in mesh.all_blocks(True):
        plot_center(ax, mesh, cellID)
        plot_edges(ax, mesh, cellID)



def visualize_data(ax, mesh):
    ax.cla()
    ax.minorticks_on()
    ax.set_xlim(-11.0, 11.0)
    ax.set_ylim(-11.0, 11.0)

    (Nx, Ny, Nz) = mesh.nCells
    data = np.zeros( (Nx, Ny) ) #xy
    #data = np.zeros( (Nz, Ny) ) #zy


    for cid in mesh.all_blocks(True):
        (i,j,k) = mesh.get_indices( cid )
        val = mesh[cid]

        data[i,j] = val[0] #xy
        #data[j,k] = val[0] #zy

        #print "({},{},{}) = {}".format(i,j,k,val[0])
        #plot_edges(ax, mesh, cid)

    extent = [-10.0, 10.0, -10.0, 10.0]
    mgrid = np.ma.masked_where(data == 0.0, data)
    ax.imshow(mgrid.T,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cm.get_cmap('plasma'),
              vmin = 0.0,
              vmax = 1.0,
              aspect='auto',
              )

