import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

from visualize import imshow



class PyMesh:
    xx = None
    yy = None
    ff = None


def get_leaf_mesh(m, args):

    rfl_max = args["rfl"]
    nx, ny, nz = m.get_size(rfl_max)
    lvlm = 2**rfl_max

    #empty arrays ready
    xx = np.zeros((nx))
    yy = np.zeros((ny))
    zz = np.zeros((nz))

    for i in range(nx):
        x,y,z = m.get_center([i,0,0], rfl_max)
        xx[i] = x
    for j in range(ny):
        x,y,z = m.get_center([0,j,0], rfl_max)
        yy[j] = y
    for k in range(nz):
        x,y,z = m.get_center([0,0,k], rfl_max)
        zz[k] = z


    if args["q"] == "mid":
        if args["dir"] == "xy":
            q = nz/2
        elif args["dir"] == "xz":
            q = ny/2
        elif args["dir"] == "yz":
            q = nx/2
    else:
        q = args["q"]

    if args["dir"] == "xy":
        ff = np.zeros((nx, ny))
    elif args["dir"] == "xz":
        ff = np.zeros((nx, nz))
    elif args["dir"] == "yz":
        ff = np.zeros((ny, nz))


    cells = m.get_cells(True)
    leafs = 0
    for cid in cells:

        if(not m.is_leaf(cid) ):
            continue

        leafs += 1

        rfl = m.get_refinement_level(cid)
        lvl = 2**rfl

        [i,j,k] = m.get_indices(cid)

        st = lvlm/lvl #stretch factor for ref lvl

        if args["dir"] == "xy" and k*st != q:
            continue
        if args["dir"] == "xz" and j*st != q:
            continue
        if args["dir"] == "yz" and i*st != q:
            continue

        val = m[i,j,k, rfl]

        i *= st
        j *= st
        k *= st

        if args["dir"] == "xy":
            for ii in range(st):
                for jj in range(st):
                    ff[i+ii, j+jj] = val
        if args["dir"] == "xz":
            for ii in range(st):
                for jj in range(st):
                    ff[i+ii, k+jj] = val
        if args["dir"] == "yz":
            for ii in range(st):
                for jj in range(st):
                    ff[j+ii, k+jj] = val


    Nc = len(cells)
    print("cells: {} / leafs {} (ratio: {}) / full {} (compression {})".format(
        Nc,
        leafs, 
        Nc/leafs,
        nx*ny*nz, 
        Nc /(nx*ny*nz) 
        ))


    m = PyMesh()
    m.xx = xx
    m.yy = yy
    m.ff = ff

    return m



def plot2DSlice(ax, vmesh, args):

    #mesh = get_mesh(m, args)
    #mesh = get_interpolated_mesh(m, args)
    mesh = get_leaf_mesh(vmesh, args)

    #normalize
    mesh.ff = mesh.ff / np.max(mesh.ff)

    imshow(ax,
            mesh.ff,
            mesh.xx[0], mesh.xx[-1],
            mesh.yy[0], mesh.yy[-1],
            vmin = 0.0,
            vmax = 1.0,
            cmap = "plasma_r",
            clip = 0.0
            )
    return









