import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
except:
    pass

from visualize import imshow



class PyMesh:
    xx = None
    yy = None
    zz = None
    ff = None


def get_leaf_mesh(m, rfl_max):

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

    ff = np.zeros((nx, ny, nz))


    cells = m.get_cells(True)
    leafs = 0
    for cid in cells:
        if(not m.is_leaf(cid) ):
            continue
        leafs += 1

        rfl = m.get_refinement_level(cid)
        lvl = 2**rfl

        [i,j,k] = m.get_indices(cid)

        st = lvlm // lvl #stretch factor for ref lvl

        val = m[i,j,k, rfl]

        i *= st
        j *= st
        k *= st

        for ii in range(st):
            for jj in range(st):
                for kk in range(st):
                    ff[i+ii, j+jj, k+kk] = val

    #Nc = len(tiles)
    #print("tiles: {} / leafs {} (ratio: {}) / full {} (compression {})".format(
    #    Nc,
    #    leafs, 
    #    1.0*Nc/(1.0*leafs),
    #    nx*ny*nz, 
    #    Nc /(1.0*nx*ny*nz) 
    #    ))

    pym = PyMesh()
    pym.xx = xx
    pym.yy = yy
    pym.zz = zz
    pym.ff = ff

    return pym

#Naive (classical) sum of the distribution
def xMarginalize(pym):
    return np.sum(pym.ff, axis=(1,2) )

#slice along the middle
def getXmid(pym):
    return np.argmin( np.abs(pym.xx) )
def getYmid(pym):
    return np.argmin( np.abs(pym.yy) )
def getZmid(pym):
    return np.argmin( np.abs(pym.zz) )

def xSliceMid(pym):
    ymid = getYmid(pym)
    zmid = getZmid(pym)
    return pym.ff[:, ymid, zmid]

def ySliceMid(pym):
    xmid = getXmid(pym)
    zmid = getZmid(pym)
    return pym.ff[xmid, :, zmid]

def zSliceMid(pym):
    xmid = getXmid(pym)
    ymid = getYmid(pym)
    return pym.ff[xmid, ymid, :]


# visualize vmesh content in x-dir
def plotXmesh(ax, n, conf, spcs, vdir):

    if vdir == "x":
        fullNvx = np.int(conf.Nvx * (2.0**conf.refinement_level))
        fullNvx = fullNvx if conf.Nvx > 1 else 2*2**conf.refinement_level
        data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, fullNvx) )
    elif vdir == "y":
        fullNvy = np.int(conf.Nvy * (2.0**conf.refinement_level))
        fullNvy = fullNvy if conf.Nvy > 1 else 2*2**conf.refinement_level
        data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, fullNvy) )
    elif vdir == "z":
        fullNvz = np.int(conf.Nvz * (2.0**conf.refinement_level))
        fullNvz = fullNvz if conf.Nvz > 1 else 2*2**conf.refinement_level
        data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, fullNvz) )


    for i in range(conf.Nx):
        cid = n.id(i)
        c = n.get_tile(cid)

        block = c.get_plasma_species(0, spcs)
        for s in range(conf.NxMesh):
            vmesh = block[s,0,0]
            pym = get_leaf_mesh(vmesh, conf.refinement_level)
            
            if vdir == "x":
                sl = xSliceMid(pym) #slice from middle
                data[ i*conf.NxMesh + s, :] = sl

                dx = pym.xx[1]-pym.xx[0]
                vmin = pym.xx[0] -dx/2.0
                vmax = pym.xx[-1]+dx/1.0
            elif vdir == "y":
                sl = ySliceMid(pym) #slice from middle
                data[ i*conf.NxMesh + s, :] = sl

                dx = pym.yy[1]-pym.yy[0]
                vmin = pym.yy[0] -dx/2.0
                vmax = pym.yy[-1]+dx/1.0
            elif vdir == "z":
                sl = zSliceMid(pym) #slice from middle

                dx = pym.zz[1]-pym.zz[0]
                vmin = pym.zz[0] -dx/2.0
                vmax = pym.zz[-1]+dx/1.0

            data[ i*conf.NxMesh + s, :] = sl

    #print(np.max(data))
    data = data/data.max()

    imshow(ax, data,
           n.get_xmin(), n.get_xmax(), 
           vmin, vmax,
           cmap = 'plasma_r',
           vmin =   0.0,
           vmax =   1.0,
           clip =   None,
           )

    if vdir == "x":
        if spcs == 0:
            ax.set_ylabel(r'$v_{x,e}$')
        if spcs == 1:
            ax.set_ylabel(r'$v_{x,p}$')
    if vdir == "y":
        if spcs == 0:
            ax.set_ylabel(r'$v_{y,e}$')
        if spcs == 1:
            ax.set_ylabel(r'$v_{y,p}$')




def plot2DSlice(ax, vmesh, args):

    rfl_max = args["rfl"]
    mesh = get_leaf_mesh(vmesh, rfl_max)

    #normalizrfl_max
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



def plot2DSliceFlux(ax, vmesh, args):

    mesh = get_leaf_mesh(vmesh, args)

    #normalize
    fmax = np.max(np.abs(mesh.ff))
    mesh.ff = mesh.ff / fmax



    imshow(ax,
            mesh.ff,
            mesh.xx[0], mesh.xx[-1],
            mesh.yy[0], mesh.yy[-1],
            vmin =-1.0,
            vmax = 1.0,
            cmap = "RdBu",
            clip = None
            )
    return









