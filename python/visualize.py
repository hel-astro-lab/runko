import matplotlib.pyplot as plt
from matplotlib import cm
import palettable as pal
palette = pal.wesanderson.Moonrise1_5.mpl_colormap

import numpy as np


Nrank = 4


##################################################
# plotting tools

# visualize matrix
def imshow(ax, 
           grid, xmin, xmax, ymin, ymax,
           cmap='plasma',
           vmin = 0.0,
           vmax = 1.0,
           clip = -1.0,
          ):

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-3.0, 3.0)

    extent = [ xmin, xmax, ymin, ymax ]

    if clip == None:
        mgrid = grid
    else:
        mgrid = np.ma.masked_where(grid <= clip, grid)
    
    mgrid = mgrid.T
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = vmin,
              vmax = vmax,
              aspect='auto',
              #vmax = Nrank,
              #alpha=0.5
              )


# Visualize current cell ownership on node
def plotNode(ax, n, conf):
    tmp_grid = np.ones( (n.getNx(), n.getNy()) ) * -1.0
    
    #for i in range(n.getNx()):
    #    for j in range(n.getNy()):
    #        cid = n.cell_id(i,j)
    #        if n.is_local(cid):
    #            tmp_grid[i,j] = 0.5


    for cid in n.getCellIds():
        c = n.getCellPtr( cid )
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

    imshow(ax, tmp_grid, 
            n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
            cmap = palette,
            vmin = 0.0,
            vmax = Nrank-1
            )


    # add text label about number of neighbors
    for cid in n.getCellIds():
        c = n.getCellPtr( cid )
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



def saveVisz(lap, n, conf):

    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/node_{}_{}.png'.format(n.rank, slap)
    plt.savefig(fname)




# visualize vmesh content in x-dir
def plotXmesh(ax, n, conf, spcs):
    data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Nvx) )


    for i in range(conf.Nx):

        cid = n.cellId(i,0)
        c = n.getCellPtr(cid)

        #dig electron population out
        pgrid = c.getPlasmaGrid()

        for s in range(conf.NxMesh):

            if spcs == 0:
                vm = pgrid.electrons[s,0,0] #electron population
            elif spcs == 1:
                vm = pgrid.positrons[s,0,0] #positron population
            else:
                raise IndexError


            vbundle = vm.getBundle(0, 0, 0) #xdir (dir = 0) @ j = 0, z = 0

            data[ i*conf.NxMesh + s, :] = vbundle.getPencil()

    #data = np.log10(data)
    imshow(ax, data,
           n.getXmin(), n.getXmax(), conf.vxmin, conf.vxmax,
           cmap = 'plasma_r',
           #vmin = -10.0,
           #vmax =  2.0,
           vmin =   0.0,
           vmax =  10.0,
           clip =   0.0,
           )

    dmax = np.max(data)
    #print(dmax)

    #imshow(ax, data,
    #       n.getXmin(), n.getXmax(), conf.vxmin, conf.vxmax,
    #       cmap = 'plasma_r',
    #       vmin =  0.0,
    #       vmax = 10.0,
    #       clip =  0.0,
    #       )


def getYee(n, conf):

    data = {'x' : np.linspace(n.getXmin(), n.getXmax(), conf.Nx*conf.NxMesh),
            'ex': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ey': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ez': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ez': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'bx': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'by': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'bz': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jx': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jy': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jz': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'rh': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           }


    for i in range(conf.Nx):
        cid = n.cellId(i,0)
        c = n.getCellPtr(cid)

        yee = c.getYee()
        for s in range(conf.NxMesh):
            indx = i*conf.NxMesh + s

            data['ex'][indx] = yee.ex[s, 0, 0]
            data['ey'][indx] = yee.ey[s, 0, 0]
            data['ez'][indx] = yee.ez[s, 0, 0]

            data['bx'][indx] = yee.bx[s, 0, 0]
            data['by'][indx] = yee.by[s, 0, 0]
            data['bz'][indx] = yee.bz[s, 0, 0]

            data['jx'][indx] = yee.jx[s, 0, 0]
            data['jy'][indx] = yee.jy[s, 0, 0]
            data['jz'][indx] = yee.jz[s, 0, 0]

            data['rh'][indx] = yee.rh[s, 0, 0]

    return data



def plotJ(ax, n, conf):
    yee = getYee(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.getXmin(), n.getXmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)

    ax.plot(yee['x'], yee['jx'], "b-")
    ax.plot(yee['x'], yee['jy'], "r--")
    ax.plot(yee['x'], yee['jz'], "g--")
    

def plotE(ax, n, conf):
    yee = getYee(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.getXmin(), n.getXmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)

    ax.plot(yee['x'], yee['ex'], "b-")
    ax.plot(yee['x'], yee['ey'], "r--")
    ax.plot(yee['x'], yee['ez'], "g--")
    






