try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
except:
    pass

import numpy as np
import sys


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
           cap = None,
           aspect = 'auto',
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
    elif type(clip) == tuple:
        cmin, cmax = clip
        print(cmin, cmax)
        mgrid = np.ma.masked_where( np.logical_and(cmin <= grid, grid <= cmax), grid)
    else:
        mgrid = np.ma.masked_where(grid <= clip, grid)

    if cap != None:
        mgrid = np.clip(mgrid, cap )
    
    mgrid = mgrid.T
    im = ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = vmin,
              vmax = vmax,
              aspect=aspect,
              #vmax = Nrank,
              #alpha=0.5
              )
    return im


# Visualize current cell ownership on grid
def plotNode(ax, n, conf):
    tmp_grid = np.ones( (n.get_Nx(), n.get_Ny()) ) * -1.0

    for cid in n.get_tile_ids():
        c = n.get_tile( cid )
        (i, j) = c.index
        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print("{}: ERROR in real cells at ({},{})".format(n.rank, i,j))
            sys.exit()
        tmp_grid[i,j] = c.communication.owner

    imshow(ax, tmp_grid, 
            n.get_xmin(), n.get_xmax(), n.get_ymin(), n.get_ymax(),
            cmap = palette,
            vmin = 0.0,
            vmax = n.size(),
            )

    # add text label about number of neighbors
    for cid in n.get_tile_ids():
        c = n.get_tile( cid )
        (i, j) = c.index
        dx = n.get_xmax() - n.get_xmin()
        dy = n.get_ymax() - n.get_ymin()

        ix = n.get_xmin() + dx*(i+0.5)/n.get_Nx()
        jy = n.get_ymin() + dy*(j+0.5)/n.get_Ny()

        Nv = c.communication.number_of_virtual_neighbors
        label = str(Nv)
        ax.text(ix, jy, label, ha='center',va='center', size=8)

    #mark boundaries with hatch
    dx = n.get_xmax() - n.get_xmin()
    dy = n.get_ymax() - n.get_ymin()
    for cid in n.get_boundary_tiles():
        c = n.get_tile( cid )
        (i, j) = c.index

        ix0 = n.get_xmin() + dx*(i+0.0)/n.get_Nx()
        jy0 = n.get_ymin() + dy*(j+0.0)/n.get_Ny()

        ix1 = n.get_xmin() + dx*(i+1.0)/n.get_Nx()
        jy1 = n.get_ymin() + dy*(j+1.0)/n.get_Ny()

        #ax.fill_between([ix0,ix1], [jy0, jy1], hatch='///', alpha=0.0)

        ax.plot([ix0, ix0],[jy0, jy1], color='k', linestyle='dotted')
        ax.plot([ix1, ix1],[jy0, jy1], color='k', linestyle='dotted')
        ax.plot([ix0, ix1],[jy0, jy0], color='k', linestyle='dotted')
        ax.plot([ix0, ix1],[jy1, jy1], color='k', linestyle='dotted')






# plot tile boundaries
def plotTileBoundaries(ax, grid, conf):

    for i in range(conf.Nx):
        for j in range(conf.Ny):
            for k in range(conf.Nz):
                try:
                    cid = grid.id(i,j)
                except:
                    cid = grid.id(i)

                c = grid.get_tile(cid)

                mins = np.array( c.mins ) 
                maxs = np.array( c.maxs )
                #lens = np.array( [conf.NxMesh+1, conf.NyMesh+1, conf.NzMesh+1] )
                #ds = (maxs - mins)/lens
                xx = np.linspace(mins[0], maxs[0], conf.NxMesh+1)
                yy = np.linspace(mins[1], maxs[1], conf.NyMesh+1)
                #zz = np.linspace(mins[2], maxs[2], conf.NzMesh+1)
                
                # inner mesh
                for y in yy:
                    ax.plot( [xx[0], xx[-1]], [y, y], "k", linestyle='dotted')
                for x in xx:
                    ax.plot( [x, x], [yy[0], yy[-1]], "k", linestyle='dotted')

                #and outer boundaries
                ax.plot([mins[0], maxs[0]], [mins[1], mins[1]], "k-") #bottom
                ax.plot([mins[0], maxs[0]], [maxs[1], maxs[1]], "k-") #top
                ax.plot([mins[0], mins[0]], [mins[1], maxs[1]], "k-") #left
                ax.plot([maxs[0], maxs[0]], [mins[1], maxs[1]], "k-") #right


def saveVisz(lap, n, conf):

    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/node_{}_{}.png'.format(n.rank(), slap)
    plt.savefig(fname)




# visualize vmesh content in x-dir
def plotXmesh(ax, n, conf, spcs):
    data = -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Nvx) )


    for i in range(conf.Nx):

        cid = n.id(i)
        c = n.get_tile(cid)

        #dig electron population out
        pgrid = c.getPlasmaGrid()

        for s in range(conf.NxMesh):

            if spcs == 0:
                vm = pgrid.electrons[s,0,0] #electron population
            elif spcs == 1:
                vm = pgrid.positrons[s,0,0] #positron population
            else:
                raise IndexError


            vbundle = vm.get_bundle(0, 0, 0) #xdir (dir = 0) @ j = 0, z = 0

            data[ i*conf.NxMesh + s, :] = vbundle.get_pencil()

    #data = np.log10(data)
    imshow(ax, data,
           n.get_xmin(), n.get_xmax(), conf.vxmin, conf.vxmax,
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
    #       n.get_xmin(), n.get_xmax(), conf.vxmin, conf.vxmax,
    #       cmap = 'plasma_r',
    #       vmin =  0.0,
    #       vmax = 10.0,
    #       clip =  0.0,
    #       )


def get_grids(n, conf):

    data = {'x' : np.linspace(n.get_xmin(), n.get_xmax(), conf.Nx*conf.NxMesh),
            'ex':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ey':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ez':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'ez':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'bx':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'by':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'bz':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jx':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jy':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jz':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'jx1':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
            'rho':  -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           }

    for i in range(conf.Nx):
        cid = n.id(i)
        c = n.get_tile(cid)

        gs = c.get_grids(0)
        for s in range(conf.NxMesh):
            indx = i*conf.NxMesh + s

            data['ex'][indx] = gs.ex[s, 0, 0]
            data['ey'][indx] = gs.ey[s, 0, 0]
            data['ez'][indx] = gs.ez[s, 0, 0]

            data['bx'][indx] = gs.bx[s, 0, 0]
            data['by'][indx] = gs.by[s, 0, 0]
            data['bz'][indx] = gs.bz[s, 0, 0]

            data['jx'][indx] = gs.jx[s, 0, 0]
            data['jy'][indx] = gs.jy[s, 0, 0]
            data['jz'][indx] = gs.jz[s, 0, 0]

            #data['jx1'][indx] = gs.jx1[s, 0, 0]

            data['rho'][indx] = gs.rho[s, 0, 0]

    return data


def getYee2D(n, conf):

    data = {'x' : np.linspace(n.get_xmin(), n.get_xmax(), conf.Nx*conf.NxMesh),
            'y' : np.linspace(n.get_ymin(), n.get_ymax(), conf.Ny*conf.NyMesh),
            'ex':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'ey':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'ez':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'ez':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'bx':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'by':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'bz':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jx':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jy':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jz':   -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jx1':  -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jy1':  -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jz1':  -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'rho':  -np.inf * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
           }

    #for cid in n.get_local_tiles():
    for cid in n.get_tile_ids():
        c = n.get_tile(cid)

        try:
            i,j = c.index
        except:
            i,j,k = c.index

        gs = c.get_grids(0)
        for r in range(conf.NyMesh):
            for q in range(conf.NxMesh):
                indx = i*conf.NxMesh + q
                jndx = j*conf.NyMesh + r

                data['ex'][indx, jndx] = gs.ex[q, r, 0]
                data['ey'][indx, jndx] = gs.ey[q, r, 0]
                data['ez'][indx, jndx] = gs.ez[q, r, 0]

                data['bx'][indx, jndx] = gs.bx[q, r, 0]
                data['by'][indx, jndx] = gs.by[q, r, 0]
                data['bz'][indx, jndx] = gs.bz[q, r, 0]
                                                    
                data['jx'][indx, jndx] = gs.jx[q, r, 0]
                data['jy'][indx, jndx] = gs.jy[q, r, 0]
                data['jz'][indx, jndx] = gs.jz[q, r, 0]

                #data['jx1'][indx, jndx] = gs.jx1[q, r, 0]
                #data['jy1'][indx, jndx] = gs.jy1[q, r, 0]
                #data['jz1'][indx, jndx] = gs.jz1[q, r, 0]

                data['rho'][indx, jndx] = gs.rho[q, r, 0]

    return data


# species-specific analysis meshes (plasma moments)
def get_analysis(n, conf, ispcs):

    data = {'x' : np.linspace(n.get_xmin(), n.get_xmax(), conf.Nx*conf.NxMesh),
           'rho':     -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'edens':   -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'temp':    -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'Vx':      -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'Vy':      -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'Vz':      -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'momx':    -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'momy':    -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'momz':    -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'pressx':  -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'pressy':  -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'pressz':  -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'shearxy': -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'shearxz': -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           'shearyz': -np.inf* np.ones( (conf.Nx*conf.NxMesh) ),
           }

    for i in range(conf.Nx):
        cid = n.id(i)
        c = n.get_tile(cid)

        analysis = c.get_analysis(ispcs)
        for s in range(conf.NxMesh):
            indx = i*conf.NxMesh + s

            data['rho'][indx] = analysis.rho[s, 0, 0]
            data['edens'][indx] = analysis.edens[s, 0, 0]
            data['temp'][indx] = analysis.temp[s, 0, 0]

            data['Vx'][indx] = analysis.Vx[s, 0, 0]
            data['Vy'][indx] = analysis.Vy[s, 0, 0]
            data['Vz'][indx] = analysis.Vz[s, 0, 0]

            data['momx'][indx] = analysis.momx[s, 0, 0]
            data['momy'][indx] = analysis.momy[s, 0, 0]
            data['momz'][indx] = analysis.momz[s, 0, 0]

            data['pressx'][indx] = analysis.pressx[s, 0, 0]
            data['pressy'][indx] = analysis.pressy[s, 0, 0]
            data['pressz'][indx] = analysis.pressz[s, 0, 0]

            data['shearxy'][indx] = analysis.shearxy[s, 0, 0]
            data['shearxz'][indx] = analysis.shearxz[s, 0, 0]
            data['shearyz'][indx] = analysis.shearyz[s, 0, 0]

    return data

def plotJ(ax, n, conf):
    gs = get_grids(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.get_xmin(), n.get_xmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)

    ax.plot(gs['x'], gs['jx'], "b-")
    #ax.plot(gs['x'], gs['jy'], "r--")
    #ax.plot(gs['x'], gs['jz'], "g--")

    ax.plot(gs['x'], gs['jx1'], "r--")

    #ratio
    #ax.plot(gs['x'], gs['jx1']/gs['jx'], "k-")

    
    ax.set_ylabel(r'$J_x$')


def plotE(ax, n, conf):
    gs = get_grids(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.get_xmin(), n.get_xmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)

    ax.plot(gs['x'], gs['ex'], "b-")
    ax.plot(gs['x'], gs['ey'], "r--")
    ax.plot(gs['x'], gs['ez'], "g--")
    
    ax.set_ylabel(r'$E_x$')


def plotDens(ax, n, conf):
    gs = get_grids(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.get_xmin(), n.get_xmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)
    #ax.set_yscale('log')

    ax.plot(gs['x'], gs['rho'], "b-")
    
    ax.set_ylabel(r'$n$')



def plot2dYee(
        ax, 
        gs, 
        n, 
        conf, 
        val = 'jx', 
        vmin=None,
        vmax=None,
        label_title=False,
        ):

    #ax.clear()
    ax.cla()
    #gs = getYee2D(n, conf)

    #filter non-existent values away
    arr = gs[val]
    arr = np.ma.masked_where(arr == -np.inf, arr)

    if vmin == None:
        vmin = np.min(arr)
    if vmax == None:
        vmax = np.max(arr)

    vminmax = np.maximum( np.abs(vmin), np.abs(vmax) )
    #print("2D {} min{} max {} minmax {}".format(val, vmin, vmax, vminmax))


    imshow(ax, arr,
           n.get_xmin(), n.get_xmax(), n.get_ymin(), n.get_ymax(),
           cmap = "RdBu",
           vmin = -vminmax,
           vmax =  vminmax,
           clip = None,
          )

    if label_title:
        ax.set_ylabel(val)
    else:
        ax.set_title(val)





