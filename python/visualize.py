try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import palettable as pal
    palette = pal.wesanderson.Moonrise1_5.mpl_colormap
except:
    pass

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
    im = ax.imshow(mgrid,
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
    return im



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

    ax.set_ylabel('node')


# plot tile boundaries
def plotTileBoundaries(ax, node, conf):

    for i in range(conf.Nx):
        for j in range(conf.Ny):
            for k in range(conf.Nz):
                cid = node.cellId(i,j)
                c = node.getCellPtr(cid)

                mins = np.array( c.mins ) 
                maxs = np.array( c.maxs )
                #lens = np.array( [conf.NxMesh+1, conf.NyMesh+1, conf.NzMesh+1] )
                #ds = (maxs - mins)/lens
                xx = np.linspace(mins[0], maxs[0], conf.NxMesh+1)
                yy = np.linspace(mins[1], maxs[1], conf.NyMesh+1)
                zz = np.linspace(mins[2], maxs[2], conf.NzMesh+1)
                
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
        cid = n.cellId(i,0)
        c = n.getCellPtr(cid)

        yee = c.getYee(0)
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

            data['jx1'][indx] = yee.jx1[s, 0, 0]

            data['rho'][indx] = yee.rho[s, 0, 0]

    return data


def getYee2D(n, conf):

    data = {'x' : np.linspace(n.getXmin(), n.getXmax(), conf.Nx*conf.NxMesh),
            'y' : np.linspace(n.getYmin(), n.getYmax(), conf.Ny*conf.NyMesh),
            'ex':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'ey':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'ez':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'ez':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'bx':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'by':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'bz':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jx':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jy':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jz':   -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'jx1':  -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
            'rho':  -1.0 * np.ones( (conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh) ),
           }

    for i in range(conf.Nx):
        for j in range(conf.Ny):
            cid = n.cellId(i,j)
            c = n.getCellPtr(cid)

            yee = c.getYee(0)
            for r in range(conf.NyMesh):
                for q in range(conf.NxMesh):

                    indx = i*conf.NxMesh + q
                    jndx = j*conf.NyMesh + r

                    data['ex'][indx, jndx] = yee.ex[q, r, 0]
                    data['ey'][indx, jndx] = yee.ey[q, r, 0]
                    data['ez'][indx, jndx] = yee.ez[q, r, 0]

                    data['bx'][indx, jndx] = yee.bx[q, r, 0]
                    data['by'][indx, jndx] = yee.by[q, r, 0]
                    data['bz'][indx, jndx] = yee.bz[q, r, 0]
                                                        
                    data['jx'][indx, jndx] = yee.jx[q, r, 0]
                    data['jy'][indx, jndx] = yee.jy[q, r, 0]
                    data['jz'][indx, jndx] = yee.jz[q, r, 0]

                    data['jx1'][indx, jndx] = yee.jx1[q, r, 0]

                    data['rho'][indx, jndx] = yee.rho[q, r, 0]

    return data


# species-specific analysis meshes (plasma moments)
def getAnalysis(n, conf, ispcs):

    data = {'x' : np.linspace(n.getXmin(), n.getXmax(), conf.Nx*conf.NxMesh),
           'rho':    -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'mgamma': -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Vx':     -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Vy':     -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Vz':     -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Tx':     -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Ty':     -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'Tz':     -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           'ekin':   -1.0 * np.ones( (conf.Nx*conf.NxMesh) ),
           }

    for i in range(conf.Nx):
        cid = n.cellId(i,0)
        c = n.getCellPtr(cid)

        analysis = c.getAnalysis(ispcs)
        for s in range(conf.NxMesh):
            indx = i*conf.NxMesh + s

            data['rho'][indx] = analysis.rho[s, 0, 0]

            data['mgamma'][indx] = analysis.mgamma[s, 0, 0]

            data['Vx'][indx] = analysis.Vx[s, 0, 0]
            data['Vy'][indx] = analysis.Vy[s, 0, 0]
            data['Vz'][indx] = analysis.Vz[s, 0, 0]

            data['Tx'][indx] = analysis.Tx[s, 0, 0]
            data['Ty'][indx] = analysis.Ty[s, 0, 0]
            data['Tz'][indx] = analysis.Tz[s, 0, 0]

            data['ekin'][indx] = analysis.ekin[s, 0, 0]

    return data

def plotJ(ax, n, conf):
    yee = getYee(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.getXmin(), n.getXmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)

    ax.plot(yee['x'], yee['jx'], "b-")
    #ax.plot(yee['x'], yee['jy'], "r--")
    #ax.plot(yee['x'], yee['jz'], "g--")

    ax.plot(yee['x'], yee['jx1'], "r--")

    #ratio
    #ax.plot(yee['x'], yee['jx1']/yee['jx'], "k-")

    
    ax.set_ylabel(r'$J_x$')


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
    
    ax.set_ylabel(r'$E_x$')


def plotDens(ax, n, conf):
    yee = getYee(n, conf)

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(n.getXmin(), n.getXmax())
    #ax.set_xlim(-3.0, 3.0)
    #ax.set_ylim(-20.0, 20.0)

    ax.plot(yee['x'], yee['rho'], "b-")
    
    ax.set_ylabel(r'$n$')



def plot2dYee(ax, n, conf, val = 'jx'):

    #ax.clear()
    ax.cla()
    yee = getYee2D(n, conf)

    vmin, vmax = np.min(yee[val]), np.max(yee[val])
    vminmax = np.maximum( np.abs(vmin), np.abs(vmax) )
    print("2D {} min{} max {} minmax {}".format(val, vmin, vmax, vminmax))

    imshow(ax, yee[val],
           n.getXmin(), n.getXmax(), n.getYmin(), n.getYmax(),
           cmap = "RdBu",
           vmin = -vminmax,
           vmax =  vminmax,
           clip = 0.0
          )
    ax.set_title(val)





