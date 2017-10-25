import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

from visualize import *
import vmesh


#Profiling / Timing
sys.path.insert(0, '../tools')
from timer import Timer



#physical "real" distribution to compare against
def physical_vel(x,y,z, amp=1.0):

    mux = 5.0
    muy = 0.0
    muz = 0.0
    sigmax = 5.0
    sigmay = 6.0
    sigmaz = 4.0

    vx = np.exp(-(x-mux)**2 / sigmax**2 )
    vy = np.exp(-(y-muy)**2 / sigmay**2 )
    #vz = np.exp(-(z-muz)**2 / sigmaz**2 )
    vz = 1.0

    return amp*vx
    #return vx*vy*vz

def randab(a, b):
    return a + (b-a)*np.random.rand()


def populate_mesh( mesh, amp):

    for k in range(mesh.Nblocks[2]):
        for j in range(mesh.Nblocks[1]):
            for i in range(mesh.Nblocks[0]):
                cid = mesh.get_block_ID([i,j,k])
                (x,y,z) = mesh.get_center( cid )

                fval = physical_vel(x,y,z, amp)
                mesh[i,j,k] = [fval, fval, fval, fval]
                #print "({},{},{}) = {}".format(i,j,k,fval)



#struct holding all velocity mesh parameters
class velParams:
    mins = [ -20.0, -20.0, -20.0 ]
    maxs = [  20.0,  20.0,  20.0 ]
    lens = None

    Nx = 50
    Ny = 2
    Nz = 2


def createVelMesh(amp):
    ################################################## 
    # set-up grid
    params = velParams()

    mesh = vmesh.vMesh()
    mesh.Nblocks = [velParams.Nx, velParams.Ny, velParams.Nz] 
    mesh.zFill(params.mins, params.maxs)
    populate_mesh( mesh, amp)

    return mesh


#struct holding all node parameters
class NodeParams:
    Nx = 50
    Ny = 1


def imshow(ax, grid, xmin, xmax, ymin, ymax):
    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    extent = [xmin, xmax, ymin, ymax]

    mgrid = np.ma.masked_where(grid == 0.0, grid)

    ax.imshow(mgrid.T,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = 'plasma_r',
              vmin = 0.0,
              vmax = 1.0,
              aspect='auto',
              #alpha=0.5
              )


def visualizeNode(ax, n, nParams, vParams):
    data = np.zeros( (nParams.Nx, vParams.Nx) )

    #loop over all velocity grids
    for i in range(nParams.Nx):
        cell = n.get_cell_index(i, 0)
        i = cell.i
        j = cell.j

        #print i,j
        for vm in [cell.getData()]:
            #vbundle = vm.get_bundle(0, 2, 2) #xdir (dir = 0) @ j = 0, z = 0
            vbundle = vm.get_bundle(0, 0, 0) #xdir (dir = 0) @ j = 0, z = 0

            #print vbundle.getPencil()
            data[i, :] = vbundle.getPencil()

    #print data
    imshow(ax, data, 0.0, 1.0, vParams.mins[0], vParams.maxs[0])
    



if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    #axs.append( plt.subplot(gs[1]) )


    nParams = NodeParams()
    vParams = velParams()

    n = vmesh.Node()

    #load cells into local node
    for i in range(nParams.Nx):
        for j in range(nParams.Ny):

            #create Cell
            c = vmesh.Cell(i, j, 0)

            #create velocity mesh

            if i == 25:
                mesh = createVelMesh(1.0)
            else:
                mesh = createVelMesh(0.0)

            mesh.clip()
            c.addData(mesh)
            c.addData(mesh)

            #push to Node
            n.add_local_cell(c)


    #visualize full node content
    visualizeNode(axs[0], n, nParams, vParams)


    #save
    lap = 0
    stri = str(lap).rjust(4, '0')
    plt.savefig("out/sSolve_"+stri+".png")


    #now loop over solving the node content
    #vsol = vmesh.sSolver()
    #vsol.setNode(n)

    vsol = vmesh.sSolver(n)
    
    for lap in range(1,50):
        print "---lap: {}".format(lap)

        #solve each cell in the full 2D grid
        for i in range(nParams.Nx):
            for j in range(nParams.Ny):
                #print "({},{})".format(i,j)

                vsol.setTargetCell(i,j)
                vsol.solve()
        #vsol.update()
        n.cycle()
        

        visualizeNode(axs[0], n, nParams, vParams)
        stri = str(lap).rjust(4, '0')
        plt.savefig("out/sSolve_"+stri+".png")


def printNode(n):
    print

